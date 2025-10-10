#!/usr/bin/env nextflow
process vcf_to_bim {
    tag "${vcf_file.baseName}"
    publishDir "data/interim/bim_files", mode: 'copy'

    input:
    tuple path(vcf_file), val(bim_file)

    output:
    path "${bim_file}.bim"

    script:
    """
    plink --vcf ${vcf_file} --make-just-bim --out ${bim_file}
    """
}

process clean_phs {
    tag "${bim_file.baseName}"
    // debug true

    input:
    val cluster_name
    path script_file
    path phs_file
    path bim_file
    val build

    output:
    path "phs601_${cluster_name}_${bim_file.baseName}.csv"

    script:
    """
    python3 ${script_file} \
        --phs_file ${phs_file} \
        --bim_file ${bim_file} \
        --build ${build} \
        --out phs601_${cluster_name}_${bim_file.baseName}.csv
    """
}

process merge_phs_outputs {
    tag "merge_$cluster_name"
    publishDir "data/interim", mode: 'copy'

    input:
    val cluster_name
    path csv_files

    output:
    path "phs601_${cluster_name}.csv"

    script:
    """
    awk 'FNR==1 && NR!=1{next} {print}' ${csv_files.join(' ')} > phs601_${cluster_name}.csv
    rm -f ${csv_files.join(' ')}
    """
}

process preprocess_metadata {
    publishDir "data/processed", mode = 'copy'

    input:
    val cluster_name
    path preprocess_script
    path 'src/config.R'
    path 'src/functions.R'

    output:
    path "${cluster_name}_phenotypes.tsv"

    script:
    """
    Rscript ${preprocess_script}
    """
}

workflow harmonize_phs {
    take:
    cluster_name
    bim_files
    
    main:
    script_file = file('src/PHS_PC.py')
    phs_file = file(params.phs_file)
    build = params.build

    harmonized_chr = clean_phs(cluster_name, script_file, phs_file, bim_files, build)

    harmonized_chr
        .collect()
        .set { all_csvs }

    merge_phs_outputs(cluster_name, all_csvs)
}

workflow {
    cluster_name = params.cluster_name

    if ( cluster_name == 'p2521' ) {
        if ( file('data/interim/bim_files/husk_all.imputed.chr1.bim').exists() ) {
            println "Skipping VCF conversion, BIM files already present"
            bim_files = Channel.fromPath('data/interim/bim_files/*.bim')
        } else {
        // Create channel of VCF files
        vcf_files = Channel.fromPath(params.vcf_file)

        // Map each VCF to its corresponding BIM filename
        bim_info = vcf_files.map { vcf_file ->
            def bim_file = vcf_file.name.replaceFirst(/\.vcf(\.gz)?$/, '')  // handles .vcf or .vcf.gz
            tuple(vcf_file, bim_file)
        }

        // Call the vcf_to_bim process
        bim_files = vcf_to_bim(bim_info)
        }
    } else {
        bim_files = Channel.fromPath(params.bim_template)
    }

    harmonize_phs(cluster_name, bim_files)

    if ( cluster_name != 'local' ) {
        preprocess_metadata(cluster_name,
  			    file(params.preprocess_script),
  			    file(params.r_config),
    			    file(params.r_functions))
    }
}
