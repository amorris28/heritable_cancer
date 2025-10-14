#!/usr/bin/env nextflow

process vcf_to_bim {
    tag "${vcf_file.baseName}"
    storeDir "${params.interim_dir}/vcf_to_bim"

    input:
    path vcf_file

    output:
    path "${vcf_file.getBaseName(2)}.bim"

    script:
    """
    plink \
        --vcf ${vcf_file} \
        --make-just-bim \
        --out ${vcf_file.getBaseName(2)}
    """
}

process harmonize_phs_by_chr {
    tag "${bim_file.baseName}"
    container params.python3_image

    input:
    path script_file
    path phs_file
    path bim_file
    val build

    output:
    path "phs601_${workflow.profile}_${bim_file.baseName}.csv", emit: harmonized_phs

    script:
    """
    python3 ${script_file} \
        --phs_file ${phs_file} \
        --bim_file ${bim_file} \
        --build ${build} \
        --out phs601_${workflow.profile}_${bim_file.baseName}.csv
    """
}

process merge_phs_outputs {
    tag "merge_$workflow.profile"
    publishDir params.processed_dir, mode: 'copy'

    input:
    path csv_files

    output:
    path "phs601_${workflow.profile}.csv", emit: csv
    path "phs601_${workflow.profile}.snps", emit: snp

    script:
    """
    awk 'FNR==1 && NR!=1{next} {print}' ${csv_files.join(' ')} > phs601_${workflow.profile}.csv
    cat phs601_${workflow.profile}.csv | awk '{print \$6}' | tail -n +1 > phs601_${workflow.profile}.snps
    """
}

process filter_bfiles {
    storeDir "${params.interim_dir}/extract_snps"

    input:
    path bfiles
    path snp_file

    output:
    path "${bfiles.getBaseName(2)}.phs601.{bed,bim,fam}", emit: bfiles
    
    script:
    """
    plink2 \
        --bfile ${input_file} \
        ${psam_flag} \
        --extract ${snp_file} \
        --make-bed \
        --out ${bfiles.getBaseName(2)}.phs601
    """
}

process filter_vcfs {
    storeDir "${params.interim_dir}/extract_snps"

    input:
    path vcf_file
    path snp_file

    output:
    path "${vcf_file.getBaseName(2)}.phs601.{bed,bim,fam}", emit: bfiles
    
    script:
    """
    plink2 \
        --vcf ${vcf_file} \
        --extract ${snp_file} \
        --make-bed \
        --out ${vcf_file.getBaseName(2)}.phs601 \
        --memory 4000
    """
}
process merge_bfiles {
    container params.gwas_image

    input:
    path merge_list
    path bfiles
     
    output:
    path "${workflow.profile}_phs601.{bed,bim,fam}", emit: bfiles

    script:
    """
    plink \
        --merge-list ${merge_list} \
        --make-bed \
        --out ${workflow.profile}_phs601
    """
}
process plink_score {
    publishDir params.processed_dir, mode: 'copy' 
    container params.gwas_image

    input:
    tuple val(base), path(bfiles)
    path phs_csv
     
    output:
    path "${workflow.profile}_phs601.profile"
    
    script:
    """
    plink \
        --bfile ${base} \
        --score ${phs_csv} 6 7 5 header sum \
        --out ${workflow.profile}_phs601
    """
}
process preprocess_metadata {
    publishDir params.processed_dir, mode: 'copy'
    container params.r_image

    input:
    path preprocess_script
    path 'src/config.R'
    path 'src/functions.R'

    output:
    path "${workflow.profile}_phenotypes.tsv"

    script:
    """
    mkdir -p data/processed
    Rscript \
        ${preprocess_script} \
        ${workflow.profile}_phenotypes.tsv
    """
}
workflow get_bim_files {
    take:
    vcf_files

    main:
    bim = params.bim_template
        ? Channel.fromPath(params.bim_template, checkIfExists:true)
        : vcf_to_bim(vcf_files)

    emit:
    bim_files = bim
}

workflow harmonize_models {
    take:
    bim_files

    main:
    harmonized_chr = harmonize_phs_by_chr(file(params.harmonize_script), 
                                          file(params.phs_file), 
                                          bim_files, 
                                          params.build)
    
    merge_phs_outputs(harmonized_chr.collect())

    emit:
    csv = merge_phs_outputs.out.csv
    snps = merge_phs_outputs.out.snp
}

workflow prepare_bfiles {
    take:
    snps
    bim_files
    vcf_files

    main:
    if (params.vcf_template) {
        filtered = filter_vcfs(vcf_files, snps)
    } else {
        filtered = filter_bfiles(bim_files)
    }
    merge_list = filtered
        .map { it.find { it.name.endsWith('.bed') }.baseName } // extract basename
        .collectFile(name: 'merge-list.txt', storeDir: params.interim_dir, newLine: true)
    merge_bfiles(merge_list, filtered.collect())

    emit:
    bfiles = merge_bfiles.out.bfiles
}

workflow calculate_phs {
    take:
    bfiles
    model_csv

    main:
    bfiles
    .map { files ->
        def base = files.find { it.name.endsWith('.bed') }.baseName
        tuple(base, files)
    }
    .set { bfiles_with_base }

    plink_score(bfiles_with_base, model_csv)

    emit:
    phs_scores = plink_score.out
}

workflow {
    if (workflow.profile == 'local') {
        //download_models()
    }
    if ( workflow.profile != 'local' ) {
       vcf_files = params.vcf_template
            ? Channel.fromPath(params.vcf_template, checkIfExists:true)
            : Channel.empty()

        bim_files = get_bim_files(vcf_files)

        model = harmonize_models(bim_files)

        bfiles = prepare_bfiles(model.snps, bim_files, vcf_files).bfiles
        
        phs_scores = calculate_phs(bfiles, model.csv)

        phenotypes = preprocess_metadata(file(params.preprocess_script),
  			                file(params.r_config),
    			            file(params.r_functions))
        // run_analysis(phs_scores, phenotypes)
    }
    if (workflow.profile == 'local') {
        // postprocess()
        // knit_manuscript()
    }
} 
