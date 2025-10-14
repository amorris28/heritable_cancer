#!/usr/bin/env nextflow
process vcf_to_bim {
    tag "${vcf_file.baseName}"
    storeDir "data/interim/vcf_to_bim"

    input:
    tuple path(vcf_file), val(bim_file)

    output:
    path "${bim_file}.bim"

    script:
    """
    plink --vcf ${vcf_file} --make-just-bim --out ${bim_file}
    """
}

process harmonize_phs {
    tag "${bim_file.baseName}"
    container params.python3_image

    input:
    path script_file
    path phs_file
    path bim_file
    val build

    output:
    path "phs601_${params.cluster_name}_${bim_file.baseName}.csv", emit: harmonized_phs

    script:
    """
    python3 ${script_file} \
        --phs_file ${phs_file} \
        --bim_file ${bim_file} \
        --build ${build} \
        --out phs601_${params.cluster_name}_${bim_file.baseName}.csv
    """
}

process merge_phs_outputs {
    tag "merge_$params.cluster_name"
    publishDir "data/interim", mode: 'copy'

    input:
    path csv_files

    output:
    path "phs601_${params.cluster_name}.csv", emit: csv
    path "phs601_${params.cluster_name}.snps", emit: snp

    script:
    """
    awk 'FNR==1 && NR!=1{next} {print}' ${csv_files.join(' ')} > phs601_${params.cluster_name}.csv
    rm -f ${csv_files.join(' ')}
    cat phs601_${params.cluster_name}.csv | awk '{print \$6}' | tail -n +1 > phs601_${params.cluster_name}.snps
    """
}

process extract_snps_p212 {
    storeDir "data/interim/extract_snps"

    input:
    path bfiles
    path snp_file
    path psam_file

    output:
    path "${input_file.getBaseName(2)}.phs601.{bed,bim,fam}", emit: bfiles
    
    script:
    """
    plink2 \
        --bfile ${input_file} \
        --psam ${psam_file} \
        --extract ${snp_file} \
        --make-bed \
        --out ${input_file.getBaseName(2)}.phs601
    """
}

process extract_snps_p2521 {
    storeDir "data/interim/extract_snps"

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
        --out ${vcf_file.getBaseName(2)}.phs601
    """
}
 process make_merge_list {
    input:
    val bfiles_basenames
     
    output:
    path "merge-list.txt"
    
    script:
    """
    for dir in ${bfiles_basenames}; do
        echo \${dir} >> merge-list.txt
    done
    """
}
process merge_chr {
    container params.gwas_image

    input:
    path merge_list
    path bfiles
     
    output:
    tuple val("${params.cluster_name}_phs601"), path("${params.cluster_name}_phs601.{bed,bim,fam}"), emit: bfiles

    script:
    """
    plink --merge-list ${merge_list} --make-bed --out ${params.cluster_name}_phs601
    """
}
process plink_score {
    publishDir "data/processed", mode: 'copy' 
    container params.gwas_image

    input:
    tuple val(base), path(bfiles)
    path phs_csv
     
    output:
    path "${params.cluster_name}_phs601_phs.profile"
    
    script:
    """
    plink \
        --bfile ${base} \
        --score ${phs_csv} 6 7 5 header sum \
        --out ${params.cluster_name}_phs601_phs
    """
}
process preprocess_metadata {
    publishDir "data/processed", mode: 'copy'
    container params.r_image

    input:
    path preprocess_script
    path 'src/config.R'
    path 'src/functions.R'

    output:
    path "data/processed/${params.cluster_name}_phenotypes.tsv"

    script:
    """
    mkdir -p data/processed
    Rscript ${preprocess_script}
    """
}

workflow calculate_phs {
    if ( params.cluster_name == 'p2521' ) {
        vcf_files = Channel.fromPath(params.vcf_template)
        bim_info = vcf_files.map { vcf_file ->
            def bim_file = vcf_file.name.replaceFirst(/\.vcf(\.gz)?$/, '')  // handles .vcf or .vcf.gz
            tuple(vcf_file, bim_file)
        }
        bim_files = vcf_to_bim(bim_info)
    } else {
        bim_files = Channel.fromPath(params.bim_template)
    }

    harmonized_chr = harmonize_phs(file(params.harmonize_script), 
                                   file(params.phs_file), 
                                   bim_files, 
                                   params.build)
    harmonized_chr
        .collect()
        .set { all_csvs }

    merge_phs_outputs(all_csvs)

    if ( params.cluster_name == 'p212' ) {
        bim_ch = Channel.fromPath(params.bim_template)

        bim_ch
            .map { bimPath ->
                def base = bimPath.baseName
                def bim  = file("${bimPath.parent}/${bimPath.baseName}.bim")
                def bed = file("${bimPath.parent}/${bimPath.baseName}.bed")
                def fam = file("${bimPath.parent}/${bimPath.baseName}.fam")
                tuple(base, bim, bed, fam)
          }
           .set { bfile_triples }
          .view()
        bfiles_ch = extract_snps_p212(bfile_triplets, merge_phs_outputs.out.snp, psam_file)
    }
    else if ( params.cluster_name == 'p2521' ) {
        vcf_files = Channel.fromPath(params.vcf_template)
        // vcf_files = vcf_files_ch.map { it.name.replaceFirst(/\.vcf(\.gz)?$/, '') }
        bfiles_ch = extract_snps_p2521(vcf_files, merge_phs_outputs.out.snp)
    }
    // Flatten the triplets to get the .bed file for each chromosome
    
    merge_list = bfiles_ch
        .map { it.find { it.name.endsWith('.bed') }.baseName } // extract basename
        .collectFile(name: 'merge-list.txt', storeDir: params.interim_dir, newLine: true)
    bfiles_ch.collect().set { bfiles }
    merge_chr(merge_list, bfiles)
    plink_score(merge_chr.out.bfiles, merge_phs_outputs.out.csv)
}

workflow {
    if ( params.cluster_name != 'local' ) {
        calculate_phs()
        preprocess_metadata(file(params.preprocess_script),
  			                file(params.r_config),
    			            file(params.r_functions))
        // run_analysis()
    }
    else if ( params.cluster_name == 'local' ) {
        // postprocess()
        // knit_manuscript()
    }
    else {
        error "Invalid environment: ${params.cluster_name}"
    }
} 
