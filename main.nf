#!/usr/bin/env nextflow
process vcf_to_bim {
    input:
    path vcf_file

    output:
    path bim_file

    script:
    """
    plink --vcf ${vcf_file} --make-just-bim --out ${bim_file}
    """
}

process clean_phs {
    tag "$cluster_name"
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
    tag "$cluster_name"
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

workflow preprocess_vcf {
    vcf_files = Channel.fromPath(params.vcf_file)
    vcf_to_bim(vcf_files)
}

workflow harmonize_phs {
    cluster_name = params.cluster_name
    script_file = file('src/PHS_PC.py')
    phs_file = file(params.phs_file)
    bim_files = Channel.fromPath(params.bim_template)
    build = params.build

    harmonized_chr = clean_phs(cluster_name, script_file, phs_file, bim_files, build)

    harmonized_chr
        .collect()
        .set { all_csvs }

    merge_phs_outputs(cluster_name, all_csvs)
}

workflow {
    if ( params.cluster_name == 'p2521' ) {
        preprocess_vcf()
    }
    harmonize_phs()
}
