#!/usr/bin/env nextflow
process harmonize_phs {
    tag "$cluster_name"
    // debug true

    input:
    tuple val(cluster_name), val(out_file), path(bim_file)
    path script_file
    path phs_file
    val build

    output:
    path out_file

    script:
    """
    python3 ${script_file} \
        --phs_file ${phs_file} \
        --bim_file ${bim_file} \
        --build ${build} \
        --out ${out_file}
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


workflow {
    output_file = "phs601_${params.cluster_name}.csv"
    script_file = file('src/01_preprocess/PHS_PC.py')

    // Expand the bim template pattern into actual files
    def bim_files = Channel.fromPath(params.bim_template)

    // Create tuples: (cluster_name, output_file, bim_file)
    cluster_info_ch = bim_files.map { bim_file ->
        // extract chromosome if present
        def chr_id = bim_file.baseName.replaceAll(/[^0-9]+/, '')
        def suffix = chr_id ? "_chr${chr_id}" : ""
        def out_file = "phs601_${params.cluster_name}${suffix}.csv"
        tuple(params.cluster_name, out_file, bim_file)
    }
    harmonized_ch = harmonize_phs(cluster_info_ch, script_file, file(params.phs_file), params.build)

    // Wait until all harmonized outputs are done, then merge
    harmonized_ch
        .collect()  // gather all files into a list
        .set { all_csvs }

    merge_phs_outputs(params.cluster_name, all_csvs)
}

