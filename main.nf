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
    path "phs601_${params.cluster}_${bim_file.baseName}.csv", emit: harmonized_phs

    script:
    """
    python3 ${script_file} \
        --phs_file ${phs_file} \
        --bim_file ${bim_file} \
        --build ${build} \
        --out phs601_${params.cluster}_${bim_file.baseName}.csv
    """
}

process merge_phs_outputs {
    tag "merge_$params.cluster"
    publishDir params.processed_dir, mode: 'copy'

    input:
    path csv_files

    output:
    path "phs601_${params.cluster}.csv", emit: csv
    path "phs601_${params.cluster}.snps", emit: snp

    script:
    """
    awk 'FNR==1 && NR!=1{next} {print}' ${csv_files.join(' ')} > phs601_${params.cluster}.csv
    cat phs601_${params.cluster}.csv | awk '{print \$6}' | tail -n +1 > phs601_${params.cluster}.snps
    """
}

process filter_bfiles {
    storeDir "${params.interim_dir}/filter_bfiles"

    input:
    tuple val(base), path(bfiles)
    path snp_file

    output:
    path "${base}.phs601.{bed,bim,fam}", emit: bfiles
    
    script:
    """
    plink2 \
        --bfile ${base} \
        ${params.psam_flag} \
        --extract ${snp_file} \
        --make-bed \
        --out ${base}.phs601
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
    path "${params.cluster}_phs601.{bed,bim,fam}", emit: bfiles

    script:
    """
    plink \
        --merge-list ${merge_list} \
        --make-bed \
        --out ${params.cluster}_phs601
    """
}
process plink_score {
    publishDir params.processed_dir, mode: 'copy' 
    container params.gwas_image

    input:
    tuple val(base), path(bfiles)
    path phs_csv
     
    output:
    path "${params.cluster}_phs601.profile"
    
    script:
    """
    plink \
        --bfile ${base} \
        --score ${phs_csv} 6 7 5 header sum \
        --out ${params.cluster}_phs601
    """
}
process preprocess_metadata {
    tag "${params.cluster}"

    publishDir params.processed_dir, 
               mode: 'copy', 
               pattern: "phenotypes.tsv"
    publishDir params.public_dir, 
               mode: 'copy', 
               pattern: "all_variants.tsv"
    publishDir params.public_dir,
               mode: 'copy',
               pattern: "${params.cluster}_age_summary.tsv"

    input:
    path preprocess_script
    path model_file
    path 'src/config.R'
    path 'src/functions.R'
    path vcf

    output:
    path "phenotypes.tsv", emit: phenotypes
    path "all_variants.tsv", optional: true, emit: all_variants
    path "${params.cluster}_age_summary.tsv", optional: true, emit: age_summary

    script:
    """
    mkdir -p data/processed
    Rscript ${preprocess_script} \
        ${model_file} \
        FALSE
    mv age_summary.tsv ${params.cluster}_age_summary.tsv
    """
}
process fit_models {
    publishDir params.public_dir, mode: 'copy'

    input:
    path metadata_file
    path 'src/config.R'
    path 'src/functions.R'
    path analysis_script

    output:
    path "${params.cluster}_phs_zph_tidy.tsv"
    path "${params.cluster}_any_pca_performance.tsv"
    path "${params.cluster}_metastatic_performance.tsv"
    path "${params.cluster}_sample_size.tsv"
    path "${params.cluster}_coxph.tsv"

    script:
    """
    Rscript ${analysis_script} \
        ${metadata_file} \
        TRUE \
        ${params.cluster} \
        ${params.bootstrap_iters}
    """
}
process plot_figures {
    publishDir params.public_dir, mode: 'copy'

    input:
    path metadata_file
    path 'src/config.R'
    path 'src/functions.R'
    path analysis_script

    output:
    path "${params.cluster}_any_pca_abs_risk.rds"
    path "${params.cluster}_any_pca_risk_table.rds"
    path "${params.cluster}_metastatic_abs_risk.rds"
    path "${params.cluster}_metastatic_risk_table.rds"
    path "${params.cluster}_zph_tt.rds"
    path "${params.cluster}_zph_tt_aggr.rds"

    script:
    """
    Rscript ${analysis_script} \
        ${metadata_file} \
        ${params.cluster}
    """
}


process run_analysis_incidence {
    publishDir params.public_dir, mode: 'copy'

    input:
    path metadata_file
    path 'src/config.R'
    path 'src/functions.R'
    path incidence_script
    path 'data/external/ASIRs_per100k_Huynh-Le_2020.xlsx'

    output:
    path "${params.cluster}_cumulative_comparison.rds"
    path "${params.cluster}_stratified_incidence.rds"
    path "${params.cluster}_cumulative_comparison_metastatic.rds"
    path "${params.cluster}_stratified_incidence_metastatic.rds"

    script:
    """
    Rscript ${incidence_script} \
        ${metadata_file} \
        ${params.cluster}
    """
}
process run_analysis_wgs {
    publishDir params.public_dir, mode: 'copy'

    input:
    path wgs_script
    path metadata_file
    path 'src/config.R'
    path 'src/functions.R'

    output:
    path "${params.cluster}_wgs_summary.tsv"
    path "${params.cluster}_fisher.tsv"
    path "${params.cluster}_wgs_phs_performance_95.tsv"
    path "${params.cluster}_binomial_phs_hoxb13.tsv"
    path "${params.cluster}_wgs.tsv"
    path "${params.cluster}_wgs_cindex.tsv"
    path "${params.cluster}_wgs_phs_abs_risk.rds"
    path "${params.cluster}_wgs_phs_risk_table.rds"
    path "${params.cluster}_wgs_hoxb13_abs_risk.rds"
    path "${params.cluster}_wgs_hoxb13_risk_table.rds"

    script:
    """
    Rscript ${wgs_script} \
        ${metadata_file} \
        TRUE \
        ${params.cluster} \
        ${params.bootstrap_iters}
    """
}
process postprocess_figures {
    publishDir params.manuscript_dir, mode: 'copy'
    container params.r_image

    input:
    path make_figures_script
    path 'data/public'
    path 'src/config.R'
    path 'src/functions.R'

    output:
    path "figures_and_tables"
//    path "km_curve.{png,pdf}"
//    path "metastatic_km_curve.{png,pdf}"
//    path "stratified_incidence_agespecific.{png,pdf}"
//    path "cumulative_comparison.{png,pdf}"
//    path "p212_phs601_wgs_km_curve.{png,pdf}"
//    path "table_1.tsv"

    script:
    """
    mkdir -p figures_and_tables
    Rscript ${make_figures_script}
    """
}
process knit_manuscript {
    publishDir params.manuscript_dir, mode: 'copy'
    container params.r_image

    input:
    path manuscript_file
    path 'src/config.R'
    path 'src/functions.R'
    path 'reference.docx'
    path 'references.bib'
    path csl
    path '_extensions'
    path 'data/public'
    path 'manuscript/figures_and_tables'
    path 'header.tex'

    output:
    path "manuscript.pdf"
    path "manuscript.docx"

    script:
    """
    quarto render ${manuscript_file} --to docx --metadata csl=${csl}
    quarto render ${manuscript_file} --to pdf --metadata csl=${csl}
    """
}
process get_hoxb13 {
    publishDir params.processed_dir, mode: 'copy'
    container params.gwas_image

    input:
    tuple path(bed), path(bim), path(fam)
    // tuple path(bim), path(bed), path(fam)

    output:
    path "tmp_rs138213197.raw"

    script:
    def prefix = bim.baseName
    """
    plink --bfile ${prefix} \
      --snp rs138213197 \
      --recode A \
      --out tmp_rs138213197
    """
    // echo rs138213197 >> rs138213197.txt
    // echo "rs138213197 C" > ref_allele.txt
    // plink2 \
    //     --bfile ${prefix} \
    //     ${params.psam_flag} \
    //     --extract rs138213197.txt \
    //     --ref-allele ref_allele.txt \
    //     --out rs138213197 \
    //     --recode vcf

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
        bfile_tuples = bim_files
            .map { bim ->
                def base = bim.baseName
                def dir = bim.parent
                def files = [
                    file("${dir}/${base}.bed"),
                    file("${dir}/${base}.bim"),
                    file("${dir}/${base}.fam")
                ]
                tuple(base, files)
            }
        filtered = filter_bfiles(bfile_tuples, snps)
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
    if ( params.cluster != 'local' ) {
       vcf_files = params.vcf_template
            ? Channel.fromPath(params.vcf_template, checkIfExists:true)
            : Channel.empty()

        bim_files = get_bim_files(vcf_files)

        model = harmonize_models(bim_files)

        bfiles = prepare_bfiles(model.snps, bim_files, vcf_files).bfiles
        
        phs_scores = calculate_phs(bfiles, model.csv)

        if ( params.cluster == 'p212' ) {
            chr17_bim = Channel.fromPath(params.bim_template)
                .filter { it.name.contains('chr17.') }
                .first()  // Ensures you get only one file
            chr17_plink_files = chr17_bim.map { bim ->
                def prefix = bim.toString().replaceAll(/\.bim$/, '')
                tuple(
                    file("${prefix}.bim"),
                    file("${prefix}.bed"),
                    file("${prefix}.fam")
                )
            }
            // hoxb13 = get_hoxb13(chr17_plink_files)
            hoxb13 = get_hoxb13(bfiles)
        } else {
            hoxb13 = Channel.of(file("$projectDir/assets/NO_FILE"))
        }
        preprocess = preprocess_metadata(file(params.preprocess_script),
                            phs_scores,
  			                file(params.r_config),
    			            file(params.r_functions),
                            hoxb13)
        fit_models(preprocess.phenotypes,
  		         file(params.r_config),
    		     file(params.r_functions),
                     file(params.fit_models_script))
        plot_figures(preprocess.phenotypes,
  			         file(params.r_config),
    			     file(params.r_functions),
                     file(params.plot_figures_script))
        run_analysis_incidence(preprocess.phenotypes,
  			                   file(params.r_config),
    			               file(params.r_functions),
                               file(params.incidence_script),
                               file(params.asir_file))
        if ( params.cluster == 'p212' ) {
            run_analysis_wgs(file(params.wgs_script),
                             preprocess.phenotypes,
  	    		             file(params.r_config),
        			         file(params.r_functions))
        }
    }
    if (params.cluster == 'local') {
        postprocess_figures(file(params.figures_script), 
                            file(params.public_dir),
  			                file(params.r_config),
    			            file(params.r_functions))

        knit_manuscript(file(params.manuscript_file),
  			            file(params.r_config),
    			        file(params.r_functions),
                        file(params.reference_docx),
                        file(params.bibtex_file),
                        file(params.csl_file),
                        file(params.extensions_dir),
                        file(params.public_dir),
                        postprocess_figures.out,
                        file(params.header_tex))
    }
} 
