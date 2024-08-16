
sumstats = readr::read_tsv('PRACTICAL_PCa_Seibert2018.GRCh37.gz')

antegenes = genio::read_plink('out/antegenes_imputed_allchr-temporary_cs_PCa')$X

sumstats[sumstats$RSID == "rs201305416", "RSID"] = "rs10905266"

readr::write_tsv(sumstats, 'PRACTICAL_PCa_Seibert2018b.GRCh37.gz')