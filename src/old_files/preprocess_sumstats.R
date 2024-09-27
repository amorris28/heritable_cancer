readr::read_tsv('raw_data/sumstats/PGS000067_hmPOS_GRCh37.txt.gz', comment = "#")
readr::write_tsv('derived_data/sumstats/PRACTICAL_PCa_Seibert2018.txt.gz')