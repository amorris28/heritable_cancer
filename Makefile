## zip public results for transfer to local
.PHONY: public
public:
	zip -r public.zip data/public/*.tsv data/public/*.rds data/public/*.pdf data/public/*.png

