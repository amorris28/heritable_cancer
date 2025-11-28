FROM ghcr.io/comorment/r

RUN R -e "install.packages(c('kableExtra'), repos='https://cloud.r-project.org')"
