FROM ghcr.io/comorment/r

RUN R -e "install.packages(c('officer', 'kableExtra', 'flextable'), repos='https://cloud.r-project.org')"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
      texlive-xetex \
      texlive-latex-recommended \
      texlive-latex-extra \
      texlive-fonts-recommended \
      texlive-lang-english \
      texlive-plain-generic \
      fonts-dejavu-core \
      fonts-dejavu-extra && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

