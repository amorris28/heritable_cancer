EXTERNAL_DATA = "data/external/"
RAW_DATA = "data/raw/"
INTERIM_DATA = "data/interim/"
PROCESSED_DATA = "data/processed/"
SENSITIVE_DATA = "data/sensitive/"
FIGS = "reports/figures/"
results = "data/results/"


pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {

  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0('%.', digits, 'f'), x)
    zzz <- paste0('0.', paste(rep('0', digits), collapse = ''))
    res[res == paste0('-', zzz)] <- zzz
    res
  }

  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit)
      if (html)
        return(sprintf('&lt; %s', format(sig.limit))) else
          return(sprintf('< %s', format(sig.limit)))
    if (x > .1)
      return(roundr(x, digits = 2)) else
        return(roundr(x, digits = digits))
  }, sig.limit = sig.limit)
}

# Constants
icd10pca = "C61"
metastase_loose = c(1,2,4,5,6,7)

