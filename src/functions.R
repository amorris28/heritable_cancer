# Function definitons
##  Data processing
import_phs_data <- function(model_file, inverse = FALSE) {
  readr::read_table(model_file) %>%
    dplyr::select(bfile_id = IID, phs = SCORESUM) %>%
    dplyr::mutate(phs = as.vector(scale(phs, center = TRUE, scale = TRUE))) ->
  phs_data
  if (inverse) {
    dplyr::mutate(phs_data, phs = phs * -1) -> phs_data
  }
  return(phs_data)
}

process_combined_data <- function(metadata, phs_data) {
  combined_data <- metadata %>%
    left_join(phs_data, by = "bfile_id") %>%
    drop_na(phs, age, status)

  return(combined_data)
}

##  Formatting
pvalr <- function(pvals, sig.limit = .001, digits = 3, html = FALSE) {
  roundr <- function(x, digits = 1) {
    res <- sprintf(paste0("%.", digits, "f"), x)
    zzz <- paste0("0.", paste(rep("0", digits), collapse = ""))
    res[res == paste0("-", zzz)] <- zzz
    res
  }

  sapply(pvals, function(x, sig.limit) {
    if (x < sig.limit) {
      if (html) {
        return(sprintf("&lt; %s", format(sig.limit)))
      } else {
        return(sprintf("< %s", format(sig.limit)))
      }
    }
    if (x > .1) {
      return(roundr(x, digits = 2))
    } else {
      return(roundr(x, digits = digits))
    }
  }, sig.limit = sig.limit)
}


