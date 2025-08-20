EXTERNAL_DATA = "data/external/"
RAW_DATA = "data/raw/"
INTERIM_DATA = "data/interim/"
PROCESSED_DATA = "data/processed/"
SENSITIVE_DATA = "data/sensitive/"
FIGS = "reports/figures/"
results = "data/results/"

# Constants
icd10pca = "C61"
metastase_loose = c(1,2,4,5,6,7)


# Function definitons
##  Data processing
import_phs_data = function(model_file, inverse = FALSE) {
    readr::read_tsv(model_file, col_names = c("bfile_id", "phs")) %>%
        dplyr::mutate(phs = as.vector(scale(phs, center = TRUE, scale = TRUE))) ->
        phs_data
    if (inverse) {dplyr::mutate(phs_data, phs = phs * -1) -> phs_data}
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

## Performance & Plotting

get_performance_metrics = function(combined_data, bootstrap_iters, model, age = "age", status = "status") {

  response =  paste0("Surv(", age, ", ", status, ")")
  fit = coxph(as.formula(sprintf('%s~%s', response, "phs")), data = combined_data)

  fit_tidy = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  HRSD = list(index = "HRSD",
              value = fit_tidy$estimate[1],
              conf.low = fit_tidy$conf.low[1],
              conf.high = fit_tidy$conf.high[1],
              iters = NULL)

  HR80_20 = get_hr(combined_data,  #default
                  CI = TRUE, 
                  boot = bootstrap_iters)

  HR20_50 = get_hr(combined_data, 
                  lower_interval = c(0.2, 0.80), # 50
                  upper_interval = c(0.0, 0.2),  # 20
                  CI = TRUE, 
                  boot = bootstrap_iters)

  HR80_50 = get_hr(combined_data, 
                  lower_interval = c(0.2, 0.80), # 50
                  upper_interval = c(0.8, 1.0),  # 80
                  CI = TRUE, 
                  boot = bootstrap_iters)

  HR95_50 = get_hr(combined_data, 
                  lower_interval = c(0.2, 0.80), # 50
                  upper_interval = c(0.95, 1.0), # 95
                  CI = TRUE, 
                  boot = bootstrap_iters)
  c_index = get_cindex(combined_data, 
                  CI = TRUE, 
                  boot = bootstrap_iters)
  HR80_20$index = "HR80_20"
  HR20_50$index = "HR20_50"
  HR80_50$index = "HR80_50"
  HR95_50$index = "HR95_50"

  indices = list(HRSD, HR80_20, HR20_50, HR80_50, HR95_50, c_index)

  perform = do.call(rbind.data.frame, lapply(indices, function(x){data.frame(index = x$index, 
                                                                  value = x$value, 
                                                                  conf.low = x$conf.low, 
                                                                  conf.high = x$conf.high)}))
  perform = cbind(model, perform)

  return(perform)
}

plot_km_curves = function(combined_data, age = "age", status = "status") {
  curves = data.frame(curve0_20   = c(0,   0.2),
                      curve20_80  = c(0.2, 0.8),
                      curve80_100 = c(0.8, 1.0))

  label_generator = function(x, y) {
    x = x * 100
    y = y * 100
    out = paste0(x, "-", y, "th")
    return(out)
  }

  km_curves = data.frame()
  for (i in seq_len(length(curves))) {
    curven <- km_curve(data = combined_data,
                       age = age,
                       status = status,
                       interval = curves[[i]],
                       age_range = 40:100,
                       scale = FALSE,
                       inverse = FALSE)
    curven$label = label_generator(curves[1, i], curves[2, i])
    km_curves = rbind(km_curves, curven)
  }


  ggplot(km_curves, aes(x = time,
                        y = estimate,
                        ymin = conf.low,
                        ymax = conf.high,
                        col = label,
                        fill = label)) +
    geom_ribbon(alpha = 0.1,
                color = 0) +
    geom_step() +
    theme_classic() +
    xlim(40, 100) +
    ylim(0, 1) +
    labs(x = "Age", y = "Disease-free Survival") +
    scale_color_brewer(palette = "Set1", direction = -1,
                      name = "PHS Centile") +
    scale_fill_brewer(palette = "Set1", direction = -1,
                      name = "PHS Centile")  +
    theme(aspect.ratio=1)
}

