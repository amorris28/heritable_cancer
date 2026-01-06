library(boot)
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

## Performance and Plotting

performance_metrics <- function(data, model, age = "age", status = "status") {
  df <- data.frame(phs = data$phs, 
                   age = data[[age]], 
                   status = data[[status]])
  print(head(df))
  # HR80_50
  boot_hr8050 <- function(data, indices) {
    boot_sample <- data[indices, ]
    calc_hr(boot_sample, upper_interval = c(0.80, 1.0), lower_interval = c(0.2, 0.80))
  }
  boot_results <- boot(df, boot_hr8050, R = 999)
  ci_results <- boot.ci(boot_results, type = "perc")
  HR80_50 <- list(
    index = "HR80_50",
    value = boot_results$t0,
    conf.low = ci_results$percent[4],
    conf.high = ci_results$percent[5]
  )

  # HR80_20
  boot_hr8020 <- function(data, indices) {
    boot_sample <- data[indices, ]
    calc_hr(boot_sample, upper_interval = c(0.80, 1.0), lower_interval = c(0.0, 0.20))
  }
  boot_results <- boot(df, boot_hr8020, R = 999)
  ci_results <- boot.ci(boot_results, type = "perc")
  HR80_20 <- list(
    index = "HR80_20",
    value = boot_results$t0,
    conf.low = ci_results$percent[4],
    conf.high = ci_results$percent[5]
  )
  # C_Index
  boot_cindex <- function(data, indices) {
    boot_sample <- data[indices, ]
    calc_cindex(boot_sample)
  }
  boot_results <- boot(df, boot_cindex, R = 999)
  ci_results <- boot.ci(boot_results, type = "perc")
  C_Index <- list(
    index = "C_Index",
    value = boot_results$t0,
    conf.low = ci_results$percent[4],
    conf.high = ci_results$percent[5]
  )
  
  indices <- list(HR80_20, HR80_50, C_Index)

  perform <- do.call(rbind.data.frame, lapply(indices, function(x) {
    data.frame(
      index = x$index,
      value = x$value,
      conf.low = x$conf.low,
      conf.high = x$conf.high
    )
  }))
  perform <- cbind(model, perform)
}

plot_km_curves <-
  function(
    combined_data, age = "age", status = "status", style = "surv"
  ) {
    # Style must be one of "surv" or "risk"
    curves <- data.frame(
      curve0_20 = c(0, 0.2),
      curve20_80 = c(0.2, 0.8),
      curve80_100 = c(0.8, 1.0)
    )

    label_generator <- function(x, y) {
      x <- x * 100
      y <- y * 100
      paste0(x, "-", y, "th")
    }

    km_curves <- data.frame()
    for (i in seq_len(length(curves))) {
      curven <- km_curve(
        data = combined_data,
        age = age,
        status = status,
        interval = curves[[i]],
        age_range = 40:100,
        scale = FALSE,
        inverse = FALSE
      )
      curven$label <- label_generator(curves[1, i], curves[2, i])
      km_curves <- rbind(km_curves, curven)
    }
    # cutoff_data <- km_curves %>%
    #   group_by(label) %>%
    #   filter(estimate == min(estimate)) %>%
    #   summarise(cutoff_time = min(time)) %>%
    #   right_join(km_curves, by = "label") %>%
    #   filter(time <= cutoff_time)
    yaxis <- "Disease-free Survival"
    # ylims = c(0, 1)
    if (style == "risk") {
      km_curves %>%
        mutate(
          estimate = (1 - estimate) * 100,
          conf.high = (1 - conf.high) * 100,
          conf.low = (1 - conf.low) * 100
        ) ->
      km_curves
      yaxis <- "Cumulative incidence, %\n(Absolute risk)"
    }


    ggplot(km_curves, aes(
      x = time,
      y = estimate,
      ymin = conf.low,
      ymax = conf.high,
      col = label,
      fill = label
    )) +
      geom_ribbon(
        alpha = 0.1,
        color = 0
      ) +
      geom_step() +
      theme_classic(base_size = 12) +
      # lims(x = c(40, 100), y = ylims) +
      lims(x = c(40, 100)) +
      labs(x = "Age", y = yaxis) +
      scale_color_brewer(
        palette = "Set1", direction = -1,
        name = "PHS Centile"
      ) +
      scale_fill_brewer(
        palette = "Set1", direction = -1,
        name = "PHS Centile"
      ) +
      theme(aspect.ratio = 1)
  }
