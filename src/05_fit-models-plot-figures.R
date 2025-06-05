library(tidyverse)
library(stringr)
library(hazrd)
source("src/00_config.R")

#########################################
# First argument: PHS file
# tab delimited file, no header, first column ID, second column PHS score
# Second argument: metadata file
# tab delimited file with column names bfile_id, age, status. Additional
# columns will be ignored.
# Third argument: 
# boolean indicating whether to inverse (x * -1) the PHS scores
# to reverse the direction of effect.

args = commandArgs(trailingOnly=TRUE)

model_file = args[1] # PHS file
# model = gsub("^.*?_", "", basename(model_file))
# model = gsub("\\_.*", "", model) # filename for output
model = gsub("(.+?)(\\_.*)", "\\1", basename(model_file))

metadata_file = args[2] # Phenotype data

inverse = args[3]

tsd_proj = args[4]

output = "reports/"

figs = "reports/figures/"

#########################################
# Import PHS data

phs = read_tsv(model_file, col_names = c("bfile_id", "phs"))

metadata = read_tsv(metadata_file) # Phenotype data

colnames(metadata) = c("bfile_id", "age", "status")

phs$phs = scale(phs$phs, center = TRUE, scale = TRUE)

if (inverse) {phs$phs = phs$phs * -1}

metadata %>%
  left_join(phs, by = "bfile_id") %>%
  filter(!(is.na(phs))) ->
  combined_data

#########################################
# Analyze PHS data

# PHS Histogram
phs_hist = 
	phs_hist(combined_data, normalize = TRUE) +
	labs(title = model) + theme(aspect.ratio=1)

ggsave(paste0(figs, paste(tsd_proj, model, "hist.png", sep = "_")), phs_hist, 
       units = "px", height = 1300, width = 1300)

# Performance metrics


HR80_20 = get_hr(combined_data, CI = TRUE, boot = 300)
ages = c(60, 70, 80)
OR80_20 = lapply(ages, function(x) {get_or(combined_data, or_age = x, CI = TRUE, boot = 300)})
c_index = get_cindex(combined_data, CI = TRUE, boot = 300)

HR80_20$age = NA
c_index$age = NA
indices = list(HR80_20, OR80_20[[1]], OR80_20[[2]], OR80_20[[3]], c_index)

perform = do.call(rbind.data.frame, lapply(indices, function(x){data.frame(index = x$index, 
                                                                 value = x$value, 
                                                                 conf.low = x$conf.low, 
                                                                 conf.high = x$conf.high, 
                                                                 age = x$age)}))
perform = cbind(model, perform)
write_tsv(perform, paste0(output, paste(tsd_proj, model, "performance.tsv", sep = "_")))


# KM Curves
curves = data.frame(curve0_20 = c(0,   0.2),
                    curve20_80 = c(0.2, 0.8),
                    curve80_100 = c(0.8, 1.0))

label_generator = function(x, y) {
  x = x * 100
  y = y * 100
  out = paste0("PHS ", x, "-", y, "th centile")
  return(out)
}

km_curves = data.frame()
for (i in seq_len(length(curves))) {
  curven <- km_curve(data = combined_data,  
                     interval = curves[[i]],
                     age_range = 40:100, 
                     scale = FALSE, 
                     inverse = FALSE)
  curven$label = label_generator(curves[1, i], curves[2, i])
  km_curves = rbind(km_curves, curven)
}


km_curve = ggplot(km_curves, aes(x = time, 
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
  labs(x = "Age", y = "Disease-free Survival", title = model) +
  scale_color_brewer(palette = "Set1", direction = -1,
                     name = "Centile") +
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name = "Centile")  + theme(aspect.ratio=1)

ggsave(paste0(figs, paste(tsd_proj, model, "km_curve.png", sep = "_")), km_curve,
       units = "px", height = 1300, width = 1600)


