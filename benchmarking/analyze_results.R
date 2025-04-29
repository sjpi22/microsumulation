###########################  Benchmarking analysis  ################
#
#  Objective: Run benchmarking simulations
########################### <<<<<>>>>> ##############################################

rm(list = ls()) # Clean environment
options(scipen = 999) # View data without scientific notation

#### 1.Libraries and functions  ==================================================

###### 1.1 Load packages
library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
library(grid)
library(patchwork)

###### 1.2 Load functions

# Load functions
distr.sources <- list.files("R", 
                            pattern="*.R$", full.names=TRUE, 
                            ignore.case=TRUE, recursive = TRUE)
sapply(distr.sources, source, .GlobalEnv)


#### 2. General parameters ========================================================

###### 2.1 General parameters
v_labels <- c("cs" = "Cross-sectional", "long" = "Longitudinal", "rcs" = "Repeated cross-sectional")
dodge_width_factor <- 0.7
error_width_factor <- 0.4
plt_size_text <- 18
plt_size_small <- 5/6*plt_size_text
path_benchmark <- "benchmarking/benchmark_output.rds" # Path to save results
path_output <- "benchmarking/result_prevalence.pdf"

#### 3. Pre-processing  ===========================================

# Load results
res <- readRDS(path_benchmark)

# Remove incidence
res$outcomes <- lapply(res$outcomes, function(x) x[method != "incidence"])
res$time <- lapply(res$time, function(x) x[!grepl("incidence", names(x))])

# Calculate maximum time
max_time <- max(unlist(lapply(res$time, function(x) do.call(rbind, x)[, "mean"])))

#### 4. Analysis  ===========================================

dt_stats <- dt_time <- plt_stats <- plt_time <- list()
for (i in 1:length(res$time)) {
  if (i == 1) {
    xlab = "Single interval from age 30 to 80"
  } else {
    xlab = "10-year intervals from age 30 to 80"
  }
  
  # Bar plot of time
  dt_time[[i]] <- do.call(rbind, res$time[[i]])
  
  # Get stats
  dt_stats[[i]] <- res$outcomes[[i]][, `:=` (age_range = age_end - age_start,
                                             age_median = (age_start + age_end)/2,
                                             mean_scaled = (mean - true) / true,
                                             sd_scaled = sd / true)]
  
  # Plot mean value with SD
  dodge_width <- dodge_width_factor*dt_stats[[i]]$age_range[1]
  error_width <- error_width_factor*dt_stats[[i]]$age_range[1]
  plt_stats[[i]] <- ggplot(dt_stats[[i]], 
                           aes(x = age_median, y = mean_scaled, color = method, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(position = position_dodge(width = dodge_width), size = 3) +
    geom_errorbar(aes(ymin = mean_scaled - sd_scaled, ymax = mean_scaled + sd_scaled), 
                  width = error_width, position = position_dodge(width = dodge_width)) +
    scale_x_continuous(breaks = seq(min(dt_stats[[i]]$age_start), max(dt_stats[[i]]$age_end), dt_stats[[i]]$age_range[1])) +
    coord_cartesian(xlim = c(min(dt_stats[[i]]$age_start), max(dt_stats[[i]]$age_end))) +
    scale_color_hue(h = c(180, 300), labels = v_labels, guide = "none") +
    scale_fill_hue(h = c(180, 300), labels = v_labels, guide = "none") +
    labs(x = "Age",
         y = "Percentage difference") +
    theme_bw() + 
    theme(plot.title = element_blank(),
          axis.text.x = element_text(size = plt_size_small),
          axis.text.y = element_text(size = plt_size_small),
          axis.title.y = element_text(size = plt_size_small),
          axis.title = element_text(size = plt_size_small),
          legend.title = element_text(size = plt_size_small),
          legend.text = element_text(size = plt_size_small))
  
  # Bar plot of time
  plt_time[[i]] <- ggplot(data.frame(dt_time[[i]], method = rownames(dt_time[[i]])), 
                          aes(x = method, y = mean, fill = method)) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.5) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_hue(h = c(180, 300), labels = v_labels) +
    coord_cartesian(ylim = c(0, max_time)) +
    labs(fill = "Formulation",
         x = xlab,
         y = "Seconds") +
    theme_bw() + 
    theme(plot.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = plt_size_small),
          axis.title.x = element_text(size = plt_size_small),
          axis.title.y = element_text(size = plt_size_small),
          legend.title = element_text(size = plt_size_small),
          legend.text = element_text(size = plt_size_small))
}


# Plot final patchwork
((wrap_elements(textGrob('Mean difference from true prevalence', gp = gpar(fontsize = plt_size_text))) /
    (plt_stats[[1]] + plt_stats[[2]] + 
       plot_layout(axis_titles = "collect", guides = "collect")) /
    wrap_elements(textGrob('Mean time', gp = gpar(fontsize = plt_size_text)))/
    (plt_time[[1]] + plt_time[[2]] + 
       plot_layout(axis_titles = "collect", guides = "collect"))) +
    plot_layout(axis_titles = "collect", guides = "collect",
                heights = c(0.2, 1, 0.2, 1)) & theme(legend.position = "bottom")) +
  plot_annotation(title = "",
                  theme = theme(plot.title = element_text(size = plt_size_text, face = "bold")))

# Save plot
ggsave(path_output, 
       width = 9, height = 8, dpi = 300, units = "in")
