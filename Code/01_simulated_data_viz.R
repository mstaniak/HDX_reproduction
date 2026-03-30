# Libraries -----
library(data.table)
library(ggplot2)

# Functions ----

plot_true_probs_in_time = function(probs_by_time_dt) {
  ggplot(probs_by_time_dt, aes(x = Time, y = Probability,
                               group = reorder(as.character(NumExchanged), NumExchanged),
                               color = reorder(as.character(NumExchanged), NumExchanged))) +
    geom_line(size = 1.2) +
    geom_point(size = 1.4) +
    scale_color_discrete(name = "no. exchanged") +
    scale_x_continuous(breaks = c(0, 5, 10)) +
    xlab("time") +
    ylab("probability") +
    facet_grid(~Segment) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          strip.text = element_text(size = 10)) 
}

plot_expected_spectra = function(expected_spectra) {
  ggplot(expected_spectra, aes(x = IntDiff, ymin = 0, ymax = Intensity)) +
    geom_linerange(size = 1.2) +
    facet_grid(Peptide ~ Time) +
    ylab("intensity") +
    xlab("mass shift") +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12),
          strip.text = element_text(size = 9)) 
}


# Simulated data ----
dataset1_probs_by_time = readRDS("./Data/Simulated/DataSet1/probs_by_time_dt.RDS")
dataset2_probs_by_time = readRDS("./Data/Simulated/DataSet2/probs_by_time_dt.RDS")
dataset3_probs_by_time = readRDS("./Data/Simulated/DataSet3/probs_by_time_dt.RDS")

dataset1_expected_spectra = readRDS("./Data/Simulated/DataSet1/expected_spectra.RDS")
dataset2_expected_spectra = readRDS("./Data/Simulated/DataSet2/expected_spectra.RDS")
dataset3_expected_spectra = readRDS("./Data/Simulated/DataSet3/expected_spectra.RDS")

plot_true_probs_in_time(dataset1_probs_by_time)
ggsave("./Figures/hdx_true_probs_smallcl.png", device = "png", width = 10, height = 10, units = "in", scale = 1, dpi = 300)
plot_true_probs_in_time(dataset2_probs_by_time)
ggsave("./Figures/hdx_true_probs_medcl.png", device = "png", width = 10, height = 10, units = "in", scale = 1, dpi = 300)
plot_true_probs_in_time(dataset3_probs_by_time)
ggsave("./Figures/hdx_true_probs_shortsegm.png", device = "png", width = 10, height = 10, units = "in", scale = 1, dpi = 300)

plot_expected_spectra(dataset1_expected_spectra)
ggsave("./Figures/simulated_spectra_smallcl.png", device = "png", width = 10, height = 10, units = "in", scale = 1, dpi = 300)
plot_expected_spectra(dataset2_expected_spectra)
ggsave("./Figures/simulated_spectra_medcl.png", device = "png", width = 10, height = 10, units = "in", scale = 1, dpi = 300)
plot_expected_spectra(dataset3_expected_spectra)
ggsave("./Figures/simulated_spectra_shortsegm.png", device = "png", width = 10, height = 10, units = "in", scale = 1, dpi = 300)
