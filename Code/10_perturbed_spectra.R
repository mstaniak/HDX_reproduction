# Libraries ----
library(igraph)
library(data.table)
library(ggplot2)
library(parallel)
# Functions ---- 
get_comparison = function(parameters, ps_m, observed_data, und, num_pars = NULL) {
  if (is.null(num_pars)) {
    num_pars = length(parameters)
  }
  comp = rbind(IsoHDX3:::getExpectedSpectra(parameters, ps_m, num_pars, observed_data, und)[, .(Peptide, Time, IntDiff, Intensity = ExpectedPeak, Type = "Fitted")],
               observed_data[, .(Peptide, Time, IntDiff, Intensity, Type = "observed")])
  comp
}
plot_comp = function(comp) {
  ggplot(comp, aes(x = reorder(as.character(IntDiff), IntDiff), y = Intensity, fill = Type, color = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(~Time) +
    theme_bw() +
    theme(legend.position = "bottom")
}


get_probs = function(pars) {
  estim_pars_by_seg = IsoHDX3:::getSegmentParametersFromBetas(pars, unique(cl_data[, .(Segment, MaxUptake)])$MaxUptake + 1)
  estim_seg_probs_by_time = IsoHDX3:::getSegmentProbabilitiesFromParams(estim_pars_by_seg, times)
  
  estim_probs_by_time_dt = rbindlist(lapply(seq_along(times), function(ith_time) {
    rbindlist(lapply(seq_along(colnames(ps_m)[-ncol(ps_m)]), function(ith_seg) {
      n_ex = seg_dt[SegmentSigned == colnames(ps_m)[ith_seg], SegmentLength + 1]
      list(Time = rep(times[ith_time], n_ex),
           Segment = rep(colnames(ps_m)[-ncol(ps_m)][ith_seg], n_ex),
           NumExchanged = seq_len(n_ex) - 1,
           Probability = estim_seg_probs_by_time[[ith_time]][[ith_seg]])
    }))
  }))
  estim_probs_by_time_dt
}

# Simulated data 
dataset2_expected_spectra = readRDS("./Data/Simulated/DataSet2/expected_spectra.RDS")
dataset2_expected_spectra_perturbed = readRDS("./Data/Simulated/DataSet2/perturbed_expected_spectra.RDS")
ps_m = readRDS("./Data/Simulated/DataSet2/ps_m.RDS")
times = readRDS("./Data/Simulated/DataSet2/times.RDS")
cl_data = readRDS("./Data/Simulated/DataSet2/cl_data.RDS")

# Simulation results

# Comparison of clean and perturbed spectra

dataset2_perturbed_results = readRDS("./RawResults/dataset2_perturbed_results.RDS")

ggplot(rbind(cbind(dataset2_expected_spectra, Spectra = "clean"),
             cbind(dataset2_expected_spectra_perturbed, Spectra = "perturbed")), 
       aes(x = reorder(as.character(IntDiff), IntDiff), y = Intensity, fill = Spectra)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Peptide ~ Time) +
  theme_bw() +
  scale_fill_manual(name = "spectra type", values = viridis::viridis(3)[-(3)]) +
  scale_x_discrete(breaks = as.character(seq(0, 13, by = 4))) +
  xlab("mass shift") +
  ylab("intensity") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Figures/hdx_noisy_spectra_comp.png", device = "png", 
       scale = 1, width = 15, height = 15, units = "in", dpi = 300)


probs_estim_list = unlist(dataset2_perturbed_results, F, F)
probs_estim_list = rbindlist(lapply(probs_estim_list, function(x) {
  if (!is.null(x$fitted[[1]]) & !is.null(x$fitted[[2]])) {
    if (max(abs(x$fitted[[1]]$fvec)) < 1e-6) {
      probs_clean = get_probs(x$fitted[[1]]$x)
      probs_clean[, Spectra := "clean"]
    } else {
      probs_clean = NULL
    }
    if (max(abs(x$fitted[[2]]$fvec)) < 1e-6) {
      probs_noisy = get_probs(x$fitted[[2]]$x)
      probs_noisy[, Spectra := "noisy"]
    } else {
      probs_noisy = NULL
    }
    probs_comp = rbind(probs_clean, probs_noisy)
    if (nrow(probs_comp) > 0) {
      probs_comp[, error_sd := x$error_sd]
      probs_comp[, iter := x$iter]
      probs_comp
    } else {
      NULL
    }
  } else {
    NULL
  }
}))


ggplot(probs_estim_list[error_sd == 0.01],
       aes(x = reorder(as.character(Time), Time), y = Probability, color = Spectra)) +
  geom_boxplot(size = 1.1) +
  geom_point(aes(x = reorder(as.character(Time), Time), y = Probability),
             data = probs_by_time_dt, color = "red", size = 1.2) +
  scale_color_manual(name = "spectra type", values = viridis::viridis(3)[-1]) +
  xlab("time") +
  ylab("exchange probability") +
  facet_grid(Segment ~ NumExchanged) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Figures/hdx_30par_noisy_001.png", device = "png", scale = 1, width = 10, height = 10, units = "in", dpi = 300)

ggplot(probs_estim_list[error_sd == 0.01][Segment == unique(Segment)[4]],
       aes(x = reorder(as.character(Time), Time), y = Probability, color = Spectra)) +
  geom_boxplot(size = 1.1) +
  geom_point(aes(x = reorder(as.character(Time), Time), y = Probability),
             data = probs_by_time_dt[Segment == unique(Segment)[4]], color = "red", size = 1.2) +
  scale_color_manual(name = "spectra type", values = viridis::viridis(3)[-1]) +
  xlab("time") +
  ylab("exchange probability") +
  facet_grid(Segment ~ NumExchanged) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
ggsave("Figures/hdx_30par_noisy_001_single_seg.png", device = "png", scale = 1, width = 10, height = 10, units = "in", dpi = 300)

conv_estim_list = rbindlist(lapply(unlist(dataset2_perturbed_results, F, F), function(x) {
  if (!is.null(x$fitted[[1]]) & !is.null(x$fitted[[2]])) {
    if (max(abs(x$fitted[[1]]$fvec)) < 1e-6) {
      conv_info_clean = data.table(NumIter = x$fitted[[1]]$iter)
      conv_info_clean[, Spectra := "clean"]
    } else {
      conv_info_clean = NULL
    }
    if (max(abs(x$fitted[[2]]$fvec)) < 1e-6) {
      conv_info_noisy = data.table(NumIter = x$fitted[[2]]$iter)
      conv_info_noisy[, Spectra := "noisy"]
    } else {
      conv_info_noisy = NULL
    }
    probs_comp = rbind(conv_info_clean, conv_info_noisy)
    if (!is.null(probs_comp) > 0) {
      if (nrow(probs_comp) > 0) {
        probs_comp[, error_sd := x$error_sd]
        probs_comp[, iter := x$iter]
        probs_comp
      } else {
        NULL
      }
    } else {
      NULL
    }
  } else {
    NULL
  }
}))

conv_estim_list[, .(convd = .N,
                    mean_iter = mean(NumIter)),
                by = c("Spectra", "error_sd")]

