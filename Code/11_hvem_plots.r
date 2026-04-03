# Libraries ---- 
library(ggplot2)
library(data.table)

# Functions ----
get_comparison = function(parameters, ps_m, observed_data, und, num_pars = NULL) {
  if (is.null(num_pars)) {
    num_pars = length(parameters)
  }
  comp = rbind(IsoHDX3:::getExpectedSpectra(parameters, ps_m, num_pars, observed_data, und)[, .(Peptide, Charge, Rep, Time, IntDiff, Intensity = ExpectedPeak, Type = "Fitted")],
               observed_data[, .(Peptide, Charge, Rep, Time, IntDiff, Intensity, Type = "observed")])
  comp
}

plot_comp = function(comp) {
  ggplot(comp, aes(x = reorder(as.character(IntDiff), IntDiff), y = Intensity, fill = Type, color = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(~round(Time, 2)) +
    theme_bw() +
    theme(legend.position = "bottom")
}

# Input data ----
hvem_input = readRDS("./Data/Real/HVEM/processed_input.RDS")
hvem_ps_m = readRDS("./Data/Real/HVEM/ps_m.RDS")
hvem_pept_cl = readRDS("./Data/Real/HVEM/peptides_cluster.RDS")
hvem_und_d = readRDS("./Data/Real/HVEM/undeuterated_dists.RDS")
# Fitted model ----
hvem_fitted = readRDS("./RawResults/hvem_fitted.RDS")

# Input data plot ----
hvem_input[ , PeptideShort := paste0(stringr::str_sub(Peptide, 1, 2), ":", stringr::str_sub(Peptide, -2, -1))]
ggplot(hvem_input, aes(x = reorder(as.character(IntDiff), IntDiff), y = Intensity)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(PeptideShort ~ round(Time, 3)) +
  scale_x_discrete(breaks = as.character(seq(0, 18, by = 3))) +
  xlab("mass shift") +
  ylab("normalized intensity") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

# Comparison of observed and predicted spectra ---- 
comp = get_comparison(hvem_fitted$x, hvem_ps_m, hvem_input, 
                      hvem_und_d, num_pars =  unique(hvem_pept_cl[order(Segment), .(Segment, MaxUptake)])[, MaxUptake + 1])
comp[, ScaledIntensity := Intensity / sum(Intensity, na.rm = T),
     by = c("Peptide", "Time", "Type")]
plot_comp(comp) +
  facet_grid(Time~Peptide)

comp[ , PeptideShort := paste0(stringr::str_sub(Peptide, 1, 2), ":", stringr::str_sub(Peptide, -2, -1))]
ggplot(comp, aes(x = reorder(as.character(IntDiff), IntDiff), y = ScaledIntensity, fill = Type, color = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(PeptideShort ~ round(Time, 3)) +
  scale_fill_manual(name = "peaks", values = viridis::viridis(3)[-1]) +
  scale_color_manual(name = "peaks", values = viridis::viridis(3)[-1]) +
  scale_x_discrete(breaks = as.character(seq(0, 18, by = 3))) +
  ylab("normalized intensity") +
  xlab("mass shift") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))


pars_to_plot_probs = hvem_fitted$x
estim_pars_by_seg = IsoHDX3:::getSegmentParametersFromBetas(pars_to_plot_probs, unique(peptides_cluster[, .(Segment, MaxUptake)])$MaxUptake + 1)
estim_seg_probs_by_time = IsoHDX3:::getSegmentProbabilitiesFromParams(estim_pars_by_seg, times)

estim_probs_by_time_dt = rbindlist(lapply(seq_along(times), function(ith_time) {
  rbindlist(lapply(seq_along(colnames(ps_m)[-1]), function(ith_seg) {
    n_ex = unique(peptides_cluster[, .(Segment, MaxUptake)])[Segment == colnames(ps_m)[-1][ith_seg], MaxUptake]
    list(Time = rep(times[ith_time], n_ex + 1),
         Segment = rep(colnames(ps_m)[-1][ith_seg], n_ex + 1),
         NumExchanged = 0:n_ex,
         Probability = estim_seg_probs_by_time[[ith_time]][[ith_seg]])
    
  }))
}))
estim_probs_by_time_dt[, Segment := factor(Segment, levels = colnames(ps_m)[-1], ordered = TRUE)]

ggplot(estim_probs_by_time_dt, aes(x = Time, y = Probability, 
                                   group = NumExchanged, color = reorder(as.character(NumExchanged), NumExchanged))) +
  geom_line(size = 2) +
  scale_color_discrete(name = "no. exchanged", palette = "viridis") +
  facet_grid(NumExchanged~Segment) +
  xlab("time") +
  ylab("mass shift probability") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))
