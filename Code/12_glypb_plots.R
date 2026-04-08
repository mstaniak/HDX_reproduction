# Libraries ----
library(data.table)
library(ggplot2)
#Functions ---
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
# Input data ----
glypb_input = readRDS("./Data/Real/GlyPb/processed_input.RDS")
glypb_fitted = readRDS("./RawResults/glypb_fitted.RDS")
num_parameters = readRDS("./Data/Real/GlyPb/num_parameters.RDS")
t0_info = readRDS("./Data/Real/GlyPb/t0_data.RDS")
und_d = readRDS("./Data/Real/GlyPb/und_d.RDS")
pept_cl = readRDS("./Data/Real/GlyPb/pept_cl.RDS")
ps_m = readRDS("./Data/Real/GlyPb/ps_m.RDS")
# Comparison of predicted and observed spectra ----
comp_g1 = get_comparison(glypb_fitted, ps_m,
                         glypb_input, und_d, num_parameters + 1)
plot_comp(comp_g1[, .(Peptide, Time, IntDiff, Intensity, Type = ifelse(Type == "Fitted", "fitted", "observed"))]) +
  facet_grid(paste(stringr::str_sub(Peptide, 1, 3), stringr::str_sub(Peptide, -3, -1), sep = ".")~Time, scales = "free_y") +
  scale_fill_manual(name = "peaks", values = viridis::viridis(3)[-1]) +
  scale_color_manual(name = "peaks", values = viridis::viridis(3)[-1]) +
  xlab("mass shift") +
  ylab("intensity") +
  scale_x_discrete(breaks = reorder(as.character(seq(1, 26, by = 4)), seq(1, 26, by = 4))) +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        strip.text = element_text(size = 10)) 
ggsave("./Figures/msHDX_spectra_comp.png", device = "png", width = 10, height = 5, units = "in", scale = 1, dpi = 300)
# Fitted segment-level probabilities ----
pars_by_seg_2 = IsoHDX3:::getSegmentParametersFromBetas(glypb_fitted$par, num_parameters + 1)
probs_by_seg_2 = IsoHDX3:::getSegmentProbabilitiesFromParams(pars_by_seg_2, unique(glypb_input$Time))

segments = as.character(unique(peptides_cluster[, .(Segment, MaxUptake)])$Segment)
segm_limits = lapply(stringr::str_split(segments, ","), function(x) stringr::str_extract_all(x, "[0-9]+", simplify = T)[, 1])
segm_order = rbindlist(lapply(seq_along(segm_limits), function(i) {
  list(Segment = segments[i],
       StartSeg = segm_limits[[i]][1],
       EndSeg = segm_limits[[i]][2])
}))[order(StartSeg, EndSeg), Segment]
seg_probs = rbindlist(lapply(seq_along(probs_by_seg_2), function(ith_time) {
  rbindlist(lapply(seq_along(probs_by_seg_2[[ith_time]]), function(ith_segment) {
    list(Segment = segments[ith_segment],
         NumExchanged = 0:(length(probs_by_seg_2[[ith_time]][[ith_segment]]) - 1),
         Probability = probs_by_seg_2[[ith_time]][[ith_segment]],
         Time = unique(glypb_input$Time)[ith_time])
  }))
}))

seg_probs[, Segment := factor(Segment, levels = segm_order, ordered = T)]
ggplot(seg_probs, aes(x = Time, y = Probability, 
                      group = NumExchanged, color = reorder(as.character(NumExchanged), NumExchanged))) +
  geom_line(size = 1.5) +
  geom_point() +
  facet_wrap(Segment~NumExchanged) +
  scale_color_manual(name = "no. exchanged", values = viridis::viridis(9)[-1]) +
  theme_bw() +
  xlab("time [h]") +
  ylab("probability") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 8),
        legend.title = element_text(size = 14),
        strip.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "bottom")
ggsave("Figures/msHDX_segm_probs.png", device = "png", width = 10, height = 10, units = "in", scale = 1, dpi = 300)

