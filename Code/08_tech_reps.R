# Libraries ----
library(data.table)
library(ggplot2)

# Functions ---- 
get_ols_info = function(results_list, true_params) {
  model_variances_homo_ols = rbindlist(lapply(results_list, function(x) {
    if (!is.null(x)) {
      if (!is.null(x$fitted_ols)) {
        if (max(abs(x$fitted_ols$fvec)) < 1e-6) {
          list(param_id = seq_along(true_params),
               true_params = true_params,
               estimated = x$fitted_ols$x,
               iter = x$sim_iter,
               error_sd = x$error_sd)
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
}

get_true_info = function(true_params_struct) {
  true_info = rbindlist(lapply(seq_along(true_params_struct),
                               function(i) {
                                 x = true_params_struct[[i]]
                                 list(segment_id = i,
                                      true_value = x,
                                      prob = c("int", 1:(length(x) - 1)))
                               }))
  true_info[, param_id := seq_along(true_value)]
  true_info
}

get_probs = function(pars, cl_data, ps_m, seg_dt, times) {
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

get_ols_probs_table = function(results_list, true_probs_table,
                               cl_data, ps_m, seg_dt, times) {
  ex_probs_tbl = rbindlist(lapply(results_list, function(x) {
    if (!is.null(x$fitted_ols)) {
      if (max(abs(x$fitted_ols$fvec)) < 1e-6) {
        probs_comp = get_probs(x$fitted_ols$x, cl_data, ps_m, seg_dt, times)
      } else {
        probs_comp = data.table()
      }
      if (nrow(probs_comp) > 0) {
        probs_comp[, error_sd := x$error_sd]
        probs_comp[, iter := x$sim_iter]
        probs_comp
      } else {
        NULL
      }
    } else {
      NULL
    }
  }))
  comp = merge(ex_probs_tbl, true_probs_table,
               by = c("Segment", "Time", "NumExchanged"))
  comp[, Segment := factor(Segment, levels = seg_dt$SegmentSigned, ordered = T)]
  
  comp[, SegmentShort := stringr::str_replace_all(Segment, "[A-Z]+", "")]
  segm_order_short = unique(comp[order(Segment), .(Segment, SegmentShort)])$SegmentShort
  comp[, SegmentShort := factor(SegmentShort, levels = segm_order_short, ordered = T)]
  comp
}

get_averages_table = function(ex_probs_tbl) {
  averages = ex_probs_tbl[, .(mean_prob = mean(Probability.x),
                              sd = sd(Probability.x),
                              n = .N,
                              true_value = unique(Probability.y)),
                          by = c("Segment", "SegmentShort", "Time", "NumExchanged", "error_sd")]
  averages[, lower_conf := mean_prob - qt(1 - 0.05/2, n - 1) * sd / sqrt(n)]
  averages[, upper_conf := mean_prob + qt(1 - 0.05/2, n - 1) * sd / sqrt(n)]
  averages
}

plot_averages = function(ave_ex_probs_tbl) {
  ggplot(ave_ex_probs_tbl) +
    # geom_point(aes(x = Time, y = mean_prob, color = error_sd)) + # ,group = error_sd, color = error_sd
    geom_line(aes(x = Time, y = mean_prob, 
                  group = as.character(error_sd), 
                  color = as.character(error_sd)),
              linewidth = 1.5) + # , group = error_sd, color = error_sd
    # geom_linerange(aes(x = Time, ymin = lower_conf, ymax = upper_conf)) +
    geom_point(aes(x = Time, y = true_value), color = "red",
               size = 1.5) +
    facet_grid(SegmentShort ~ paste("no. ex.:", NumExchanged), scales = "free") +
    scale_color_manual(name = "std. dev.", values = viridis::viridis(3)[-1]) +
    xlab("time") +
    ylab("exchange probability") +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 8),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)) 
}

plot_boxplots = function(ex_probs_table) {
  ggplot(ex_probs_table,
         aes(x = reorder(as.character(Time), Time), y = Probability.x, 
             fill = reorder(as.character(error_sd), error_sd))) +
    geom_boxplot(size = 1.1) +
    geom_point(aes(x = reorder(as.character(Time), Time), y = Probability.y),
               color = "red", size = 1.2, inherit.aes = F) +
    scale_fill_manual(name = "std. dev", values = viridis::viridis(3)[-1]) +
    xlab("time") +
    ylab("exchange probability") +
    facet_wrap(SegmentShort ~ paste("no. ex.:", NumExchanged),
               scales = "free_y") +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
}


# Simulated data ----
dataset1_cl_data = readRDS("./Data/Simulated/DataSet1/cl_data.RDS")
dataset1_times = readRDS("./Data/Simulated/DataSet1/times.RDS")
dataset1_seg_dt = readRDS("./Data/Simulated/DataSet1/seg_dt.RDS")
dataset1_ps_m = readRDS("./Data/Simulated/DataSet1/ps_m.RDS")
dataset1_probs_by_time_dt = readRDS("./Data/Simulated/DataSet1/probs_by_time_dt.RDS")
dataset1_true_params = readRDS("./Data/Simulated/DataSet1/true_params.RDS")
dataset1_true_params_struct = readRDS("./Data/Simulated/DataSet1/true_params_struct.RDS")

# Results of simulations ---- 

res_1rep = readRDS("./RawResults/dataset1_homoscedastic.RDS")
res_2rep = readRDS("./RawResults/dataset1_2rep_homoscedastic.RDS")
res_2rep = lapply(res_2rep, function(x) {
  names(x)[c(2, 3)] = c("fitted_ols", "sim_iter")
  x
})
res_5rep = readRDS("./RawResults/dataset1_5rep_homoscedastic.RDS")
res_5rep = lapply(res_5rep, function(x) {
  names(x)[c(2, 3)] = c("fitted_ols", "sim_iter")
  x
})

# Exchange probabilities in every repetition ----

dataset1_1rep_homo_probs = get_ols_probs_table(res_1rep, 
                                               dataset1_probs_by_time_dt,
                                               dataset1_cl_data, dataset1_ps_m, 
                                               dataset1_seg_dt, dataset1_times)
dataset1_1rep_homo_probs[, NumRep := "1"]
dataset1_2rep_homo_probs = get_ols_probs_table(res_2rep, 
                                               dataset1_probs_by_time_dt,
                                               dataset1_cl_data, dataset1_ps_m, 
                                               dataset1_seg_dt, dataset1_times)
dataset1_2rep_homo_probs[, NumRep := "2"]
dataset1_5rep_homo_probs = get_ols_probs_table(res_5rep, 
                                               dataset1_probs_by_time_dt,
                                               dataset1_cl_data, dataset1_ps_m, 
                                               dataset1_seg_dt, dataset1_times)
dataset1_5rep_homo_probs[, NumRep := "5"]
dataset1_all_reps_homo_probs = rbind(dataset1_1rep_homo_probs,
                                     dataset1_2rep_homo_probs,
                                     dataset1_5rep_homo_probs)
dataset1_all_reps_homo_probs = merge(dataset1_all_reps_homo_probs,
                                     dataset1_probs_by_time_dt,
                                     by = c("Time", "Segment", "NumExchanged"),
                                     all.x = T, all.y = T)

# Overall exchange distributions 
dataset1_homo_1rep_ave_probs = get_averages_table(dataset1_1rep_homo_probs)
dataset1_homo_1rep_ave_probs[, NumRep := "1"]
dataset1_homo_2rep_ave_probs = get_averages_table(dataset1_2rep_homo_probs)
dataset1_homo_2rep_ave_probs[, NumRep := "2"]
dataset1_homo_5rep_ave_probs = get_averages_table(dataset1_5rep_homo_probs)
dataset1_homo_5rep_ave_probs[, NumRep := "5"]

dataset1_homo_all_reps_ave_probs = rbind(dataset1_homo_1rep_ave_probs,
                                         dataset1_homo_2rep_ave_probs,
                                         dataset1_homo_5rep_ave_probs)

byrep_comp_mean = dataset1_all_reps_homo_probs[, .(AveProb = mean(Probability.x),
                                                   NumConv = .N,
                                                   ProbSD = sd(Probability.x)),
                                               by = c("error_sd", "NumRep", "Segment", "Time", "NumExchanged")][order(error_sd, Segment, Time, NumExchanged, NumRep)]

byrep_comp_mean = dataset1_all_reps_homo_probs[, .(AveProb = mean(Probability.x),
                                                   TrueProb = unique(Probability.y),
                                                   NumConv = .N,
                                                   ProbSD = sd(Probability.x)),
                                               by = c("error_sd", "NumRep", "Segment", "Time", "NumExchanged")][order(error_sd, Segment, Time, NumExchanged, NumRep)]

dcast(byrep_comp_mean[error_sd == 0.5, .(NumRep, Segment, Time, NumExchanged, 
                                         TrueProb, AveProb, ProbSD)],
      Segment + Time + NumExchanged + TrueProb ~ NumRep, value.var = c("AveProb", "ProbSD"))


byrep_comp_mean[error_sd == 0.5, .(NumRep, Segment, Time, NumExchanged, 
                                   TrueProb, AveProb, ProbSD)]

byrep_comp_chisq = dataset1_all_reps_homo_probs[, .(ChiSqPerIter = sum((Probability.x - Probability.y)^2 / Probability.y)),
                                                by = c("error_sd", "NumRep", "Segment", "Time", "iter")]
byrep_comp_chisq_mean = byrep_comp_chisq[, .(AveChiSq = mean(ChiSqPerIter)), 
                                         by = c("error_sd", "NumRep", "Segment", "Time")]
dcast(byrep_comp_chisq_mean[, .(AveAveChiSq = mean(AveChiSq)), by = c("error_sd", "Segment", "NumRep")][order(error_sd, Segment, NumRep)],
      error_sd + Segment ~ NumRep, value.var = "AveAveChiSq")

byrep_comp_chisq_mean_time = byrep_comp_chisq[, .(AveChiSq = mean(ChiSqPerIter)), 
                                              by = c("error_sd", "NumRep", "Segment")]

byrep_comp_chisq_med = byrep_comp_chisq[, .(AveChiSq = median(ChiSqPerIter)), 
                                        by = c("error_sd", "NumRep", "Segment", "Time")]

ggplot(byrep_comp_chisq_med,
       aes(x = reorder(as.character(Time), Time), y = AveChiSq, 
           color = reorder(as.character(NumRep), NumRep))) +
  geom_point(size = 1.5) +
  scale_color_manual(name = "no. tech. rep.",
                     values = viridis::viridis(4)[-1]) +
  facet_grid(as.numeric(Segment) ~ error_sd) +
  theme_bw() +
  xlab("time") +
  ylab("average chi-squared distance") +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
ggsave("Figures/byrep_chisq_med.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)


# ggplot(byrep_comp_mean, aes(x = reorder(as.character(Time), Time),
#                             color = reorder(NumRep, as.numeric(NumRep)),
#                             y = AveProb)) +
#   geom_line() +
#   facet_wrap(Segment ~ NumExchanged, scales = "free") +
#   scale_color_manual(name = "no. tech. rep.", values = viridis::viridis(4)[-1]) +
#   theme_bw() +
#   ylab("exchange probability") +
#   xlab("time") +
#   theme(legend.position = "bottom")


ggplot(dataset1_all_reps_homo_probs, aes(x = reorder(as.character(Time), Time),
                                         color = reorder(NumRep, as.numeric(NumRep)),
                                         y = Probability.x)) +
  geom_boxplot() +
  facet_grid(Segment ~ NumExchanged) +
  scale_color_manual(name = "no. tech. rep.", values = viridis::viridis(4)[-1]) +
  theme_bw() +
  ylab("exchange probability") +
  xlab("time") +
  theme(legend.position = "bottom")
ggsave("Figures/byrep_probs_all_segments.png", device = "png", scale = 1, width = 10, height = 10, units = "in", dpi = 300)

ggplot(dataset1_all_reps_homo_probs[Segment == unique(Segment)[1]], 
       aes(x = reorder(as.character(Time), Time),
           color = reorder(NumRep, as.numeric(NumRep)),
           y = Probability.x)) +
  geom_boxplot() +
  geom_point(aes(x = reorder(as.character(Time), Time),
                 y = Probability.y),
             inherit.aes = F, color = "red", size = 1) +
  facet_grid(Segment ~ NumExchanged) +
  scale_color_manual(name = "no. tech. rep.", values = viridis::viridis(4)[-1]) +
  theme_bw() +
  ylab("exchange probability") +
  xlab("time") +
  theme(legend.position = "bottom")
ggsave("Figures/byrep_probs_single_segment.png", device = "png", scale = 1, width = 10, height = 10, units = "in", dpi = 300)

# Bias and variance
param_info_1rep_homo = get_ols_info(res_1rep, dataset1_true_params)
param_info_1rep_homo[, NumRep := "1"]
param_info_2rep_homo = get_ols_info(res_2rep, dataset1_true_params)
param_info_2rep_homo[, NumRep := "2"]
param_info_5rep_homo = get_ols_info(res_5rep, dataset1_true_params)
param_info_5rep_homo[, NumRep := "5"]
param_info_all_homo = rbind(param_info_1rep_homo, param_info_2rep_homo,
                            param_info_5rep_homo)

true_info = get_true_info(dataset1_true_params_struct)

ggplot(merge(param_info_all_homo,
             true_info, by = "param_id")[prob != "int", .(MeanRelBias = mean((estimated - true_params) / true_params)),
                                         by = c("param_id", "error_sd", "NumRep")][MeanRelBias < 100],
       aes(x = param_id, y = MeanRelBias, color = NumRep)) +
  geom_point(size = 3) +
  facet_grid(~error_sd) +
  scale_color_manual(name = "no. tech. replicates", values = viridis::viridis(4)[-1]) +
  xlab("parameter ID") +
  ylab("relative bias") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
ggsave("Figures/byrep_relbias.png", device = "png", 
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)

ggplot(merge(param_info_all_homo, true_info, by = "param_id")[, .(variance = var(estimated, na.rm = T)),
                                                              by = c("param_id", "error_sd", "NumRep")][variance < 100],
       aes(x = param_id, y = variance, color = NumRep)) +
  geom_point(size = 3) +
  facet_grid(~error_sd) +
  scale_color_manual(name = "no. tech. replicates", values = viridis::viridis(4)[-1]) +
  xlab("parameter ID") +
  ylab("empirical variance") +
  theme_bw() +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12)) 
ggsave("Figures/byrep_variance.png", device = "png", 
       scale = 1, width = 10, height = 10, units = "in", dpi = 300)

