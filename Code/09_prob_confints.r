library(data.table)
library(ggplot2)

dataset1_res_hetero = readRDS("./RawResults/dataset1_heteroscedastic.RDS")
dataset1_res_homo = readRDS("./RawResults/dataset1_homoscedastic.RDS")
dataset1_num_parameters = sapply(readRDS("./Data/Simulated/DataSet1/true_params_struct.RDS"), length)
dataset1_num_parameters = dataset1_num_parameters - 1
dataset1_ps_m = readRDS("./Data/Simulated/DataSet1/ps_m.RDS")
dataset1_cl_dt = readRDS("./Data/Simulated/DataSet1/cl_data.RDS")
dataset1_seg_dt = readRDS("./Data/Simulated/DataSet1/seg_dt.RDS")
dataset1_probs_by_time_dt = readRDS("./Data/Simulated/DataSet1/probs_by_time_dt.RDS")
dataset1_times = readRDS("./Data/Simulated/DataSet1/times.RDS")

get_jacobian_probs = function(segment_id, num_parameters, time, parameters) {
  num_parameters_segment = (num_parameters + 1)[segment_id]
  segment_params = IsoHDX3:::getSegmentParametersFromBetas(parameters, num_parameters + 1)
  
  derivatives = lapply(seq_len(num_parameters_segment), function(parameter_id) {
    sapply(seq_len(num_parameters_segment) - 1, function(num_exchanged) {
      get_derivative(segment_params[[segment_id]], time = time, parameter_id = parameter_id, num_exchanged = num_exchanged,
                     max_exchanged = num_parameters[segment_id])
    })
  })
  do.call("cbind", derivatives)
}

get_derivative = function(betas_segment, time, parameter_id, num_exchanged, max_exchanged) {
  betas = betas_segment
  betas[length(betas)] = exp(betas[length(betas)])
  intercept = betas[1]
  betas = betas[-1]
  total = 1 + sum(exp(betas * time + intercept))
  if (parameter_id == 1) {
    if (num_exchanged == 0) {
      -(total - 1) / (total^2)
    } else {
      exp(betas[num_exchanged] * time + intercept) / (total^2)
    }
  } else if (parameter_id > 1 & parameter_id < max_exchanged + 1) {
    if (num_exchanged == 0) { # parameter_id > 1 <=> num_exchanged + 1 != parameter_id
      - time * exp(betas[parameter_id - 1] * time + intercept) / (total^2)
    } else if (((num_exchanged + 1) == parameter_id) & parameter_id > 1) {
      current_exp = exp(betas[num_exchanged] * time + intercept)
      (total - current_exp) * time * current_exp / (total ^ 2)
    } else if (((num_exchanged + 1) != parameter_id) & parameter_id > 1) {
      current_exp = exp(betas[parameter_id - 1] * time + intercept)
      exchanged_exp = exp(betas[num_exchanged] * time + intercept)
      -time * exchanged_exp * current_exp / (total ^ 2)
    }
  } else { # parameter_id == max_exchanged + 1
    if (num_exchanged == 0) {
      - time * exp(betas[length(betas)] * time + betas_segment[length(betas_segment)] + intercept) / (total^2)
    } else if (((num_exchanged + 1) == parameter_id) & parameter_id > 1) {
      current_exp = exp(betas[num_exchanged] * time + intercept)
      ((total - current_exp) * time * current_exp * exp(betas_segment[length(betas_segment)])) / (total ^ 2)
    } else if (((num_exchanged + 1) != parameter_id) & parameter_id > 1) {
      current_exp = exp(betas[parameter_id - 1] * time + intercept)
      exchanged_exp = exp(betas[num_exchanged] * time + intercept)
      -time * exchanged_exp * current_exp * exp(betas_segment[length(betas_segment)]) / (total ^ 2)
    }
  }
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


get_plgls_probs_table = function(results_list, true_probs_table,
                                 cl_data, ps_m, seg_dt, times) {
  ex_probs_tbl = rbindlist(lapply(results_list, function(x) {
    if (!is.null(x$fitted_plgls)) {
      if (max(abs(x$fitted_plgls$fvec)) < 1e-6) {
        probs_comp = get_probs(x$fitted_plgls$x, cl_data, ps_m, seg_dt, times)
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

dataset1_ps_m

all_probs_tb = get_plgls_probs_table(dataset1_res_hetero, 
                                     dataset1_probs_by_time_dt,
                                     dataset1_cl_dt, dataset1_ps_m, 
                                     dataset1_seg_dt, dataset1_times)
all_probs_tb_homo = get_plgls_probs_table(dataset1_res_homo, 
                                          dataset1_probs_by_time_dt,
                                          dataset1_cl_dt, dataset1_ps_m, 
                                          dataset1_seg_dt, dataset1_times)


# All segments

num_parameters = dataset1_num_parameters + 1
segments = colnames(dataset1_ps_m)[-ncol(dataset1_ps_m)]
segment_ranges = lapply(seq_along(segments), function(i) {
  if (i == 1) {
    subset_betas = 1:num_parameters[i]
  } else {
    subset_betas = (sum(num_parameters[1:(i - 1)]) + 1):(sum(num_parameters[1:i]))
  }
  subset_betas
})
segment_ranges

all_segments_hetero = lapply(
  seq_along(segments), function(ith_segment) {
    # print(ith_segment)
    lapply(c(0.01, 0.05), function(hetero_sd) {
      # print(hetero_sd)
      lapply(unique(dataset1_times), function(time) {
        # print(time)
        probs_tbl_subset = all_probs_tb[Segment == segments[ith_segment] & Time == time & error_sd == hetero_sd]
        lapply(probs_tbl_subset[, unique(iter)], function(i) {
          # print(i)
          res_this_iter = dataset1_res_hetero[sapply(dataset1_res_hetero, function(x) if (is.null(x)) FALSE else (x$sim_iter == i & x$error_sd == hetero_sd))][[1]]
          this_iter_pars = res_this_iter$fitted_plgls$x
          probs_jac = get_jacobian_probs(ith_segment, dataset1_num_parameters, time, this_iter_pars)
          plgls_cov_i1 = res_this_iter$cov
          jac_transf = probs_jac %*% plgls_cov_i1[segment_ranges[[ith_segment]], segment_ranges[[ith_segment]]] %*% t(probs_jac)
          this_iter_vars = diag(jac_transf)
          
          this_iter_info = probs_tbl_subset[iter == i][order(NumExchanged)]
          
          list(segment = segments[ith_segment],
               time = time,
               error_sd = hetero_sd,
               iter = res_this_iter$sim_iter,
               num_exchanged = this_iter_info$NumExchanged,
               # det = det(res_this_iter$cov),
               covered = this_iter_info$Probability.x - qnorm(1 - (0.05/2)) * sqrt(this_iter_vars) <= this_iter_info$Probability.y & this_iter_info$Probability.x + qnorm(1 - (0.05/2)) * sqrt(this_iter_vars) >= this_iter_info$Probability.y)
        })
        
      })
    })
  }
)
all_segments_hetero_tbl = data.table::rbindlist(unlist(unlist(unlist(all_segments_hetero, F, F), F, F), F, F))

table(all_segments_hetero_tbl$covered)
coverages_hetero = all_segments_hetero_tbl[, .(coverage = mean(covered),
                                               upper_ci = qbeta(0.05/2, sum(covered), .N - sum(covered) + 1),
                                               lower_ci = qbeta(1 - 0.05/2, sum(covered) + 1, .N - sum(covered))),
                                           # upper_ci = mean(covered) + qnorm(1 - (0.05/2)) * sqrt(mean(covered) * (1 - mean(covered)) / .N),
                                           # lower_ci = mean(covered) - qnorm(1 - (0.05/2)) * sqrt(mean(covered) * (1 - mean(covered)) / .N)),
                                           by = c("segment", "time", "error_sd", "num_exchanged")]

coverages_hetero[, segment := factor(segment, levels = segments, ordered = T)]
ggplot(coverages_hetero, 
       aes(x = reorder(as.character(time), time), y = coverage, color = reorder(as.character(num_exchanged), num_exchanged))) +
  geom_point(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_linerange(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0.95, color = "red") +
  facet_grid(segment ~ error_sd) +
  xlab("time") +
  ylab("coverage of the 95% CI") +
  coord_cartesian(ylim = c(0.6, 1.05)) +
  theme_bw() +
  scale_color_discrete(name = "no. exchanged", palette = "viridis") +
  # xlim(c(0, 0.1)) +
  # ylim(c(0, 0.1)) +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 10))
ggsave("./Figures/prob_coverage_hetero.png", device = "png", scale = 1,
       width = 10, height = 10, units = "in", dpi = 300)


all_segments_homo = lapply(
  seq_along(segments), function(ith_segment) {
    # print(ith_segment)
    lapply(c(0.5, 1), function(hetero_sd) {
      # print(hetero_sd)
      lapply(unique(dataset1_times), function(time) {
        # print(time)
        probs_tbl_subset = all_probs_tb_homo[Segment == segments[ith_segment] & Time == time & error_sd == hetero_sd]
        lapply(probs_tbl_subset[, unique(iter)], function(i) {
          # print(i)
          res_this_iter = dataset1_res_homo[sapply(dataset1_res_homo, function(x) if (is.null(x)) FALSE else (x$sim_iter == i & x$error_sd == hetero_sd))][[1]]
          this_iter_pars = res_this_iter$fitted_plgls$x
          probs_jac = get_jacobian_probs(ith_segment, dataset1_num_parameters, time, this_iter_pars)
          plgls_cov_i1 = res_this_iter$cov
          jac_transf = probs_jac %*% plgls_cov_i1[segment_ranges[[ith_segment]], segment_ranges[[ith_segment]]] %*% t(probs_jac)
          this_iter_vars = diag(jac_transf)
          
          this_iter_info = probs_tbl_subset[iter == i][order(NumExchanged)]
          
          list(segment = segments[ith_segment],
               time = time,
               error_sd = hetero_sd,
               iter = res_this_iter$sim_iter,
               num_exchanged = this_iter_info$NumExchanged,
               # det = det(res_this_iter$cov),
               covered = this_iter_info$Probability.x - qnorm(1 - (0.05/2)) * sqrt(this_iter_vars) <= this_iter_info$Probability.y & this_iter_info$Probability.x + qnorm(1 - (0.05/2)) * sqrt(this_iter_vars) >= this_iter_info$Probability.y)
        })
        
      })
    })
  }
)
all_segments_homo_tbl = data.table::rbindlist(unlist(unlist(unlist(all_segments_homo, F, F), F, F), F, F))

# table(all_segments_hetero_tbl$covered)

coverages_homo = all_segments_homo_tbl[, .(coverage = mean(covered),
                                           upper_ci = qbeta(0.05/2, sum(covered), .N - sum(covered) + 1),
                                           lower_ci = qbeta(1 - 0.05/2, sum(covered) + 1, .N - sum(covered))),
                                           # upper_ci = mean(covered) + qnorm(1 - (0.05/2)) * sqrt(mean(covered) * (1 - mean(covered)) / .N),
                                           # lower_ci = mean(covered) - qnorm(1 - (0.05/2)) * sqrt(mean(covered) * (1 - mean(covered)) / .N)),
                                       by = c("segment", "time", "error_sd", "num_exchanged")]

coverages_homo[, segment := factor(segment, levels = segments, ordered = T)]
coverages_homo

ggplot(coverages_homo, 
       aes(x = reorder(as.character(time), time), y = coverage, color = reorder(as.character(num_exchanged), num_exchanged))) +
  geom_point(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_linerange(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0.95, color = "red") +
  facet_grid(segment ~ error_sd) +
  xlab("time") +
  ylab("coverage of the 95% CI") +
  coord_cartesian(ylim = c(0.6, 1.05)) +
  theme_bw() +
  scale_color_discrete(name = "no. exchanged", palette = "viridis") +
  # xlim(c(0, 0.1)) +
  # ylim(c(0, 0.1)) +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 10))
ggsave("./Figures/prob_coverage_homo.png", device = "png", scale = 1,
       width = 10, height = 10, units = "in", dpi = 300)

prob_empirical_variances_hetero = all_probs_tb[, .(var = var(Probability.x, na.rm = T)),
                                               by = c("Segment", "Time", "NumExchanged", "error_sd")]
prob_empirical_variances_homo = all_probs_tb_homo[, .(var = var(Probability.x, na.rm = T)),
                                                  by = c("Segment", "Time", "NumExchanged", "error_sd")]

model_variances_raw_hetero = lapply(
  seq_along(segments), function(ith_segment) {
    # print(ith_segment)
    lapply(c(0.01, 0.05), function(hetero_sd) {
      # print(hetero_sd)
      lapply(unique(dataset1_times), function(time) {
        # print(time)
        probs_tbl_subset = all_probs_tb[Segment == segments[ith_segment] & Time == time & error_sd == hetero_sd]
        lapply(probs_tbl_subset[, unique(iter)], function(i) {
          # print(i)
          res_this_iter = dataset1_res_hetero[sapply(dataset1_res_hetero, function(x) if (is.null(x)) FALSE else (x$sim_iter == i & x$error_sd == hetero_sd))][[1]]
          this_iter_pars = res_this_iter$fitted_plgls$x
          probs_jac = get_jacobian_probs(ith_segment, dataset1_num_parameters, time, this_iter_pars)
          plgls_cov_i1 = res_this_iter$cov
          jac_transf = probs_jac %*% plgls_cov_i1[segment_ranges[[ith_segment]], segment_ranges[[ith_segment]]] %*% t(probs_jac)
          this_iter_vars = diag(jac_transf)
          
          this_iter_info = probs_tbl_subset[iter == i][order(NumExchanged)]
          
          list(segment = segments[ith_segment],
               time = time,
               error_sd = hetero_sd,
               iter = res_this_iter$sim_iter,
               num_exchanged = this_iter_info$NumExchanged,
               variance = this_iter_vars)
        })
        
      })
    })
  }
)
model_variances_hetero = data.table::rbindlist(unlist(unlist(unlist(model_variances_raw_hetero, F, F), F, F), F, F))


model_variances_homo_raw = lapply(
  seq_along(segments), function(ith_segment) {
    # print(ith_segment)
    lapply(c(0.5, 1), function(hetero_sd) {
      # print(hetero_sd)
      lapply(unique(dataset1_times), function(time) {
        # print(time)
        probs_tbl_subset = all_probs_tb_homo[Segment == segments[ith_segment] & Time == time & error_sd == hetero_sd]
        lapply(probs_tbl_subset[, unique(iter)], function(i) {
          # print(i)
          res_this_iter = dataset1_res_homo[sapply(dataset1_res_homo, function(x) if (is.null(x)) FALSE else (x$sim_iter == i & x$error_sd == hetero_sd))][[1]]
          this_iter_pars = res_this_iter$fitted_plgls$x
          probs_jac = get_jacobian_probs(ith_segment, dataset1_num_parameters, time, this_iter_pars)
          plgls_cov_i1 = res_this_iter$cov
          jac_transf = probs_jac %*% plgls_cov_i1[segment_ranges[[ith_segment]], segment_ranges[[ith_segment]]] %*% t(probs_jac)
          this_iter_vars = diag(jac_transf)
          
          this_iter_info = probs_tbl_subset[iter == i][order(NumExchanged)]
          
          list(segment = segments[ith_segment],
               time = time,
               error_sd = hetero_sd,
               iter = res_this_iter$sim_iter,
               num_exchanged = this_iter_info$NumExchanged,
               variance = this_iter_vars)
        })
        
      })
    })
  }
)
model_variances_homo = data.table::rbindlist(unlist(unlist(unlist(model_variances_homo_raw, F, F), F, F), F, F))
model_variances_homo_mean = model_variances_homo[, .(variance = mean(variance, na.rm = T)),
                                                 by = c("segment", "time", "error_sd", "num_exchanged")]
class(model_variances_hetero$time)


limits_tbl = data.table(
  time = c((0.1 + 0.5) / 2, (0.5 + 1) / 2, (1 + 2.5) / 2, (2.5 + 5) / 2, (5 + 10) / 2))


all_vars_homo = rbind(cbind(model_variances_homo_mean, vartype = "mean model-based"),
                      cbind(prob_empirical_variances_homo[, .(segment = Segment, time = Time,
                                                              error_sd, num_exchanged = NumExchanged, variance = var)], vartype = "empirical"))
all_vars_homo[, time := factor(as.character(time),
                               levels = sort(c(unique(time), unique(limits_tbl$time))),
                               ordered = T)]
all_vars_homo[ , unique(time)]
limits_tbl[, time := factor(as.character(time), 
                            levels = levels(all_vars_homo$time),
                            ordered = T)]

library(patchwork)
segments

by_segment_plot = lapply(segments, function(seg) {
  ggplot(all_vars_homo[segment == seg],
         aes(x = time, y = variance, 
             shape = vartype, 
             group = paste(error_sd, segment, num_exchanged, time),
             color = reorder(as.character(num_exchanged), num_exchanged))) +
    # scale_x_discrete(breaks = levels(all_vars_homo$time),
    #                  labels = levels(all_vars_homo$time)) +
    geom_vline(aes(xintercept = time),
               data = data.frame(time = c(1.5, 2.5, 3.5, 4.5, 5.5)),
               linetype = 2, color = "grey") +
    geom_point(size = 2, position = position_dodge(width = 0.9)) +
    facet_grid(segment ~ error_sd, scales = "free_y") +
    xlab("time") +
    ylab("variance") +
    theme_bw() +
    scale_color_discrete(name = "no. exchanged", palette = "viridis") +
    scale_shape_discrete(name = "variance") 
  # xlim(c(0, 0.1)) +
  # ylim(c(0, 0.1)) 
})
by_segment_plot[1:5] = lapply(by_segment_plot[1:5], function(x) {
  x + theme(legend.position = "none",
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            strip.text = element_text(size = 10))
  
})
by_segment_plot[[6]] = by_segment_plot[[6]] +     theme(legend.position = "bottom",
                                                        axis.text = element_text(size = 10),
                                                        axis.title = element_text(size = 10),
                                                        legend.text = element_text(size = 10),
                                                        legend.title = element_text(size = 10),
                                                        strip.text = element_text(size = 10))


length(by_segment_plot)
by_segment_plot[[1]] + by_segment_plot[[2]] + by_segment_plot[[3]] + by_segment_plot[[4]] + by_segment_plot[[5]] + by_segment_plot[[6]] + plot_layout(ncol = 1)
ggsave("./Figures/prob_var_emp_mod_homo.png", device = "png", scale = 1,
       width = 10, height = 10, units = "in", dpi = 300)

# Hetero error ----
model_variances_hetero_raw = lapply(
  seq_along(segments), function(ith_segment) {
    # print(ith_segment)
    lapply(c(0.01, 0.05), function(hetero_sd) {
      # print(hetero_sd)
      lapply(unique(dataset1_times), function(time) {
        # print(time)
        probs_tbl_subset = all_probs_tb[Segment == segments[ith_segment] & Time == time & error_sd == hetero_sd]
        lapply(probs_tbl_subset[, unique(iter)], function(i) {
          # print(i)
          res_this_iter = dataset1_res_hetero[sapply(dataset1_res_hetero, function(x) if (is.null(x)) FALSE else (x$sim_iter == i & x$error_sd == hetero_sd))][[1]]
          this_iter_pars = res_this_iter$fitted_plgls$x
          probs_jac = get_jacobian_probs(ith_segment, dataset1_num_parameters, time, this_iter_pars)
          plgls_cov_i1 = res_this_iter$cov
          jac_transf = probs_jac %*% plgls_cov_i1[segment_ranges[[ith_segment]], segment_ranges[[ith_segment]]] %*% t(probs_jac)
          this_iter_vars = diag(jac_transf)
          
          this_iter_info = probs_tbl_subset[iter == i][order(NumExchanged)]
          
          list(segment = segments[ith_segment],
               time = time,
               error_sd = hetero_sd,
               iter = res_this_iter$sim_iter,
               num_exchanged = this_iter_info$NumExchanged,
               variance = this_iter_vars)
        })
        
      })
    })
  }
)
model_variances_hetero = data.table::rbindlist(unlist(unlist(unlist(model_variances_hetero_raw, F, F), F, F), F, F))
model_variances_hetero_mean = model_variances_hetero[, .(variance = mean(variance, na.rm = T)),
                                                 by = c("segment", "time", "error_sd", "num_exchanged")]
limits_tbl = data.table(
  time = c((0.1 + 0.5) / 2, (0.5 + 1) / 2, (1 + 2.5) / 2, (2.5 + 5) / 2, (5 + 10) / 2))

all_vars_hetero = rbind(cbind(model_variances_hetero_mean, vartype = "mean model-based"),
                      cbind(prob_empirical_variances_hetero[, .(segment = Segment, time = Time,
                                                              error_sd, num_exchanged = NumExchanged, variance = var)], vartype = "empirical"))
all_vars_hetero[, time := factor(as.character(time),
                               levels = sort(c(unique(time), unique(limits_tbl$time))),
                               ordered = T)]
all_vars_hetero[ , unique(time)]
limits_tbl[, time := factor(as.character(time), 
                            levels = levels(all_vars_hetero$time),
                            ordered = T)]

by_segment_plot_hetero = lapply(segments, function(seg) {
  ggplot(all_vars_hetero[segment == seg],
         aes(x = time, y = variance, 
             shape = vartype, 
             group = paste(error_sd, segment, num_exchanged, time),
             color = reorder(as.character(num_exchanged), num_exchanged))) +
    # scale_x_discrete(breaks = levels(all_vars_hetero$time),
    #                  labels = levels(all_vars_hetero$time)) +
    geom_vline(aes(xintercept = time),
               data = data.frame(time = c(1.5, 2.5, 3.5, 4.5, 5.5)),
               linetype = 2, color = "grey") +
    geom_point(size = 2, position = position_dodge(width = 0.9)) +
    facet_grid(segment ~ error_sd, scales = "free_y") +
    xlab("time") +
    ylab("variance") +
    theme_bw() +
    scale_color_discrete(name = "no. exchanged", palette = "viridis") +
    scale_shape_discrete(name = "variance") 
  # xlim(c(0, 0.1)) +
  # ylim(c(0, 0.1)) 
})
by_segment_plot_hetero[1:5] = lapply(by_segment_plot_hetero[1:5], function(x) {
  x + theme(legend.position = "none",
            axis.text = element_text(size = 10),
            axis.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10),
            strip.text = element_text(size = 10))
  
})
by_segment_plot_hetero[[6]] = by_segment_plot_hetero[[6]] +     theme(legend.position = "bottom",
                                                        axis.text = element_text(size = 10),
                                                        axis.title = element_text(size = 10),
                                                        legend.text = element_text(size = 10),
                                                        legend.title = element_text(size = 10),
                                                        strip.text = element_text(size = 10))


by_segment_plot_hetero[[1]] + by_segment_plot_hetero[[2]] + by_segment_plot_hetero[[3]] + by_segment_plot_hetero[[4]] + by_segment_plot_hetero[[5]] + by_segment_plot_hetero[[6]] + plot_layout(ncol = 1)
ggsave("./Figures/prob_var_emp_mod_hetero.png", device = "png", scale = 1,
       width = 10, height = 10, units = "in", dpi = 300)

