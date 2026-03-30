# Libraries ---- 
library(data.table)
library(ggplot2)

# Functions ----
# TODO: OLS for all datasets based on the code below:
# num_hessians_homo_error = readRDS("num_hessians_homo_error_smallcl.RDS")
# inverse_hessians = lapply(num_hessians_homo_error, function(x) {
#   if (!is.null(x)) {
#     x$inv_hess = solve(x$num_hessian) 
#   } else {
#     x$inv_hess = NULL
#   }
#   x
# })
# 
# sigma_sq_estim = lapply(res_raw, function(x) {
#   if (!is.null(x)) {
#     if (!is.null(x$fitted_ols)) {
#       if (max(abs(x$fitted_ols$fvec)) < 1e-6) {
#         predicted = IsoHDX3:::getExpectedSpectra(x$fitted_ols$x, ps_m, num_parameters + 1, x$data, und_d)
#         comp = merge(x$data, predicted, by = c("Rep", "Peptide", "Time", "IntDiff"))
#         x$comp = comp
#         x$sigma_sq = (1 / (nrow(x$data) - length(true_params))) * comp[, sum((Intensity - ExpectedPeak)^2)]
#         return(x)
#       } else {
#         NULL
#       }
#     } else {
#       NULL
#     }   
#   } else {
#     NULL
#   }
# })
# 
# varcov_estimates_data = lapply(seq_along(inverse_hessians), function(i) {
#   x = sigma_sq_estim[[i]]
#   x$inv_hess = inverse_hessians[[i]]$inv_hess
#   x
# })
# varcov_estimates = lapply(varcov_estimates_data, function(x) {
#   x$varcov = x$sigma_sq * x$inv_hess
#   x
# })
# varcov_estimates_data[[1]]
# varcov_estimates[[1]]
# 
# model_variances_homo = rbindlist(lapply(varcov_estimates, function(x) {
#   if (!is.null(x)) {
#     if (!is.null(x$fitted_ols)) {
#       if (max(abs(x$fitted_ols$fvec)) < 1e-6) {
#         list(param_id = seq_along(true_params),
#              true_params = true_params,
#              estimated = x$fitted_ols$x,
#              model_variance = diag(x$varcov),
#              iter = x$sim_iter,
#              error_sd = x$error_sd)
#       } else {
#         NULL
#       }
#     } else {
#       NULL
#     }
#   } else {
#     NULL
#   }
# }))
# 
# comp_all_homo_ols = model_variances_homo[, .(mean_model_var = mean(model_variance),
#                                              empirical_var = var(estimated),
#                                              true_value = unique(true_params),
#                                              mean_estimated = mean(estimated),
#                                              n = .N),
#                                          by = c("param_id", "error_sd")]

get_plgls_model_variances = function(results_list, true_params) {
  rbindlist(lapply(results_list, function(x) {
    if (!is.null(x)) {
      if (!is.null(x$fitted_plgls)) {
        if (max(abs(x$fitted_plgls$fvec)) < 1e-6) {
          list(param_id = seq_along(true_params),
               true_params = true_params,
               estimated = x$fitted_plgls$x,
               model_variance = diag(x$cov),  
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

get_coverages = function(model_variances, true_params_struct) {
  true_info = get_true_info(true_params_struct)
  
  model_variances_with_true = merge(model_variances,
                                    true_info,
                                    by = c("param_id"))
  model_variances_with_true[, lower := estimated - qnorm(1 - 0.05/2) * sqrt(model_variance)]
  model_variances_with_true[, upper := estimated + qnorm(1 - 0.05/2) * sqrt(model_variance)]
  
  
  coverages = model_variances_with_true[, .(coverage = 100 * mean(lower <= true_params & true_params <= upper),
                                            se = 100 * sqrt(mean(lower <= true_params & true_params <= upper) * (1 - mean(lower <= true_params & true_params <= upper)) / .N)),
                                        by = c("param_id", "error_sd", "segment_id")][order(error_sd, segment_id, param_id), .(error_sd, segment_id, param_id, coverage, se)]
  coverages
}


# Simulated data and results ----
dataset1_true_params_struct = readRDS("./Data/Simulated/DataSet1/true_params_struct.RDS")
dataset1_true_params = readRDS("./Data/Simulated/DataSet1/true_params.RDS")

dataset1_hetero_res = readRDS("./RawResults/dataset1_heteroscedastic.RDS")
dataset1_homo_res = readRDS("./RawResults/dataset1_homoscedastic.RDS")

dataset2_true_params_struct = readRDS("./Data/Simulated/DataSet2/true_params_struct.RDS")
dataset2_true_params = readRDS("./Data/Simulated/DataSet2/true_params.RDS")

dataset2_hetero_res = readRDS("./RawResults/dataset2_heteroscedastic.RDS")
dataset2_homo_res = readRDS("./RawResults/dataset2_homoscedastic.RDS")

dataset3_true_params_struct = readRDS("./Data/Simulated/DataSet3/true_params_struct.RDS")
dataset3_true_params = readRDS("./Data/Simulated/DataSet3/true_params.RDS")

dataset3_hetero_res = readRDS("./RawResults/dataset3_heteroscedastic.RDS")
dataset3_homo_res = readRDS("./RawResults/dataset3_homoscedastic.RDS")

# Model variances in each repetition ----

dataset1_plgls_model_variances_hetero = get_plgls_model_variances(dataset1_hetero_res,
                                                                  dataset1_true_params)
dataset1_plgls_model_variances_homo = get_plgls_model_variances(dataset1_homo_res,
                                                                dataset1_true_params)


dataset2_plgls_model_variances_hetero = get_plgls_model_variances(dataset2_hetero_res,
                                                                  dataset2_true_params)
dataset2_plgls_model_variances_homo = get_plgls_model_variances(dataset2_homo_res,
                                                                dataset2_true_params)

dataset3_plgls_model_variances_hetero = get_plgls_model_variances(dataset3_hetero_res,
                                                                  dataset3_true_params)
dataset3_plgls_model_variances_homo = get_plgls_model_variances(dataset3_homo_res,
                                                                dataset3_true_params)

# Empirical coverages for model parameters based on PL-GLS ----
dataset1_plgls_coverage_hetero = get_coverages(dataset1_plgls_model_variances_hetero, 
                                               dataset1_true_params_struct)
dataset1_plgls_coverage_homo = get_coverages(dataset1_plgls_model_variances_homo, 
                                               dataset1_true_params_struct)

dataset2_plgls_coverage_hetero = get_coverages(dataset2_plgls_model_variances_hetero, 
                                               dataset2_true_params_struct)
dataset2_plgls_coverage_homo = get_coverages(dataset2_plgls_model_variances_homo, 
                                             dataset2_true_params_struct)

dataset3_plgls_coverage_hetero = get_coverages(dataset3_plgls_model_variances_hetero, 
                                               dataset3_true_params_struct)
dataset3_plgls_coverage_homo = get_coverages(dataset3_plgls_model_variances_homo, 
                                             dataset3_true_params_struct)

dataset1_plgls_coverage_hetero[, VarType := "hetero"]
dataset1_plgls_coverage_homo[, VarType := "homo"]
dataset2_plgls_coverage_hetero[, VarType := "hetero"]
dataset2_plgls_coverage_homo[, VarType := "homo"]
dataset3_plgls_coverage_hetero[, VarType := "hetero"]
dataset3_plgls_coverage_homo[, VarType := "homo"]

# Final tables ----

print(xtable::xtable(dcast(dataset1_plgls_coverage_hetero,
                     segment_id + param_id ~ error_sd, value.var = "coverage"), digits = 0), include.rownames = F)
print(xtable::xtable(dcast(dataset1_plgls_coverage_homo,
                           segment_id + param_id ~ error_sd, value.var = "coverage"), digits = 0), include.rownames = F)

print(xtable::xtable(dcast(rbind(dataset1_plgls_coverage_hetero, dataset1_plgls_coverage_homo),
                           segment_id + param_id ~ error_sd + VarType, value.var = "coverage"), digits = 0), include.rownames = F)


print(xtable::xtable(dcast(dataset2_plgls_coverage_hetero,
                           segment_id + param_id ~ error_sd, value.var = "coverage"), digits = 0), include.rownames = F)
print(xtable::xtable(dcast(dataset2_plgls_coverage_homo,
                           segment_id + param_id ~ error_sd, value.var = "coverage"), digits = 0), include.rownames = F)

print(xtable::xtable(dcast(rbind(dataset2_plgls_coverage_hetero, dataset2_plgls_coverage_homo),
                           segment_id + param_id ~ error_sd + VarType, value.var = "coverage"), digits = 0), include.rownames = F)


print(xtable::xtable(dcast(dataset3_plgls_coverage_hetero,
                           segment_id + param_id ~ error_sd, value.var = "coverage"), digits = 0), include.rownames = F)
print(xtable::xtable(dcast(dataset3_plgls_coverage_homo,
                           segment_id + param_id ~ error_sd, value.var = "coverage"), digits = 0), include.rownames = F)

print(xtable::xtable(dcast(rbind(dataset3_plgls_coverage_hetero, dataset3_plgls_coverage_homo),
                           segment_id + param_id ~ error_sd + VarType, value.var = "coverage"), digits = 0), include.rownames = F)


dataset1_plgls_coverage_hetero[, .(AveSE = mean(se)),
                               by = c("error_sd", "VarType")]
dataset1_plgls_coverage_homo[, .(AveSE = mean(se)),
                               by = c("error_sd", "VarType")]

dataset2_plgls_coverage_hetero[, .(AveSE = mean(se)),
                               by = c("error_sd", "VarType")]
dataset2_plgls_coverage_homo[, .(AveSE = mean(se)),
                             by = c("error_sd", "VarType")]


dataset3_plgls_coverage_hetero[, .(AveSE = mean(se)),
                               by = c("error_sd", "VarType")]
dataset3_plgls_coverage_homo[, .(AveSE = mean(se)),
                             by = c("error_sd", "VarType")]

# TODO: OLS vs PL-GLS variance plots for all data sets based on the following
# comp_all_homo_plgls = model_variances_homo_plgls[, .(mean_model_var = mean(model_variance),
#                                                      empirical_var = var(estimated),
#                                                      true_value = unique(true_params),
#                                                      mean_estimated = mean(estimated),
#                                                      n = .N),
#                                                  by = c("param_id", "error_sd")]
# ggplot(rbind(
#   cbind(comp_all_homo_plgls, method = "PL-GLS"),
#   cbind(comp_all_homo_ols, method = "OLS"))[mean_model_var > 0], aes(x = empirical_var, y = mean_model_var, color = method)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   facet_wrap(~error_sd, scales = "free") +
#   xlab("empirical variance") +
#   ylab("mean model-based variance") +
#   theme_bw() +
#   scale_color_discrete(palette = "viridis") +
#   # xlim(c(0, 0.1)) +
#   # ylim(c(0, 0.1)) +
#   theme(legend.position = "bottom",
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 14),
#         legend.text = element_text(size = 14),
#         legend.title = element_text(size = 14),
#         strip.text = element_text(size = 14))
# ggsave("medcl_homo_var_comp.png", device = "png", scale = 1,
#        width = 7, height = 5, units = "in", dpi = 300)

