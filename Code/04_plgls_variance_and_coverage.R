# Libraries ---- 
library(data.table)
library(ggplot2)

# Functions ----
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
