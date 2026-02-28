# Libraries ---- 
library(data.table)

# Functions ---- 
get_plgls_info = function(results_list, true_params) {
  model_variances_homo_plgls = rbindlist(lapply(results_list, function(x) {
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

get_wide_form_data = function(comparison_table, true_info) {
  summary_table = comparison_table[, .(mean_model_var = mean(model_variance),
                                                 empirical_var = var(estimated),
                                                 true_value = unique(true_params),
                                                 mean_estimated = mean(estimated),
                                                 n = .N),
                                             by = c("param_id", "error_sd")]
  wide_form = dcast(summary_table, param_id + true_value ~ error_sd, 
                    value.var = c("n", "mean_estimated", "mean_model_var", "empirical_var"))
  wide_form = merge(wide_form, true_info, by = "param_id")
  wide_form
}

get_summary_table = function(raw_results, true_params_struct, true_params) {
  comparison_table = get_plgls_info(raw_results, true_params)
  true_info = get_true_info(true_params_struct)
  wide_form_table = get_wide_form_data(comparison_table, true_info)
  wide_form_table
}

# Data and results for each case ---- 

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

# Raw tables ----

dataset1_hetero_tbl = get_summary_table(dataset1_hetero_res,
                                        dataset1_true_params_struct,
                                        dataset1_true_params)
dataset1_homo_tbl = get_summary_table(dataset1_homo_res,
                                      dataset1_true_params_struct,
                                      dataset1_true_params)

dataset2_hetero_tbl = get_summary_table(dataset2_hetero_res,
                                        dataset2_true_params_struct,
                                        dataset2_true_params)
dataset2_homo_tbl = get_summary_table(dataset2_homo_res,
                                      dataset2_true_params_struct,
                                      dataset2_true_params)

dataset3_hetero_tbl = get_summary_table(dataset3_hetero_res,
                                        dataset3_true_params_struct,
                                        dataset3_true_params)
dataset3_homo_tbl = get_summary_table(dataset3_homo_res,
                                      dataset3_true_params_struct,
                                      dataset3_true_params)

# Final tables to export ---- 


# Data set 1 
print(xtable::xtable(
  dataset1_hetero_tbl[, .(segment_id, prob, true_value.x,
                          mean_estimated_0.01, 
                          rel_error_0.01 = 100 * (mean_estimated_0.01 - true_value.x) / true_value.x, 
                          1e4 * empirical_var_0.01, 1e4 * mean_model_var_0.01,
                          mean_estimated_0.05, 
                          rel_error_0.05 = 100 * (mean_estimated_0.05 - true_value.x) / true_value.x,
                          1e4 * empirical_var_0.05, 1e4 * mean_model_var_0.05)],
  digits = 3), include.rownames = F)
print(xtable::xtable(
  dataset1_homo_tbl[, .(segment_id, prob, true_value.x,
                        mean_estimated_0.5, 
                        rel_error_0.5 = 100 * (mean_estimated_0.5 - true_value.x) / true_value.x, 
                        empirical_var_0.5, mean_model_var_0.5,
                        mean_estimated_1, 
                        rel_error_1 = 100 * (mean_estimated_1 - true_value.x) / true_value.x,
                        empirical_var_1, mean_model_var_1)],
  digits = 3), include.rownames = F)


print(xtable::xtable(
  dataset2_hetero_tbl[, .(segment_id, prob, true_value.x,
                          mean_estimated_0.01, 
                          rel_error_0.01 = 100 * (mean_estimated_0.01 - true_value.x) / true_value.x, 
                          1e4 * empirical_var_0.01, 1e4 * mean_model_var_0.01,
                          mean_estimated_0.05, 
                          rel_error_0.05 = 100 * (mean_estimated_0.05 - true_value.x) / true_value.x,
                          1e4 * empirical_var_0.05, 1e4 * mean_model_var_0.05)],
  digits = 3), include.rownames = F)
print(xtable::xtable(
  dataset2_homo_tbl[, .(segment_id, prob, true_value.x,
                        mean_estimated_0.5, 
                        rel_error_0.5 = 100 * (mean_estimated_0.5 - true_value.x) / true_value.x, 
                        empirical_var_0.5, mean_model_var_0.5,
                        mean_estimated_1, 
                        rel_error_1 = 100 * (mean_estimated_1 - true_value.x) / true_value.x,
                        empirical_var_1, mean_model_var_1)],
  digits = 3), include.rownames = F)

print(xtable::xtable(
  dataset3_hetero_tbl[, .(segment_id, prob, true_value.x,
                          mean_estimated_0.01, 
                          rel_error_0.01 = 100 * (mean_estimated_0.01 - true_value.x) / true_value.x, 
                          1e2 * empirical_var_0.01, 1e2 * mean_model_var_0.01,
                          mean_estimated_0.05, 
                          rel_error_0.05 = 100 * (mean_estimated_0.05 - true_value.x) / true_value.x,
                          1e2 * empirical_var_0.05, 1e2 * mean_model_var_0.05)],
  digits = 3), include.rownames = F)
print(xtable::xtable(
  dataset3_homo_tbl[, .(segment_id, prob, true_value.x,
                        mean_estimated_0.5, 
                        rel_error_0.5 = 100 * (mean_estimated_0.5 - true_value.x) / true_value.x, 
                        empirical_var_0.5, mean_model_var_0.5,
                        mean_estimated_1, 
                        rel_error_1 = 100 * (mean_estimated_1 - true_value.x) / true_value.x,
                        empirical_var_1, mean_model_var_1)],
  digits = 3), include.rownames = F)
