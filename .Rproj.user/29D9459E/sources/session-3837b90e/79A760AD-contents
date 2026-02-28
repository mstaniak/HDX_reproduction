# Libraries ----
library(data.table)

# Functions ----

get_theta_estimates = function(results_list) {
  rbindlist(lapply(results_list, function(x) {
    if (!is.null(x)) {
      if (!is.null(x$fitted_plgls)) {
        if (max(abs(x$fitted_plgls$fvec)) < 1e-6) {
          list(theta = x$theta,
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

summarize_theta_estimates_per_dataset = function(theta_estimates_hetero,
                                                 theta_estimates_homo) {
  rbind(cbind(theta_estimates_hetero, VarType = "hetero"),
        cbind(theta_estimates_homo, VarType = "homo"))[, .(avetheta = mean(theta, na.rm = T),
                                                           vartheta = var(theta, na.rm = T)),
                                                       by = c("VarType", "error_sd")]
}

# Results of simulations ----
dataset1_hetero_res = readRDS("./RawResults/dataset1_heteroscedastic.RDS")
dataset1_homo_res = readRDS("./RawResults/dataset1_homoscedastic.RDS")

dataset2_hetero_res = readRDS("./RawResults/dataset2_heteroscedastic.RDS")
dataset2_homo_res = readRDS("./RawResults/dataset2_homoscedastic.RDS")

dataset3_hetero_res = readRDS("./RawResults/dataset3_heteroscedastic.RDS")
dataset3_homo_res = readRDS("./RawResults/dataset3_homoscedastic.RDS")

# Estimates of theta ----

dataset1_thetas_hetero = get_theta_estimates(dataset1_hetero_res)
dataset1_thetas_homo = get_theta_estimates(dataset1_homo_res)

dataset2_thetas_hetero = get_theta_estimates(dataset2_hetero_res)
dataset2_thetas_homo = get_theta_estimates(dataset2_homo_res)

dataset3_thetas_hetero = get_theta_estimates(dataset3_hetero_res)
dataset3_thetas_homo = get_theta_estimates(dataset3_homo_res)

# Summaries ----
dataset1_theta_summaries = summarize_theta_estimates_per_dataset(dataset1_thetas_hetero,
                                                                 dataset1_thetas_homo)
dataset2_theta_summaries = summarize_theta_estimates_per_dataset(dataset2_thetas_hetero,
                                                                 dataset2_thetas_homo)
dataset3_theta_summaries = summarize_theta_estimates_per_dataset(dataset3_thetas_hetero,
                                                                 dataset3_thetas_homo)
# Final table ----

theta_summary = cbind(cluster = rep(c("Data set 1",
                      "Data set 2",
                      "Data set 3"),
                    each = 4),
      rbind(dataset1_theta_summaries,
            dataset2_theta_summaries,
            dataset3_theta_summaries))

print(xtable::xtable(dcast(theta_summary,
      cluster ~ VarType + error_sd,
      value.var = c("avetheta", "vartheta"))[, .(cluster, avetheta_homo_0.5, vartheta_homo_0.5,
                                                 avetheta_homo_1, vartheta_homo_1,
                                                 avetheta_hetero_0.01, vartheta_hetero_0.01,
                                                 avetheta_hetero_0.05, vartheta_hetero_0.05)],
      digits = 4), include.rownames = F)
