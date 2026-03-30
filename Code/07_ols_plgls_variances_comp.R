# Libraries ----
library(data.table)
library(ggplot2)

# Functions ----
get_ols_model_variances = function(numerical_hessians, estimation_results,
                                   ps_m, num_parameters, true_params, 
                                   undeuterated_dists) {
  inverse_hessians = lapply(numerical_hessians, function(x) {
    if (!is.null(x)) {
      if (!is.null(x$num_hessian)) {
        x$inv_hess = solve(x$num_hessian) 
      }
    }
    x
  })
  
  sigma_sq_estim = lapply(estimation_results, function(x) {
    if (!is.null(x)) {
      if (!is.null(x$fitted_ols)) {
        if (max(abs(x$fitted_ols$fvec)) < 1e-6) {
          predicted = IsoHDX3:::getExpectedSpectra(x$fitted_ols$x, ps_m, num_parameters + 1, x$data, undeuterated_dists)
          comp = merge(x$data, predicted, by = c("Rep", "Peptide", "Time", "IntDiff"))
          x$comp = comp
          x$sigma_sq = (1 / (nrow(x$data) - length(true_params))) * comp[, sum((Intensity - ExpectedPeak)^2)]
          return(x)
        } else {
          NULL
        }
      } else {
        NULL
      }   
    } else {
      NULL
    }
  })
  
  varcov_estimates_data = lapply(seq_along(inverse_hessians), function(i) {
    x = sigma_sq_estim[[i]]
    x$inv_hess = inverse_hessians[[i]]$inv_hess
    x
  })
  
  varcov_estimates = lapply(varcov_estimates_data, function(x) {
    x$varcov = x$sigma_sq * x$inv_hess
    x
  })
  
  
  model_variances_tbl = rbindlist(lapply(varcov_estimates, function(x) {
    if (!is.null(x)) {
      if (!is.null(x$fitted_ols)) {
        if (max(abs(x$fitted_ols$fvec)) < 1e-6) {
          list(param_id = seq_along(true_params),
               true_params = true_params,
               estimated = x$fitted_ols$x,
               model_variance = diag(x$varcov),
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
  model_variances_tbl
}

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

get_variances_comparison = function(model_variances_tbl, var_type, method) {
  res = model_variances_tbl[, .(mean_model_var = mean(model_variance),
                                empirical_var = var(estimated),
                                true_value = unique(true_params),
                                mean_estimated = mean(estimated),
                                n = .N),
                            by = c("param_id", "error_sd")]
  res[, VarType := var_type]
  res[, method := method]
  res
}

plot_variances_comparison = function(model_variances_ols,
                                     model_variances_plgls) {
  ggplot(rbind(model_variances_ols,
               model_variances_plgls)[mean_model_var > 0], aes(x = empirical_var, y = mean_model_var, color = method)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    facet_wrap(~error_sd, scales = "free") +
    xlab("empirical variance") +
    ylab("mean model-based variance") +
    theme_bw() +
    scale_color_discrete(palette = "viridis") +
    # xlim(c(0, 0.1)) +
    # ylim(c(0, 0.1)) +
    theme(legend.position = "bottom",
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          strip.text = element_text(size = 14))
}


# Simulated data ----

dataset1_ps_m = readRDS("./Data/Simulated/DataSet1/ps_m.RDS")
dataset1_und_d = readRDS("./Data/Simulated/DataSet1/undeuterated_dists.RDS")
dataset1_num_pars = sapply(readRDS("./Data/Simulated/DataSet1/true_params_struct.RDS"),
                           length) - 1
dataset1_true_params = readRDS("./Data/Simulated/DataSet1/true_params.RDS")

dataset2_ps_m = readRDS("./Data/Simulated/DataSet2/ps_m.RDS")
dataset2_und_d = readRDS("./Data/Simulated/DataSet2/undeuterated_dists.RDS")
dataset2_num_pars = sapply(readRDS("./Data/Simulated/DataSet2/true_params_struct.RDS"),
                           length) - 1
dataset2_true_params = readRDS("./Data/Simulated/DataSet2/true_params.RDS")

dataset3_ps_m = readRDS("./Data/Simulated/DataSet3/ps_m.RDS")
dataset3_und_d = readRDS("./Data/Simulated/DataSet3/undeuterated_dists.RDS")
dataset3_num_pars = sapply(readRDS("./Data/Simulated/DataSet3/true_params_struct.RDS"),
                            length) - 1
dataset3_true_params = readRDS("./Data/Simulated/DataSet3/true_params.RDS")

# Results of simulation and numerical Hessians 

dataset1_results_homo = readRDS("./RawResults/dataset1_homoscedastic.RDS")
dataset1_results_hetero = readRDS("./RawResults/dataset1_heteroscedastic.RDS")

dataset1_num_hessians_homo = readRDS("./RawResults/dataset1_num_hessians_homo.RDS")
dataset1_num_hessians_hetero = readRDS("./RawResults/dataset1_num_hessians_hetero.RDS")

dataset2_results_homo = readRDS("./RawResults/dataset2_homoscedastic.RDS")
dataset2_results_hetero = readRDS("./RawResults/dataset2_heteroscedastic.RDS")

dataset2_num_hessians_homo = readRDS("./RawResults/dataset2_num_hessians_homo.RDS")
dataset2_num_hessians_hetero = readRDS("./RawResults/dataset2_num_hessians_hetero.RDS")

dataset3_results_homo = readRDS("./RawResults/dataset3_homoscedastic.RDS")
dataset3_results_hetero = readRDS("./RawResults/dataset3_heteroscedastic.RDS")

dataset3_num_hessians_homo = readRDS("./RawResults/dataset3_num_hessians_homo.RDS")
dataset3_num_hessians_hetero = readRDS("./RawResults/dataset3_num_hessians_hetero.RDS")

# Results

dataset1_ols_model_variances_homo = get_ols_model_variances(dataset1_num_hessians_homo,
                                                            dataset1_results_homo,
                                                            dataset1_ps_m,
                                                            dataset1_num_pars,
                                                            dataset1_true_params,
                                                            dataset1_und_d)
dataset1_plgls_model_variances_homo = get_plgls_model_variances(dataset1_results_homo,
                                                                dataset1_true_params)

dataset1_ols_model_variances_comp_homo = get_variances_comparison(dataset1_ols_model_variances_homo,
                                                                  "homoscedastic", "OLS")
dataset1_plgls_model_variances_comp_homo = get_variances_comparison(dataset1_plgls_model_variances_homo,
                                                                    "homoscedastic", "PLG-GLS")

dataset1_ols_model_variances_hetero = get_ols_model_variances(dataset1_num_hessians_hetero,
                                                              dataset1_results_hetero,
                                                              dataset1_ps_m,
                                                              dataset1_num_pars,
                                                              dataset1_true_params,
                                                              dataset1_und_d)
dataset1_plgls_model_variances_hetero = get_plgls_model_variances(dataset1_results_hetero,
                                                                  dataset1_true_params)

dataset1_ols_model_variances_comp_hetero = get_variances_comparison(dataset1_ols_model_variances_hetero,
                                                                    "heteroscedastic", "OLS")
dataset1_plgls_model_variances_comp_hetero = get_variances_comparison(dataset1_plgls_model_variances_hetero,
                                                                      "heteroscedastic", "PLG-GLS")

dataset2_ols_model_variances_homo = get_ols_model_variances(dataset2_num_hessians_homo,
                                                            dataset2_results_homo,
                                                            dataset2_ps_m,
                                                            dataset2_num_pars,
                                                            dataset2_true_params,
                                                            dataset2_und_d)
dataset2_plgls_model_variances_homo = get_plgls_model_variances(dataset2_results_homo,
                                                                dataset2_true_params)

dataset2_ols_model_variances_comp_homo = get_variances_comparison(dataset2_ols_model_variances_homo,
                                                                  "homoscedastic", "OLS")
dataset2_plgls_model_variances_comp_homo = get_variances_comparison(dataset2_plgls_model_variances_homo,
                                                                    "homoscedastic", "PLG-GLS")

dataset2_ols_model_variances_hetero = get_ols_model_variances(dataset2_num_hessians_hetero,
                                                              dataset2_results_hetero,
                                                              dataset2_ps_m,
                                                              dataset2_num_pars,
                                                              dataset2_true_params,
                                                              dataset2_und_d)
dataset2_plgls_model_variances_hetero = get_plgls_model_variances(dataset2_results_hetero,
                                                                  dataset2_true_params)

dataset2_ols_model_variances_comp_hetero = get_variances_comparison(dataset2_ols_model_variances_hetero,
                                                                    "heteroscedastic", "OLS")
dataset2_plgls_model_variances_comp_hetero = get_variances_comparison(dataset2_plgls_model_variances_hetero,
                                                                      "heteroscedastic", "PLG-GLS")


dataset3_ols_model_variances_homo = get_ols_model_variances(dataset3_num_hessians_homo,
                                                       dataset3_results_homo,
                                                       dataset3_ps_m,
                                                       dataset3_num_pars,
                                                       dataset3_true_params,
                                                       dataset3_und_d)
dataset3_plgls_model_variances_homo = get_plgls_model_variances(dataset3_results_homo,
                                                                dataset3_true_params)

dataset3_ols_model_variances_comp_homo = get_variances_comparison(dataset3_ols_model_variances_homo,
                                                                  "homoscedastic", "OLS")
dataset3_plgls_model_variances_comp_homo = get_variances_comparison(dataset3_plgls_model_variances_homo,
                                                                    "homoscedastic", "PLG-GLS")

dataset3_ols_model_variances_hetero = get_ols_model_variances(dataset3_num_hessians_hetero,
                                                            dataset3_results_hetero,
                                                            dataset3_ps_m,
                                                            dataset3_num_pars,
                                                            dataset3_true_params,
                                                            dataset3_und_d)
dataset3_plgls_model_variances_hetero = get_plgls_model_variances(dataset3_results_hetero,
                                                                dataset3_true_params)

dataset3_ols_model_variances_comp_hetero = get_variances_comparison(dataset3_ols_model_variances_hetero,
                                                                  "heteroscedastic", "OLS")
dataset3_plgls_model_variances_comp_hetero = get_variances_comparison(dataset3_plgls_model_variances_hetero,
                                                                    "heteroscedastic", "PLG-GLS")


plot_variances_comparison(dataset1_ols_model_variances_comp_homo,
                          dataset1_plgls_model_variances_comp_homo)
ggsave("./Figures/smallcl_ols_plgls_var_comp_homo.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)
plot_variances_comparison(dataset1_ols_model_variances_comp_hetero,
                          dataset1_plgls_model_variances_comp_hetero)
ggsave("./Figures/smallcl_ols_plgls_var_comp_hetero.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)

plot_variances_comparison(dataset2_ols_model_variances_comp_homo,
                          dataset2_plgls_model_variances_comp_homo)
ggsave("./Figures/medcl_ols_plgls_var_comp_homo.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)
plot_variances_comparison(dataset2_ols_model_variances_comp_hetero,
                          dataset2_plgls_model_variances_comp_hetero)
ggsave("./Figures/medcl_ols_plgls_var_comp_hetero.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)


plot_variances_comparison(dataset3_ols_model_variances_comp_homo,
                          dataset3_plgls_model_variances_comp_homo)
ggsave("./Figures/shortsegm_ols_plgls_var_comp_homo.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)
plot_variances_comparison(dataset3_ols_model_variances_comp_hetero,
                          dataset3_plgls_model_variances_comp_hetero)
ggsave("./Figures/shortsegm_ols_plgls_var_comp_hetero.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)

