library(data.table)
library(ggplot2)


get_ols_gls_comparison = function(results_list, true_params) {
  pointwise_ols = rbindlist(lapply(results_list, function(x) {
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
  
  summary_ols = pointwise_ols[, .(true_value = unique(true_params),
                                  mean_estimated = mean(estimated),
                                  n = .N),
                              by = c("param_id", "error_sd")]
  
  
  pointwise_gls = rbindlist(lapply(results_list, function(x) {
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
  summary_gls = pointwise_gls[, .(true_value = unique(true_params),
                                  mean_estimated = mean(estimated),
                                  n = .N),
                              by = c("param_id", "error_sd")]
  
  ols_gls_comparison = rbind(
    cbind(summary_gls, method = "PL-GLS"),
    cbind(summary_ols, method = "OLS"))
  ols_gls_comparison
}

plot_pointwise_comparison = function(pointwise_comparison) {
  ggplot(pointwise_comparison,
         aes(x = true_value, y = mean_estimated, color = method, shape = method)) +
    geom_point(size = 3)  +
    # geom_label(aes(label = param_id)) +
    # scale_color_discrete(palette = "viridis") +
    scale_color_manual(values = rev(viridis::viridis(3)[-1])) +
    geom_abline(slope = 1, intercept = 0) +
    xlab("true value of the parameter") +
    ylab("average estimate") +
    # geom_point(aes(y = true_value), color = "red", size = 1.5) +
    facet_grid(~error_sd) +
    theme_bw() +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)) 
}


dataset1_true_params = readRDS("./Data/Simulated/DataSet1/true_params.RDS")
dataset1_hetero_res = readRDS("./RawResults/dataset1_heteroscedastic.RDS")
dataset1_homo_res = readRDS("./RawResults/dataset1_homoscedastic.RDS")

dataset2_true_params = readRDS("./Data/Simulated/DataSet2/true_params.RDS")
dataset2_hetero_res = readRDS("./RawResults/dataset2_heteroscedastic.RDS")
dataset2_homo_res = readRDS("./RawResults/dataset2_homoscedastic.RDS")

dataset3_true_params = readRDS("./Data/Simulated/DataSet3/true_params.RDS")
dataset3_hetero_res = readRDS("./RawResults/dataset3_heteroscedastic.RDS")
dataset3_homo_res = readRDS("./RawResults/dataset3_homoscedastic.RDS")

dataset1_hetero_comparison = get_ols_gls_comparison(dataset1_hetero_res, dataset1_true_params)
dataset1_homo_comparison = get_ols_gls_comparison(dataset1_homo_res, dataset1_true_params)

plot_pointwise_comparison(dataset1_hetero_comparison)
ggsave("./Figures/smallcl_hetero_pointwise_comp.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)
plot_pointwise_comparison(dataset1_homo_comparison)
ggsave("./Figures/smallcl_homo_pointwise_comp.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)

dataset2_hetero_comparison = get_ols_gls_comparison(dataset2_hetero_res, dataset2_true_params)
dataset2_homo_comparison = get_ols_gls_comparison(dataset2_homo_res, dataset2_true_params)

plot_pointwise_comparison(dataset2_hetero_comparison)
ggsave("./Figures/medcl_hetero_pointwise_comp.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)
plot_pointwise_comparison(dataset2_homo_comparison)
ggsave("./Figures/medcl_homo_pointwise_comp.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)

dataset3_hetero_comparison = get_ols_gls_comparison(dataset3_hetero_res, dataset3_true_params)
dataset3_homo_comparison = get_ols_gls_comparison(dataset3_homo_res, dataset3_true_params)

plot_pointwise_comparison(dataset3_hetero_comparison)
ggsave("./Figures/shortsegm_hetero_pointwise_comp.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)
plot_pointwise_comparison(dataset3_homo_comparison)
ggsave("./Figures/shortsegm_homo_pointwise_comp.png", device = "png", scale = 1,
       width = 7, height = 5, units = "in", dpi = 300)

# ggsave("medcl_homo_estim_comp.png", device = "png", scale = 1,
#        width = 7, height = 5, units = "in", dpi = 300)

