# Code and data required to reproduce "A Statistical Model to Resolve High-Resolution HDX Dynamics through Isotopic Distributions of Overlapping Peptides" analyses

This repository contains data and R scripts required to reproduces results presented in the article "A Statistical Model to Resolve High-Resolution HDX Dynamics through Isotopic
Distributions of Overlapping Peptides" (Staniak, Claesen, Burzykowski, 2026). 

1. Simulated and real data stored in the Data/ folder were pre-processed to fit into a format required by the IsoHDX package. Raw data are available on request. 
Similarly, the repository stores fitted models. Code used to generate those fits is available on request. 

2. The repository consists of the following files:
  
  * 00_setup.R - installation of R packages necessary to run the analyses,
  * 01_simulated_data_viz.R - plots of true spectra and exchange probablities for simulated data,
  * 02_summary_tables_plgls.R - tables and plots describing bias and variance of PL-GLS-based estimates,
  * 03_ols_gls_pointwise_comp.R - comparisons of point-wise estimates based on OLS and PL-GLS approaches,
  * 04_plgls_variance_and_coverage.R - coverages of CIs for model parameters based on PL-GLS,
  * 05_exchange_probabilities.R - plots of point estimates of exchange probabilities,
  * 06_theta_estimation_plgls.R - table describing estimation of the $\theta$ parameter in the PL-GLS approach,
  * 07_ols_plgls_variances_comp.R - comparisons of model-variances to empirical variances for OLS and PL-GLS approaches,
  * 08_tech_reps.R - plots and tables describing estimation with additional technical replicates,
  * 09_prob_confints.R - coverages of CIs for exchange probabilities,
  * 10_perturbed_spectra.R - results of estimation with spurious peaks,
  * 11_hvem_plots.R - visualization of the model fitted for the HVEM case study,
  * 12_glypb_plots.R - visualization of the model fitted for the GlyPb case study.
