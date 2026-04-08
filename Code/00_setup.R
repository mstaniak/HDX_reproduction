if (!require("devtools")) {
  install.packages("devtools")
}
if (!require("BiocManager")) {
  install.packages("BiocManager")
}
install.packages(c("data.table", "ggplot2", "viridis"))
BiocManager::install("BRAIN")
devtools::install("mstaniak/IsoHDX", ref = "polynom-implementation")