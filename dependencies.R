
# Installation instructions for roplspvs (Ropls-with-permutated-variable-selection)
#This package was created for r 3.6.2 using rstudio 1.1.463 but has been tested
#for R 3.6.2, 4.0.0, 4.0.3, 4.0.5 and 4.1.2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")
install.packages("tools")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("kableExtra")
install.packages("gridExtra")
install.packages("ggpubr")
install.packages("matrixStats")
install.packages("stringr")
install.packages("tryCatchLog")
install.packages("devtools")
install.packages("DescTools")
install.packages("precrec")
install.packages("pROC")
install.packages("rstatix")
install.packages("rmarkdown")

