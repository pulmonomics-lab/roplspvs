#This package was created for r 3.6.2 using rstudio 1.1.463
#Run this code line by line
install.packages("ropls")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("kableExtra")
install.packages("gridExtra")
install.packages("ggpubr")
install.packages("matrixStats")
install.packages("stringr")
install.packages("tryCatchLog")

#if install.packages("ropls") does not work run:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ropls")

