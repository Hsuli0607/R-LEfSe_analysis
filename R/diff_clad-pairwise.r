### Load Required Libraries ###
# This script assumes all required packages have been installed.
# Use install.packages("package_name") or BiocManager::install("package_name") if any are missing.

library(phyloseq)
library(microbiomeMarker) # For LEfSe analysis
library(dplyr)
library(data.table)       # For the fread function
library(tidyr)
library(ape)




lefse_result <- lefse_analysis(physeq_object = physeq)





### Plot Diff_Cladogram Function
