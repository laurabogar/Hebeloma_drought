# Preliminary analysis of 3' Tag Seq data
# 12 Sept 2022

# Step 2: Attempt differential expression analysis

library(DESeq2)
library(tidyverse)

samples = read_excel("sample_to_condition_table.xlsx")
tpm = read_csv("kallisto_tpm_output.csv")
est_counts = read_csv("kallisto_est_counts.csv")

