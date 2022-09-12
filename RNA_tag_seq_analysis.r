# Preliminary analysis of 3' Tag Seq data
# 12 Sept 2022

# Step 1: Get kallisto tpm estimates into a single table

library(tidyverse)
library(DESeq2)

kallisto_output = "kallisto_output/"

list.files(kallisto_output)

samples = read_excel("sample_to_condition_table.xlsx")

head(samples)
