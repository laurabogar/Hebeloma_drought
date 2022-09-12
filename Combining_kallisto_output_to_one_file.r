# Preliminary analysis of 3' Tag Seq data
# 12 Sept 2022

# Step 1: Get kallisto tpm estimates into a single table

library(tidyverse)
library(readxl)

kallisto_output = "kallisto_output/"

samples = list.files(kallisto_output)

firstoutput = read_tsv("kallisto_output/HC007/abundance.tsv")

# samplename <- function(df, index) {
#   mutate(df, "" := mean({{col}}))
# }

tpm = select(firstoutput, target_id)
for (i in 1:length(samples)) {
  newsamplepath = (paste("kallisto_output", samples[i], "abundance.tsv", sep = "/"))
  newsample_ID = samples[i]
  abund = as.data.frame(read_table(newsamplepath))
  tojoin = as.data.frame(cbind(abund$target_id, abund$tpm)) # can't use tibble for dynamic column naming
  colnames(tojoin)[1] = "target_id" # base R naming so each column has own sample ID
  colnames(tojoin)[2] = newsample_ID # here's the dynamic naming
  tpm = left_join(tpm, tojoin) # pasting it back together into a tibble for output
}

est_counts = select(firstoutput, target_id)
for (i in 1:length(samples)) {
  newsamplepath = (paste("kallisto_output", samples[i], "abundance.tsv", sep = "/"))
  newsample_ID = samples[i]
  abund = as.data.frame(read_table(newsamplepath))
  tojoin = as.data.frame(cbind(abund$target_id, abund$est_counts)) # can't use tibble for dynamic column naming
  colnames(tojoin)[1] = "target_id" # base R naming so each column has own sample ID
  colnames(tojoin)[2] = newsample_ID # here's the dynamic naming
  est_counts = left_join(est_counts, tojoin) # pasting it back together into a tibble for output
}

write_csv(tpm, "kallisto_tpm_output.csv")
write_csv(est_counts, "kallisto_est_counts.csv")

