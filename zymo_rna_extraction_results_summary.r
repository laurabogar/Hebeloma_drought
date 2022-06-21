# Zymo RNA extraction results

library("tidyverse")
library("cowplot")

results = read_csv("Hebeloma RNA quantification and plan.csv")

write_csv(results, "Hebeloma_RNA_quant_and_plan_original.csv")

ggplot(subset(results, `% loss to clean and concentrate kit` >=0), aes(x = `% loss to clean and concentrate kit`)) +
  geom_histogram() +
  theme_cowplot() +

