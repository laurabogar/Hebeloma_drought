# Making some preliminary figures for the Hebeloma data

library(tidyverse)
library(cowplot)
library(lme4)
library(lmerTest)

fulldata = read_csv("Merged Hebeloma 11Apr2022 OT.csv")

fulldata$percent_col = 100*(fulldata$colonized_tips/fulldata$total_tips)

fulldata$water_level = as.factor(fulldata$water_level)

# Put water levels in correct order for plotting
fulldata$water_level = factor(fulldata$water_level, levels = c("L", "M", "H"))

weird_point = fulldata[fulldata$water_level == "L" & 
                         fulldata$n == "Minus" &
                         fulldata$percent_col <100,]
# I had no plants in this project that were Minus N and also low water.
# what's up with this point?
# Checked original data sheet photo (from 2/14) and found that it was mis-transcribed.
# this is plant HC004

fulldata$n[fulldata$seedling == "HC004"] = "Plus"

fulldata$minutes_in_cooler = difftime(fulldata$approximate_measurement_time,
                                   fulldata$time,
                                   units = "mins")

fulldata$harvest_day = as.factor(fulldata$harvest_day)
fulldata$drydown_day = recode_factor(fulldata$harvest_day,
                              "1" = "1",
                              "2" = "3",
                              "3" = "7",
                              "4" = "10")

# fulldata$drydown_day = numeric(nrow(fulldata))
# for (i in 1:nrow(fulldata)) {
#   if (fulldata$harvest_day[i] == "1") {
#     fulldata$drydown_day[i] = 1
#   } else if (fulldata$harvest_day[i] == "2") {
#     fulldata$drydown_day[i] = 3
#   } else if (fulldata$harvest_day[i] == "3") {
#     fulldata$drydown_day[i] = 7
#   } else if (fulldata$harvest_day[i] == "4") {
#     fulldata$drydown_day[i] = 10
#   }
# }


colonized = subset(fulldata, colonized == "Y")

# Percent colonization by water treatment
ggplot(data = colonized) +
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(0.9),
               aes(x = n, y = percent_col,
                   color = water_level)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.2),
              aes(x = n, y = percent_col,
                  color = water_level)) +
  xlab("Nitrogen treatment") +
  ylab("Percent ectomycorrhizal colonization")

# water potentials
ggplot(data = fulldata) +
  theme_cowplot() +
  geom_jitter(width = 0.25,
              aes(x = drydown_day, y = water_potential_m_pa,
                 color = water_level,
                 shape = colonized)) +
  xlab("Days of drydown") +
  ylab("Water potential")
# Hmm not a lot of variation over time.

# Checking to see if plant water potentials got more negative if the plants sat
# for a long time in the cooler before measurement.

# Let's just check one day at a time, since I think I did days 1 and 2 as
# one big batch and after that I did two measurement batches per day.
# I started measuring on Day 1 around 4PM; similar for day 2.

ggplot(data = subset(fulldata, harvest_day == 1)) +
  theme_cowplot() +
  geom_point(aes(x = time, y = water_potential_m_pa,
                  color = water_level,
                  shape = colonized)) +
  xlab("Time at harvest") +
  ylab("Water potential (MPa)")

# This looks okay to me, although it's possible things got drier as they sat
# (points on the left side are a bit drier)

# What about day 2?

ggplot(data = subset(fulldata, harvest_day == 2)) +
  theme_cowplot() +
  geom_point(aes(x = time, y = water_potential_m_pa,
                 color = water_level,
                 shape = colonized)) +
  xlab("Time at harvest") +
  ylab("Water potential (MPa)")

# Let's look at all days at once:

ggplot(data = fulldata) +
  theme_cowplot() +
  geom_point(aes(x = minutes_in_cooler, 
                 y = water_potential_m_pa,
                 color = water_level,
                 shape = colonized)) +
  xlab("Minutes since harvest") +
  ylab("Water potential (MPa)")

# This is encouraging. I don't see a pattern of drier plants after more time in the cooler.

# Formalizing with stats:

cooler_check = lmer(water_potential_m_pa ~ minutes_in_cooler + (1|colonized) + (1|n) + (1|harvest_day) + (1|water_level), 
                    data = fulldata)
summary(cooler_check)

cooler_check_simpler = lm(water_potential_m_pa ~ minutes_in_cooler, 
                    data = fulldata)
summary(cooler_check_simpler)
# Okay, this is kind of a bummer, but both models agree:
# The longer a shoot spent in the cooler, the drier it was.
# However, the effect size estimate is MINISCULE:
# Plants became 0.001 MPa drier with each minute in the cooler.
# (LME says -0.0015; normal lm says -0.0012.)
# (Similar effect predicted by both models.)
# So, in an hour, we'd see a difference of 0.06 MPa.
# Our absolute longest-sitting samples had sat for almost 400 minutes
# so they may have appeared to be about 0.4 MPa drier than they really were.
# This is good to know, but I don't think it torpedoes our analysis.
# As long as two plants are more than half an MPa different,
# that distinction likely reflects real physiology.

# Would it make sense to correct the water potential estimates with this simple model?

fulldata$corrected_phi = fulldata$water_potential_m_pa - 0.001561*fulldata$minutes_in_cooler

ggplot(data = fulldata) +
  theme_cowplot() +
  geom_jitter(width = 0.25,
              aes(x = drydown_day, 
                  y = corrected_phi,
                  color = water_level,
                  shape = colonized)) +
  xlab("Days of drydown") +
  ylab("Water potential (MPa)\n(corrected for time spent in cooler)")




# Example ggplot code from an old project below:

massplot = ggplot(data = bio_and_col_onlyclean) +
  geom_boxplot(outlier.alpha = 0,
               aes(x = Fungi, y = total_biomass)) +
  geom_jitter(width = 0.20,
              aes(x = Fungi, y = total_biomass)) +
  facet_grid(. ~ N_level, labeller = labeller(N_level = labels)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Total plant biomass (g)") +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  geom_text(data = anothertry, aes(x, y, label = labs)) +
  xlab("Fungi on roots at harvest")
