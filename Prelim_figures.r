# Making some preliminary figures for the Hebeloma data

library(tidyverse)
library(cowplot)
library(lme4)
library(lmerTest)

# Constants:

PHI_LOSS_PER_MIN = -0.001561 # from lmer below, to correct water potential
APPROX_ROOT_WET_TO_DRY_CONVERSION = 0.261 # conversion factor from small test of 10 root systems,
# in Google Sheets as "wet to dry root mass conversion"

fulldata = read_csv("Merged Hebeloma 11Apr2022 OT.csv")
shootmass = read_csv("LB Hebeloma Seedling Shoot Dry Mass.xlsx - Sheet1.csv", skip = 2)

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

fulldata$corrected_phi = fulldata$water_potential_m_pa - PHI_LOSS_PER_MIN*fulldata$minutes_in_cooler
fulldata = mutate(fulldata, ID = seedling)
fulldata = left_join(fulldata, shootmass)

fulldata$percent_col[fulldata$colonized == "N"] = 0 #NAs were making analysis difficult

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

colmodel = lmer(percent_col ~ n * water_level + (1|drydown_day),
                data = colonized)
summary(colmodel)
anova(colmodel)
# water level is significant, n level no.

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
  geom_point(aes(x = as.numeric(minutes_in_cooler), 
                 y = water_potential_m_pa,
                 color = water_level,
                 shape = colonized)) +
  geom_smooth(aes(x = as.numeric(minutes_in_cooler), 
                y = water_potential_m_pa),
              formula = y ~ x, 
              method = "lm") +
  xlab("Minutes since harvest") +
  ylab("Water potential (MPa)")

# This is encouraging. I don't see a pattern of drier plants after more time in the cooler.
# buuut adding the smooth geom really undermined that impression.

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
# Added constant and column above to do just that.


ggplot(data = fulldata) +
  theme_cowplot() +
  geom_jitter(width = 0.25,
              aes(x = drydown_day, 
                  y = corrected_phi,
                  color = water_level,
                  shape = colonized)) +
  xlab("Days of drydown") +
  ylab("Water potential (MPa)\n(corrected for time spent in cooler)")

ggplot(data = colonized) +
  theme_cowplot() +
  geom_point(aes(x = percent_col, 
                 y = corrected_phi,
                 color = water_level)) +
  geom_smooth(aes(x = percent_col, 
                  y = corrected_phi),
              formula = y ~ x, 
              method = "lm") +
  xlab("Percent colonization") +
  ylab("Corrected water potential (MPa)")

write_csv(fulldata, "data_for_picking_samples.csv")

# Looking at biomass

fulldata$`mass (g)`[fulldata$seedling == "MM037"] = 0.077 # was mis-typed originally as 0.77


fulldata = mutate(fulldata, sketchy_mass_total = root_mass_mg + `mass (g)`)
# Need to account for the fact that your root masses are wet
# and your shoot masses are dry. For now, though, I guess this is fine.

fulldata = mutate(fulldata, est_dry_root_mass = root_mass_mg * APPROX_ROOT_WET_TO_DRY_CONVERSION)
fulldata = mutate(fulldata, total_plant_mass = est_dry_root_mass + `mass (g)`)
fulldata$harv_day_cont = as.numeric(fulldata$harvest_day)


ggplot(data = fulldata) +
  theme_cowplot() +
  geom_point(aes(x = percent_col, 
                 y = total_plant_mass,
                 color = water_level, 
                 shape = n)) +
  geom_smooth(aes(x = percent_col, 
                  y = total_plant_mass),
              formula = y ~ x, 
              method = "lm") +
  xlab("Percent colonization") +
  ylab("Approximate total biomass (g)")

labels = c(Minus = "No N", Plus = "With N")
ggplot(data = fulldata) +
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(0.9),
               aes(x = water_level, y = total_plant_mass,
                   color = colonized)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.2),
             aes(x = water_level, y = total_plant_mass,
                 color = colonized)) +
  facet_grid(. ~ n, labeller = labeller(n = labels)) +
  xlab("Water level") +
  ylab("Approximate total dry biomass (g)")

massmodel = lm(total_plant_mass ~ water_level * n * colonized, data = fulldata)
summary(massmodel)
anova(massmodel)

# Water potential boxplot

labels = c(Minus = "No N", Plus = "With N")
ggplot(data = fulldata) +
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(0.9),
               aes(x = water_level, y = corrected_phi,
                   color = colonized)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.2),
             aes(x = water_level, y = corrected_phi,
                 color = colonized, alpha = harv_day_cont)) +
  facet_grid(. ~ n, labeller = labeller(n = labels)) +
  xlab("Water level") +
  ylab("Water potential (MPa)")

fulldata$corrected_phi = as.numeric(fulldata$corrected_phi)

phimodel = lm(corrected_phi ~ water_level * n * colonized * harvest_day, data = fulldata)
plot(phimodel)
summary(phimodel)
original_anova = aov(phimodel)
summary(original_anova)
TukeyHSD(original_anova)

phimodel_lme = lmer(corrected_phi ~ water_level * n * colonized +(1|harvest_day), data = fulldata)
summary(phimodel_lme)

labels = c(L = "Low water", M = "Medium water", H = "High water")
panel1 = ggplot(data = subset(fulldata, n == "Plus")) +
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(0.9),
               aes(x = drydown_day, y = corrected_phi,
                   color = colonized)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.2),
             aes(x = drydown_day, y = corrected_phi,
                 color = colonized)) +
  scale_color_manual(values = c("cadetblue3", "cadetblue4")) +
  facet_grid(. ~ water_level, labeller = labeller(water_level = labels)) +
  guides(color = "none") +
  xlab("Days since watering") +
  ylab("Water potential (MPa)") +
  theme(plot.margin = unit(c(1.5,3,1,2), "lines"))


panel2 = ggplot(data = subset(fulldata, n == "Minus")) +
  theme_cowplot() +
  theme(legend.position = "left") +
  geom_boxplot(outlier.alpha = 0,
               position = position_dodge(0.9),
               aes(x = drydown_day, y = corrected_phi,
                   color = colonized)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.9,
                                             jitter.width = 0.2),
             aes(x = drydown_day, y = corrected_phi,
                 color = colonized)) +
  scale_color_manual(name = "Ectomycorrhizal colonization   ", values = c("cadetblue3", "cadetblue4")) +
  facet_grid(. ~ water_level, labeller = labeller(water_level = labels)) +
  xlab("Days since watering") +
  ylab("Water potential (MPa)") +
  theme(plot.margin = unit(c(2,3,1,2), "lines"))

water_potential_plot_for_poster = plot_grid(panel1, panel2, 
                                            labels = c("A: With N", "B: No N"),
          ncol = 1)

save_plot("plots/water_potential_multipanel.pdf", 
          water_potential_plot_for_poster,
          base_height = 6.5)

# What can I say on the poster that is true and interesting?

fulldata$days_since_water = as.numeric(as.character(fulldata$drydown_day))

just_highwater = subset(fulldata, water_level == "H")

highwater_lm = lm(corrected_phi ~ days_since_water*colonized, data = subset(just_highwater, n = "Plus"))
summary(highwater_lm)

test = lm(corrected_phi ~ harvest_day*colonized*n, data = just_highwater)
anova(test)

test = lm(corrected_phi ~ harvest_day*colonized, data = subset(just_highwater, n = "Plus"))
myanova = aov(test)
TukeyHSD(myanova)

phimodel_again = lm(corrected_phi ~ water_level * n * colonized * days_since_water, data = fulldata)
plot(phimodel_again)
summary(phimodel_again)
original_anova_again = aov(phimodel_again)
summary(original_anova_again)
TukeyHSD(original_anova_again)

phimodel = lm(corrected_phi ~ water_level * n * colonized * harvest_day, data = fulldata)
plot(phimodel)
summary(phimodel)
original_anova = aov(phimodel)
summary(original_anova)
TukeyHSD(original_anova)

# only marginally sigificant interactions are water_level:harvest_day and n:colonized:harvest_day, so for Tukey simplicity I'm going to build a simpler model.

phimodel_simpler = lm(corrected_phi ~ water_level + n + colonized + harvest_day + water_level:harvest_day + n:colonized:harvest_day, data = fulldata)
phimodel_simpler_anova = aov(phimodel_simpler)
TukeyHSD(phimodel_simpler_anova) # still hard to interpret.

# let's eliminate all non-significant terms.

phimodel_simplest = lm(corrected_phi ~ water_level + colonized + harvest_day, data = fulldata)
simplest_anova = aov(phimodel_simplest)
TukeyHSD(simplest_anova)

# This output is pretty hard to interpret. Maybe I am more interested in whether having fungi
# changed the rate at which you dried down? That is did colonization change the slope of the line
# relation harvest day to shoot water potential?
# Actually that IS what this is showing I think. It's just not quite an ANOVA question.


ggplot(data = fulldata) +
  geom_jitter(aes(x = as.numeric(as.character(drydown_day)), y = corrected_phi,
                 color = water_level,
                 shape = colonized))

phithroughtime = lm(corrected_phi ~ harv_day_cont * colonized * water_level, data = fulldata)
summary(phithroughtime)

# Example ggplot code from an old project below (using to scavenge ggplot syntax):

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
