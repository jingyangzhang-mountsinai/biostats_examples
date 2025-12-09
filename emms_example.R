# Author: Jingyang (Judy) Zhang
# Date: Dec. 9th, 2025
# Purpose: an example of estimating marginal means.


library(emmeans)
data("warpbreaks")

# Fit an anova model.
mod_anova <- aov(breaks ~ wool * tension, data = warpbreaks)
summary(mod_anova)
mod_lm <- lm(breaks ~ wool * tension, data = warpbreaks)
summary(mod_lm)

# Estimated marginal means for wool type. 
emmeans(mod_anova, ~wool)

# Estimated marginal means for wool at each tension level.
emmeans(mod, ~wool | tension)
