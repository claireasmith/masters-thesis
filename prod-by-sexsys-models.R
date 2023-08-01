# Claire Smith
# 26 July 2023
# Exploring the relationship between sex system and pollen production

library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)

### Read data
sizenum_raw <- read.csv("processed-data/size-prod-norep.csv", stringsAsFactors = T)

### Prepare individual-level data for model
sizenum <- sizenum_raw %>% 
  # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(Sex_sys, Species) %>% 
  filter(!is.na(Avg_diam) & !is.na(Avg_pol_anth))

sizenum$Log_avg_pol_anth <- log(sizenum$Avg_pol_anth)

### Get species-level data to get species' means for plotting
sizenum_sp <- sizenum %>% 
  group_by(Sex_sys, Species, source) %>% 
  summarize(Avg_pol_anth = mean(Avg_pol_anth),
            Avg_diam = mean(Avg_diam),
            N = n()) %>%
  filter(N>=10) %>% # Keep only species with at least 10 individuals
  droplevels() %>% 
  as.data.frame()

# For the mixed-effect model, I'll only take the species with at least 10 individuals measured per group. 
spp_over10 <- sizenum %>% group_by(Sex_sys, Species) %>% summarize(n=n()) %>% filter(n>=10)
keep_spp <- spp_over10$Species
sizenum_filt <- sizenum %>% group_by(Sex_sys, Species) %>% 
  filter(Species %in% keep_spp)

### Fit the model

#### Option 1: mixed effects model with random intercept for species  
prodsize_mod <- lmer(Log_avg_pol_anth ~ Avg_diam + Sex_sys + (1|Species), data = sizenum_filt)
#### Option 2: linear model with species nested within sex system
# prodsize_mod <- lm(Log_avg_pol_anth ~ Avg_diam + Sex_sys + Sex_sys:Species, data = sizenum_filt)

summary(prodsize_mod)
anova(prodsize_mod) 

### Calculate marginal means, test for pairwise differences between sex systems

# emmeans(prodsize_mod , list(pairwise ~ Sex_sys), adjust = "tukey")
emmeans(prodsize_mod , list(pairwise ~ Sex_sys), adjust = "tukey")$`pairwise differences of Sex_sys`
# plot(emmeans(prodsize_mod , ~Sex_sys))
(em_prodsize <- emmeans(prodsize_mod , ~Sex_sys) %>%  as.data.frame())

### Plot estimated marginal means for pollen production: 
ggplot(data=em_prodsize, aes(x=Sex_sys, y=emmean)) + 
  geom_point(size=3) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  # Add letters corresponding to Tukey's HSD
geom_text(aes(x=1, y=9.2), label="AB") + 
  geom_text(aes(x=2, y=7.2), label="A") + 
  geom_text(aes(x=3, y=8.4), label="B") + 
  scale_y_continuous(name="Pollen production per anther",
                     breaks = log(c(250, 500,1000, 2500, 5000, 10000, 20000)),
                     labels = c(250, 500,1000, 2500, 5000, 10000, 20000),
                     limits = c(log(250), NA)) +
  scale_x_discrete(name = "Sex system", labels = c("Dioecious", "Monoecious", "Hermaphroditic")) + 
  theme_bw()

### Plot estimated marginal means alongside raw data and species' means
ggplot(data=em_prodsize, aes(x=Sex_sys, y=emmean)) + 
  geom_point(size=3) + 
  geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL), width=0.1) +
  # Add letters corresponding to Tukey's HSD
  geom_text(aes(x=1, y=9.2), label="AB") + 
  geom_text(aes(x=2, y=7.2), label="A") + 
  geom_text(aes(x=3, y=8.4), label="B") + 
  # Raw data (individual means) with jittered points
  gghalves::geom_half_point(data=sizenum, aes(x=Sex_sys, y=log(Avg_pol_anth)),
                            ## draw on the left
                            side = "l",
                            range_scale = .4,
                            alpha = .3,
                            size=1) +
  # Species means
  gghalves::geom_half_point(data=sizenum_sp, aes(x=Sex_sys, y=log(Avg_pol_anth)),
                            ## draw on the right
                            side = "r",
                            transformation = position_identity(),
                            # transformation = position_jitter(w = 0, h = 0),  # we do NOT want vertical jitter!
                            alpha = .7,
                            size=5,
                            shape=95) + # species' means as bars
  scale_y_continuous(name="Pollen production per anther",
                     breaks = log(c(50, 100, 250, 500,1000, 2500, 5000, 10000, 20000)),
                     labels = c(50, 100, 250, 500,1000, 2500, 5000, 10000, 20000),
                     limits = c(log(50), NA)) +
  scale_x_discrete(name = "Sex system", labels = c("Dioecious", "Monoecious", "Hermaphroditic")) + 
  theme_bw()
