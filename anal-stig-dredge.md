Model selection for stigma pollen loads
================
Claire Smith
2023-06-19

``` r
# Load packages
library(tidyverse)
library(lme4) # for linear mixed effects models
library(lmerTest)
library(emmeans)
library(car)
library(knitr)

library(MuMIn) # for model selection

# Source files
source("clean-stig-01.R")
source("clean-stig-02.R")
source("theme_cs.R")
options(na.action = "na.fail")
```

### Prepare data

``` r
# Load data
stig <- read.csv("processed-data/stig-all.csv", stringsAsFactors = T)
```

``` r
stig$Date <- as.Date(stig$Date)
stig$Sex_sys <- as.factor(stig$Sex_sys)
stig$source <- as.factor(stig$source)
# I'll go through each species and run a model predicting pollen capture by flower height
# Since many likely have log-normal distributions of pollen capture, I'll create a log flw pollen variable now to # save time -- but will still check each species individually 
stigl <- stig %>% 
  mutate(Log_flw_pollen=log(Flw_pollen + 1)) %>% 
  rowwise %>% 
  mutate(Avg_dist = mean(c_across(c('D1', 'D2', 'D3', 'D4', 'D5')), na.rm=TRUE)) %>% 
  mutate(Inv_avg_dist = 1/Avg_dist) %>% 
  # Take only entries stigma length data
  filter(!is.na(Flw_pollen) & !is.na(Infl_max) & !is.na(Avg_dist) &!is.na(Stigma_length) & !is.na(Date)) %>%
  droplevels()
stigl <- stigl %>% filter((Species=="Carex hirtifolia" & Stigma_length>0.5) | Species != "Carex hirtifolia") # remove C. hirtifolia outlier with stigma length = 0.2
```

``` r
# How many species have max infl height and stigma length? 
unique(stigl$Species) #  13 species 
```

    ##  [1] Ambrosia artemisiifolia Chenopodium album       Carex communis         
    ##  [4] Carex hirtifolia        Carex pedunculata       Carex plantaginea      
    ##  [7] Carex stipata           Plantago lanceolata     Rumex acetosella       
    ## [10] Rumex crispus           Scirpus microcarpus     Schizacne purpurascens 
    ## [13] Thalictrum dioicum     
    ## 13 Levels: Ambrosia artemisiifolia Carex communis ... Thalictrum dioicum

``` r
# [1] Ambrosia artemisiifolia Chenopodium album       Carex communis          Carex hirtifolia       
#  [5] Carex pedunculata       Carex plantaginea       Carex stipata           Plantago lanceolata    
#  [9] Rumex acetosella        Rumex crispus           Scirpus microcarpus     Schizacne purpurascens 
# [13] Thalictrum dioicum

# Create unique ID for each plant within a species (might be repeat labels within each site/date category)
stigl$Ind_ID=paste(stigl$Date, stigl$Site, stigl$Plant, sep = "-")

#How many unique individuals per species? 
stigl %>% group_by(Species) %>%
  distinct(Plant) %>% 
  summarize(n=n())
```

    ## # A tibble: 13 × 2
    ##    Species                     n
    ##    <fct>                   <int>
    ##  1 Ambrosia artemisiifolia    30
    ##  2 Carex communis             22
    ##  3 Carex hirtifolia           25
    ##  4 Carex pedunculata          15
    ##  5 Carex plantaginea          19
    ##  6 Carex stipata              30
    ##  7 Chenopodium album          29
    ##  8 Plantago lanceolata        30
    ##  9 Rumex acetosella           30
    ## 10 Rumex crispus              25
    ## 11 Schizacne purpurascens     30
    ## 12 Scirpus microcarpus        30
    ## 13 Thalictrum dioicum         20

``` r
# #How many total flowers? There are multiple flowers measured per plant
stigl %>% group_by(Species) %>%
  summarize(n=n())
```

    ## # A tibble: 13 × 2
    ##    Species                     n
    ##    <fct>                   <int>
    ##  1 Ambrosia artemisiifolia   233
    ##  2 Carex communis            201
    ##  3 Carex hirtifolia          206
    ##  4 Carex pedunculata          80
    ##  5 Carex plantaginea          96
    ##  6 Carex stipata             320
    ##  7 Chenopodium album         149
    ##  8 Plantago lanceolata       272
    ##  9 Rumex acetosella          261
    ## 10 Rumex crispus             179
    ## 11 Schizacne purpurascens    246
    ## 12 Scirpus microcarpus       268
    ## 13 Thalictrum dioicum        110

``` r
species <- levels(stigl$Species) # list of all species in data
# for (spec in species){print(spec)}
# Filter for each species
spec = species[13]
dat <- stigl %>% filter(Species==spec) %>% 
  as_tibble() %>% 
  # Use "scale" to standardize predictor variables to z-scores
  mutate(across(c("Infl_max","Stigma_length","Date", "Inv_avg_dist", "Avg_dist"), scale)) %>% 
  mutate(across(c("Infl_max","Stigma_length","Date", "Inv_avg_dist", "Avg_dist"), as.vector)) %>%
  select(Species, Plant, Flower, Flw_pollen, Log_flw_pollen, Infl_max, Stigma_length, Date, Inv_avg_dist, Avg_dist)

#Take a look at the data
hist(dat$Flw_pollen) 
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-1.png)<!-- -->

``` r
hist(dat$Log_flw_pollen) 
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-2.png)<!-- -->

``` r
hist(dat$Infl_max) 
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-3.png)<!-- -->

``` r
hist(dat$Stigma_length)
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-4.png)<!-- -->

``` r
hist(dat$Avg_dist)
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-5.png)<!-- -->

``` r
hist(dat$Inv_avg_dist)
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-6.png)<!-- -->

``` r
table(dat$Date)
```

    ## 
    ##  -1.21775701398301 -0.712272970442717 -0.206788926902428  0.298695116637862 
    ##                 20                 20                 25                 25 
    ##   1.81514724725873 
    ##                 20

``` r
ggplot(data=dat, aes(x=Infl_max, y=Log_flw_pollen)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Plant height") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
    geom_smooth(method="lm", se=F) + 
  theme_cs()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-dredge_files/figure-gfm/fit%20models-7.png)<!-- -->

``` r
ggplot(data=dat, aes(x=Stigma_length, y=Log_flw_pollen)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Stigma length") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  geom_smooth(method="lm", se=F) + 
  theme_cs()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-dredge_files/figure-gfm/fit%20models-8.png)<!-- -->

``` r
ggplot(data=dat, aes(x=Inv_avg_dist, y=Log_flw_pollen)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  geom_smooth(method="lm", se=F) + 
  theme_cs()
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-dredge_files/figure-gfm/fit%20models-9.png)<!-- -->

``` r
ggplot(data=dat, aes(x=Date, y=Log_flw_pollen)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Date") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  geom_smooth(method="lm", se=F) + 
  theme_cs() # dates are in number format
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-dredge_files/figure-gfm/fit%20models-10.png)<!-- -->

``` r
# Fit full model
dat_mod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=dat)
summary(dat_mod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: dat
    ## 
    ## REML criterion at convergence: 242.3
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7575 -0.5647 -0.0172  0.6252  2.1184 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.1236   0.3515  
    ##  Residual             0.4037   0.6354  
    ## Number of obs: 110, groups:  Plant, 20
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)    1.99623    0.10001 12.17903  19.960 1.13e-10 ***
    ## Infl_max       0.14445    0.09856 12.61341   1.466   0.1672    
    ## Stigma_length  0.18562    0.07384 89.94800   2.514   0.0137 *  
    ## Date           0.28539    0.06248 89.75537   4.568 1.56e-05 ***
    ## Inv_avg_dist   0.03448    0.10263 11.60809   0.336   0.7429    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.036                     
    ## Stigm_lngth  0.021  0.011              
    ## Date         0.017 -0.016 -0.029       
    ## Inv_avg_dst  0.027 -0.025 -0.060 -0.016

``` r
anova(dat_mod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##               Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Infl_max      0.8672  0.8672     1 12.613  2.1482   0.16722    
    ## Stigma_length 2.5507  2.5507     1 89.948  6.3184   0.01373 *  
    ## Date          8.4232  8.4232     1 89.755 20.8655 1.564e-05 ***
    ## Inv_avg_dist  0.0456  0.0456     1 11.608  0.1129   0.74291    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Model selection
dat_mod_dredge <- MuMIn::dredge(dat_mod)
```

    ## Warning in MuMIn::dredge(dat_mod): comparing models fitted by REML

    ## Fixed term is "(Intercept)"

``` r
head(dat_mod_dredge) # Take a look
```

    ## Global model call: lmer(formula = Log_flw_pollen ~ Infl_max + Stigma_length + Date + 
    ##     Inv_avg_dist + (1 | Plant), data = dat)
    ## ---
    ## Model selection table 
    ##    (Int)    Dat Inf_max Inv_avg_dst Stg_lng df   logLik  AICc delta weight
    ## 2  1.983 0.2915                              4 -120.108 248.6  0.00  0.459
    ## 10 1.990 0.2872                      0.1855  5 -119.529 249.6  1.04  0.273
    ## 4  1.986 0.2893  0.1384                      5 -120.214 251.0  2.41  0.138
    ## 12 1.994 0.2857  0.1446              0.1790  6 -119.825 252.5  3.87  0.066
    ## 6  1.984 0.2909             0.04886          5 -121.491 253.6  4.96  0.038
    ## 14 1.991 0.2869             0.03803  0.1911  6 -120.805 254.4  5.83  0.025
    ## Models ranked by AICc(x) 
    ## Random terms (all models): 
    ##   1 | Plant

``` r
results <- get.models(dat_mod_dredge, subset = delta <= 2) # Top models have lowest AICc, models within 2 AICc are indistinguishable
results
```

    ## $`2`
    ## Linear mixed model fit by REML ['lmerModLmerTest']
    ## Formula: Log_flw_pollen ~ Date + (1 | Plant)
    ##    Data: dat
    ## REML criterion at convergence: 240.2155
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  Plant    (Intercept) 0.2266  
    ##  Residual             0.6756  
    ## Number of obs: 110, groups:  Plant, 20
    ## Fixed Effects:
    ## (Intercept)         Date  
    ##      1.9832       0.2915  
    ## 
    ## $`10`
    ## Linear mixed model fit by REML ['lmerModLmerTest']
    ## Formula: Log_flw_pollen ~ Date + Stigma_length + (1 | Plant)
    ##    Data: dat
    ## REML criterion at convergence: 239.0589
    ## Random effects:
    ##  Groups   Name        Std.Dev.
    ##  Plant    (Intercept) 0.3507  
    ##  Residual             0.6365  
    ## Number of obs: 110, groups:  Plant, 20
    ## Fixed Effects:
    ##   (Intercept)           Date  Stigma_length  
    ##        1.9898         0.2872         0.1855  
    ## 
    ## attr(,"rank")
    ## function (x) 
    ## do.call("rank", list(x))
    ## <environment: 0x000002c48e763720>
    ## attr(,"call")
    ## AICc(x)
    ## attr(,"class")
    ## [1] "function"     "rankFunction"
    ## attr(,"beta")
    ## [1] "none"

``` r
## If just one top model: 
# Get summary of top model results
topmod <- get.models(dat_mod_dredge, subset = delta<=2)[[2]]
modsum <- summary(topmod)
modsum
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Date + Stigma_length + (1 | Plant)
    ##    Data: dat
    ## 
    ## REML criterion at convergence: 239.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.70069 -0.59501  0.02313  0.59528  2.01312 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.1230   0.3507  
    ##  Residual             0.4052   0.6365  
    ## Number of obs: 110, groups:  Plant, 20
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)    1.98983    0.09982 13.67409  19.935 1.69e-11 ***
    ## Date           0.28722    0.06257 89.13555   4.591 1.44e-05 ***
    ## Stigma_length  0.18548    0.07379 95.03156   2.514   0.0136 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Date  
    ## Date         0.018       
    ## Stigm_lngth  0.022 -0.030

``` r
anova(topmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##               Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Date          8.5380  8.5380     1 89.136 21.0735 1.442e-05 ***
    ## Stigma_length 2.5603  2.5603     1 95.032  6.3193   0.01362 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
## If more than one top model: 
# Average model results, if needed:
# summary(model.avg(dat_mod_dredge, subset = delta <= 2))

#normality of residuals
hist(residuals(topmod)) 
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-11.png)<!-- -->

``` r
qqnorm(residuals(topmod))
qqline(residuals(topmod))
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-12.png)<!-- -->

``` r
# homogeneity
plot(topmod)
```

![](anal-stig-dredge_files/figure-gfm/fit%20models-13.png)<!-- -->

``` r
# No obvious deviations from normality - variance in residuals does seem to decrease with increased fitted values. 
```

``` r
options(na.action = "na.omit") # put back NA settings the way they were...
```
