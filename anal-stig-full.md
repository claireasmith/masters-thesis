Full pollen load model, including stig length
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
# Source files
source("clean-stig-01.R")
source("clean-stig-02.R")
source("theme_cs.R")
# Load data
stig <- read.csv("processed-data/stig-all.csv", stringsAsFactors = T)
stig$Date <- as.Date(stig$Date)
stig$Sex_sys <- as.factor(stig$Sex_sys)
stig$source <- as.factor(stig$source)
summary(stig)
```

    ##                 Species          Date                 Site          Plant     
    ##  Thalictrum dioicum :1046   Min.   :2001-07-09   BI     : 666   N14    : 135  
    ##  Festuca campestris : 748   1st Qu.:2002-07-09   FBF    : 297   N10    : 130  
    ##  Phleum pratense    : 534   Median :2004-05-28   CT     : 255   N11    : 130  
    ##  Plantago lanceolata: 404   Mean   :2007-10-29   KMW    : 113   N13    : 127  
    ##  Bromus inermis     : 391   3rd Qu.:2004-09-03   KME    : 112   N12    : 125  
    ##  Rumex acetosella   : 378   Max.   :2021-09-20   (Other): 292   N18    : 116  
    ##  (Other)            :3252   NA's   :93           NA's   :5018   (Other):5990  
    ##      Flower       Flw_pollen      Stigmas_per_flw Flower_height   
    ##  2      :1955   Min.   :   0.00   Min.   :1.000   Min.   :  4.00  
    ##  1      :1867   1st Qu.:   2.00   1st Qu.:1.000   1st Qu.: 21.50  
    ##  3      :1458   Median :   8.00   Median :2.000   Median : 39.00  
    ##  A      : 201   Mean   :  35.67   Mean   :1.747   Mean   : 40.83  
    ##  B      : 183   3rd Qu.:  38.00   3rd Qu.:2.000   3rd Qu.: 56.00  
    ##  (Other):1085   Max.   :2200.00   Max.   :4.000   Max.   :126.00  
    ##  NA's   :   4   NA's   :21        NA's   :2788    NA's   :6104    
    ##     Infl_max         Infl_min            D1               D2         
    ##  Min.   :  5.00   Min.   :  2.50   Min.   :  0.50   Min.   :   1.00  
    ##  1st Qu.: 32.38   1st Qu.: 38.00   1st Qu.: 12.00   1st Qu.:  17.00  
    ##  Median : 53.50   Median : 52.30   Median : 23.00   Median :  32.00  
    ##  Mean   : 56.82   Mean   : 56.42   Mean   : 46.26   Mean   :  61.59  
    ##  3rd Qu.: 75.50   3rd Qu.: 71.00   3rd Qu.: 46.00   3rd Qu.:  66.00  
    ##  Max.   :158.00   Max.   :108.50   Max.   :950.00   Max.   :1000.00  
    ##  NA's   :509      NA's   :4111     NA's   :200      NA's   :200      
    ##        D3                D4                D5          Stigma_length  
    ##  Min.   :   3.00   Min.   :   2.50   Min.   :   2.00   Min.   :0.200  
    ##  1st Qu.:  20.00   1st Qu.:  25.00   1st Qu.:  30.00   1st Qu.:0.950  
    ##  Median :  39.00   Median :  48.00   Median :  53.50   Median :1.800  
    ##  Mean   :  73.27   Mean   :  88.03   Mean   :  95.47   Mean   :1.898  
    ##  3rd Qu.:  78.00   3rd Qu.: 101.00   3rd Qu.: 111.00   3rd Qu.:2.500  
    ##  Max.   :1000.00   Max.   :1500.00   Max.   :1000.00   Max.   :6.500  
    ##  NA's   :209       NA's   :211       NA's   :249       NA's   :3966   
    ##     source               Sex_sys    
    ##  CS2021:1735   dioecious     :1424  
    ##  JF2001:2230   hermaphroditic:3421  
    ##  JF2004:2788   monoecious    :1908  
    ##                                     
    ##                                     
    ##                                     
    ## 

I think stigma length is probably an important covariate for determining
pollen load. I’ll pull out only the species that have stigma length
measurements. I’ll run a full model to predict pollen load: pollen load
~ stigma length + density + height. These will all be JF species so the
measure for height will be inflorescence max height.

### Prepare data

``` r
# I'll go through each species and run a model predicting pollen capture by flower height

# Since many likely have log-normal distributions of pollen capture, I'll create a log flw pollen variable now to # save time -- but will still check each species individually 
stigl <- stig %>% 
  mutate(Log_flw_pollen=log(Flw_pollen + 1)) %>% 
  rowwise %>% 
  mutate(Avg_dist = mean(c_across(c('D1', 'D2', 'D3', 'D4', 'D5')), na.rm=TRUE)) %>% 
  # Take only entries stigma length data
  filter(!is.na(Flw_pollen) & !is.na(Infl_max) & !is.na(Avg_dist) &!is.na(Stigma_length) & !is.na(Date))
```

### Ambrosia artemisiifolia (n=30)

``` r
#model for height
AMARmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=AMAR)
```

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

``` r
summary(AMARmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: AMAR
    ## 
    ## REML criterion at convergence: 623.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -3.09737 -0.61356  0.05233  0.70478  2.27772 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.2166   0.4654  
    ##  Residual             0.6935   0.8328  
    ## Number of obs: 233, groups:  Plant, 30
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)
    ## (Intercept)   -1.491e+02  2.842e+02  2.124e+02  -0.525    0.600
    ## Infl_max      -7.064e-04  3.969e-03  2.389e+01  -0.178    0.860
    ## Stigma_length  2.181e-01  1.377e-01  2.085e+02   1.585    0.115
    ## Date           1.209e-02  2.246e-02  2.124e+02   0.538    0.591
    ## Inv_avg_dist   4.668e+00  4.345e+00  2.430e+01   1.074    0.293
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.029                     
    ## Stigm_lngth  0.438  0.041              
    ## Date        -1.000 -0.030 -0.438       
    ## Inv_avg_dst -0.058 -0.152 -0.136  0.058
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling

``` r
kable(anova(AMARmod))
```

|               |    Sum Sq |   Mean Sq | NumDF |     DenDF |   F value |   Pr(\>F) |
|:--------------|----------:|----------:|------:|----------:|----------:|----------:|
| Infl_max      | 0.0219603 | 0.0219603 |     1 |  23.88595 | 0.0316654 | 0.8602647 |
| Stigma_length | 1.7412397 | 1.7412397 |     1 | 208.54791 | 2.5107575 | 0.1145868 |
| Date          | 0.2008478 | 0.2008478 |     1 | 212.40695 | 0.2896098 | 0.5910341 |
| Inv_avg_dist  | 0.8004318 | 0.8004318 |     1 |  24.30424 | 1.1541720 | 0.2932282 |

``` r
# No significant effect of density or max inflorescence height - significant effect of stigma length. Longer stigmas capture more pollen. 
# plot(AMARmod)
```

No significant effect of anything on pollen capture

``` r
AMARmod.fit <- AMAR
fit.AMARmod <- fitted(AMARmod)
AMARmod.fit <- mutate(AMARmod.fit, fit.AMARmod)
emmsAMAR <- emmip(AMARmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

Fitted values with line of best fit:

``` r
# emmip_ggplot(emms) #plots predicted slope
# ggplot() + 
#   geom_point(data = AMARmod.fit, aes(x=Stigma_length, y=fit.AMARmod)) + 
#   geom_smooth(data=AMARmod.fit, aes(x=Stigma_length, y=fit.AMARmod),
#               method="lm",
#               formula = y ~ x, linetype=1, colour="blue", se=F, linewidth=0.7) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(log(10)+1,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)",
#                      breaks = log(c(0.5, 1, 2, 5)),
#                      labels = c(0.5, 1, 2, 5),
#                      limits=c(0.5,5)) +
#   theme_cs()
```

Raw data with emmip slope:

``` r
# ggplot() + 
#   geom_point(data = AMAR, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsAMAR, aes(xvar, yvar), linewidth=0.7) + # emmip slope
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   # scale_x_continuous(name = "Stigma length (cm)",
#   #                    breaks = log(c(0, 0.1, 0.2, 0.5 , 1, 2, 5)),
#   #                    labels = c(0, 0.1, 0.2, 0.5, 1, 2, 5)) + 
#     scale_x_continuous(name = "Stigma length (cm)") + 
# 
#                      # breaks = log(c(0.5, 1, 2, 5)),
#                      # labels = c(0.5, 1, 2, 5),
#                      # limits=c(0.5,5)) +
#   theme_cs()
# 
# # max(AMAR$Stigma_length)
```

### Carex communis (n=22)

``` r
CACO <- stigl %>% 
  filter(Species=="Carex communis") %>% 
  as.data.frame()
# summary(CACO) 
CACO$Inv_avg_dist <- 1/(CACO$Avg_dist)
```

``` r
# fit full model with all terms
CACOmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=CACO)
summary(CACOmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: CACO
    ## 
    ## REML criterion at convergence: 540.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8123 -0.7090 -0.0014  0.6028  2.6241 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.1713   0.4139  
    ##  Residual             0.7532   0.8678  
    ## Number of obs: 201, groups:  Plant, 22
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -3.973e+03  4.538e+02  1.850e+02  -8.755 1.27e-15 ***
    ## Infl_max      -1.866e-02  3.719e-02  1.952e+01  -0.502  0.62142    
    ## Stigma_length  4.591e-01  1.490e-01  1.942e+02   3.082  0.00235 ** 
    ## Date           3.167e-01  3.618e-02  1.850e+02   8.755 1.27e-15 ***
    ## Inv_avg_dist   6.326e+00  1.223e+01  1.892e+01   0.517  0.61089    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.024                     
    ## Stigm_lngth  0.173 -0.007              
    ## Date        -1.000 -0.025 -0.174       
    ## Inv_avg_dst -0.041  0.015  0.020  0.040

``` r
kable(anova(CACOmod))
```

|               |     Sum Sq |    Mean Sq | NumDF |     DenDF |    F value |   Pr(\>F) |
|:--------------|-----------:|-----------:|------:|----------:|-----------:|----------:|
| Infl_max      |  0.0625773 |  0.0625773 |     1 |  19.49623 |  0.0830873 | 0.7761985 |
| Stigma_length | 14.4001727 | 14.4001727 |     1 | 193.73236 | 19.1198881 | 0.0000200 |
| Date          |  0.6986871 |  0.6986871 |     1 |  27.54485 |  0.9276847 | 0.3438453 |
| Inv_avg_dist  |  0.0373594 |  0.0373594 |     1 |  18.87446 |  0.0496041 | 0.8261461 |

Significant effect of stigma length on pollen capture.

``` r
CACOmod.fit <- CACO
fit.CACOmod <- fitted(CACOmod)
CACOmod.fit <- mutate(CACOmod.fit, fit.CACOmod)
emmsCACO <- emmip(CACOmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CACOmod.fit, aes(x=Stigma_length, y=fit.CACOmod)) + 
  geom_smooth(data=CACOmod.fit, aes(x=Stigma_length, y=fit.CACOmod),
              method="lm",
              formula = y ~ x, linetype=1, colour="blue", se=F, linewidth=0.7) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CACO%20stig%20length%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = CACO, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsCACO, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = CACO, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Carex hirtifolia (n=25)

``` r
CAHI <- stigl %>% 
  filter(Species=="Carex hirtifolia") %>% 
  as.data.frame()
# summary(CAHI) 
CAHI$Inv_avg_dist <- 1/(CAHI$Avg_dist)
```

``` r
# fit full model with all terms
CAHImod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=CAHI)
summary(CAHImod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: CAHI
    ## 
    ## REML criterion at convergence: 692
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.02736 -0.66031 -0.09946  0.54818  3.07024 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.0461   0.2147  
    ##  Residual             1.6416   1.2812  
    ## Number of obs: 206, groups:  Plant, 25
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)   
    ## (Intercept)   -728.62812 1102.20480  191.59372  -0.661  0.50937   
    ## Infl_max         0.05673    0.04135   21.17098   1.372  0.18448   
    ## Stigma_length    0.80908    0.28891  192.41850   2.800  0.00562 **
    ## Date             0.05797    0.08783  191.60151   0.660  0.51003   
    ## Inv_avg_dist     6.26855    8.00363   20.17346   0.783  0.44260   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.024                     
    ## Stigm_lngth  0.279 -0.050              
    ## Date        -1.000 -0.024 -0.280       
    ## Inv_avg_dst -0.103  0.186 -0.004  0.103

``` r
# anova(CAHImod)
kable(anova(CAHImod))
```

|               |     Sum Sq |    Mean Sq | NumDF |     DenDF |   F value |   Pr(\>F) |
|:--------------|-----------:|-----------:|------:|----------:|----------:|----------:|
| Infl_max      |  3.1629189 |  3.1629189 |     1 |  21.15254 | 1.9267574 | 0.1795615 |
| Stigma_length | 15.8272733 | 15.8272733 |     1 | 193.05616 | 9.6415107 | 0.0021886 |
| Date          |  2.1193655 |  2.1193655 |     1 |  38.98982 | 1.2910553 | 0.2627894 |
| Inv_avg_dist  |  0.8541586 |  0.8541586 |     1 |  19.88517 | 0.5203283 | 0.4790978 |

Significant effect of stigma length on pollen capture.

``` r
CAHImod.fit <- CAHI
fit.CAHImod <- fitted(CAHImod)
CAHImod.fit <- mutate(CAHImod.fit, fit.CAHImod)
emmsCAHI <- emmip(CAHImod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CAHImod.fit, aes(x=Stigma_length, y=fit.CAHImod)) + 
  geom_smooth(data=CAHImod.fit, aes(x=Stigma_length, y=fit.CAHImod),
              method="lm",
              formula = y ~ x, linetype=1, colour="blue", se=F, linewidth=0.7) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CAHI%20stig%20length%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = CAHI, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsCAHI, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = CAHI, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Carex pedunculata (n=15)

``` r
CAPE <- stigl %>% 
  filter(Species=="Carex pedunculata") %>% 
  as.data.frame()
# summary(CAPE) 
CAPE$Inv_avg_dist <- 1/(CAPE$Avg_dist)
```

``` r
# fit full model with all terms
CAPEmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=CAPE)
summary(CAPEmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: CAPE
    ## 
    ## REML criterion at convergence: 205
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2177 -0.4934 -0.1183  0.6099  2.2360 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.4132   0.6428  
    ##  Residual             0.6070   0.7791  
    ## Number of obs: 80, groups:  Plant, 15
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -5.854e+03  8.976e+02  6.847e+01  -6.522 9.93e-09 ***
    ## Infl_max      -2.892e-02  8.132e-02  9.254e+00  -0.356    0.730    
    ## Stigma_length  1.332e-01  2.002e-01  7.075e+01   0.666    0.508    
    ## Date           4.671e-01  7.165e-02  6.848e+01   6.520 1.00e-08 ***
    ## Inv_avg_dist   1.426e+01  3.460e+01  1.082e+01   0.412    0.688    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.076                     
    ## Stigm_lngth  0.567 -0.007              
    ## Date        -1.000 -0.077 -0.568       
    ## Inv_avg_dst -0.049 -0.145 -0.049  0.048

``` r
# anova(CAPEmod)
kable(anova(CAPEmod))
```

|               |     Sum Sq |    Mean Sq | NumDF |    DenDF |    F value |   Pr(\>F) |
|:--------------|-----------:|-----------:|------:|---------:|-----------:|----------:|
| Infl_max      |  0.0022283 |  0.0022283 |     1 |  9.17922 |  0.0036708 | 0.9529870 |
| Stigma_length | 15.0011272 | 15.0011272 |     1 | 70.16413 | 24.7117351 | 0.0000045 |
| Date          |  0.1326355 |  0.1326355 |     1 | 13.47623 |  0.2184938 | 0.6476585 |
| Inv_avg_dist  |  0.0301563 |  0.0301563 |     1 | 10.78844 |  0.0496773 | 0.8277894 |

Significant effect of stigma length on pollen capture.

``` r
CAPEmod.fit <- CAPE
fit.CAPEmod <- fitted(CAPEmod)
CAPEmod.fit <- mutate(CAPEmod.fit, fit.CAPEmod)
emmsCAPE <- emmip(CAPEmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CAPEmod.fit, aes(x=Stigma_length, y=fit.CAPEmod)) + 
  # geom_smooth(data=CAPEmod.fit, aes(x=Stigma_length, y=fit.CAPEmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="blue", se=F, linewidth=0.7) + 
  geom_line(data= emmsCAPE, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CAPE%20stig%20length%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = CAPE, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsCAPE, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = CAPE, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Carex plantaginea (n=19)

``` r
CAPL <- stigl %>% 
  filter(Species=="Carex plantaginea") %>% 
  as.data.frame()
# summary(CAPL) 
CAPL$Inv_avg_dist <- 1/(CAPL$Avg_dist)
```

``` r
# fit full model with all terms
CAPLmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=CAPL)
summary(CAPLmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: CAPL
    ## 
    ## REML criterion at convergence: 251
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.7375 -0.5978 -0.1032  0.4501  2.8760 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.3285   0.5732  
    ##  Residual             0.6050   0.7778  
    ## Number of obs: 96, groups:  Plant, 19
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)   
    ## (Intercept)   -644.52863  213.69582   85.43338  -3.016  0.00337 **
    ## Infl_max         0.02824    0.04644   16.45280   0.608  0.55148   
    ## Stigma_length    0.15033    0.11720   90.53656   1.283  0.20289   
    ## Date             0.05143    0.01704   85.43269   3.019  0.00334 **
    ## Inv_avg_dist   -12.63261   14.61086   18.13878  -0.865  0.39855   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max    -0.033                     
    ## Stigm_lngth  0.030 -0.181              
    ## Date        -1.000  0.029 -0.032       
    ## Inv_avg_dst -0.051 -0.077  0.272  0.050

``` r
# anova(CAPLmod)
kable(anova(CAPLmod))
```

|               |    Sum Sq |   Mean Sq | NumDF |    DenDF |   F value |   Pr(\>F) |
|:--------------|----------:|----------:|------:|---------:|----------:|----------:|
| Infl_max      | 0.2236608 | 0.2236608 |     1 | 16.45280 | 0.3696757 | 0.5514825 |
| Stigma_length | 0.9953538 | 0.9953538 |     1 | 90.53656 | 1.6451611 | 0.2028929 |
| Date          | 5.5140976 | 5.5140976 |     1 | 85.43269 | 9.1139238 | 0.0033438 |
| Inv_avg_dist  | 0.4522760 | 0.4522760 |     1 | 18.13878 | 0.7475400 | 0.3985465 |

Significant effect of date on pollen capture.

``` r
CAPLmod.fit <- CAPL
fit.CAPLmod <- fitted(CAPLmod)
CAPLmod.fit <- mutate(CAPLmod.fit, fit.CAPLmod)
emmsCAPL <- emmip(CAPLmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)

CAPLmod.fit2 <- CAPL
fit.CAPLmod2 <- fitted(CAPLmod)
CAPLmod.fit2 <- mutate(CAPLmod.fit2, fit.CAPLmod2)
emmsCAPL2 <- emmip(CAPLmod, ~ Date, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CAPLmod.fit2, aes(x=Date, y=fit.CAPLmod2)) + 
  # geom_smooth(data=CAPLmod.fit2, aes(x=Date, y=fit.CAPLmod2),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsCAPL2, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Date", labels=as.character(unique(CAPLmod.fit2$Date)),
                     breaks=unique(CAPLmod.fit2$Date)) +
  theme_cs() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](anal-stig-full_files/figure-gfm/CAPL%20date%20fitted%20values%20plot-1.png)<!-- -->

``` r
# # emmip_ggplot(emms) #plots predicted slope
# ggplot() + 
#   geom_point(data = CAPLmod.fit, aes(x=Stigma_length, y=fit.CAPLmod)) + 
#   # geom_smooth(data=CAPLmod.fit, aes(x=Stigma_length, y=fit.CAPLmod),
#   #             method="lm",
#   #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
#   geom_line(data= emmsCAPL, aes(xvar, yvar), linewidth=0.7, colour="blue") +
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
#   scale_x_continuous(name = "Stigma length (cm)") +
#   theme_cs()
```

``` r
# ggplot() + 
#   geom_point(data = CAPL, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsCAPL, aes(xvar, yvar), linewidth=0.7, colour="blue") +
#   # geom_smooth(method="lm", data = CAPL, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Carex stipata (n=30)

``` r
CAST <- stigl %>% 
  filter(Species=="Carex stipata") %>% 
  as.data.frame()
# summary(CAST) 
```

``` r
# fit full model with all terms
CASTmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist + (1|Plant), data=CAST)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
summary(CASTmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist +  
    ##     (1 | Plant)
    ##    Data: CAST
    ## 
    ## REML criterion at convergence: 899.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8678 -0.6668 -0.1696  0.4190  3.8804 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.0000   0.0000  
    ##  Residual             0.8982   0.9478  
    ## Number of obs: 320, groups:  Plant, 30
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -2.548e+03  5.479e+02  3.150e+02  -4.651 4.86e-06 ***
    ## Infl_max       2.093e-02  9.539e-03  3.150e+02   2.194   0.0289 *  
    ## Stigma_length  7.331e-01  1.736e-01  3.150e+02   4.222 3.17e-05 ***
    ## Date           2.027e-01  4.360e-02  3.150e+02   4.650 4.90e-06 ***
    ## Avg_dist      -1.751e-03  7.093e-04  3.150e+02  -2.469   0.0141 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max    -0.028                     
    ## Stigm_lngth  0.170 -0.194              
    ## Date        -1.000  0.028 -0.170       
    ## Avg_dist     0.018  0.220 -0.030 -0.018
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
# anova(CASTmod)
kable(anova(CASTmod))
```

|               |    Sum Sq |   Mean Sq | NumDF |    DenDF |   F value |   Pr(\>F) |
|:--------------|----------:|----------:|------:|---------:|----------:|----------:|
| Infl_max      |  3.817916 |  3.817916 |     1 | 315.0000 |  4.250434 | 0.0400604 |
| Stigma_length | 23.235518 | 23.235518 |     1 | 315.0000 | 25.867787 | 0.0000006 |
| Date          |  7.200490 |  7.200490 |     1 | 314.9957 |  8.016208 | 0.0049343 |
| Avg_dist      |  5.474096 |  5.474096 |     1 | 315.0000 |  6.094237 | 0.0140928 |

boundary (singular) fit: see help(‘isSingular’)… Significant effect of
height, date, stigma length, and avg distance - but is it valid?

``` r
CASTmod.fit <- CAST
fit.CASTmod <- fitted(CASTmod)
CASTmod.fit <- mutate(CASTmod.fit, fit.CASTmod)
emmsCAST_stiglength <- emmip(CASTmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
emmsCAST_date <- emmip(CASTmod, ~ Date, cov.reduce=range,
              type="response", plotit=F)
emmsCAST_avgdist <- emmip(CASTmod, ~ Avg_dist, cov.reduce=range,
              type="response", plotit=F)
emmsCAST_inflmax <- emmip(CASTmod, ~ Infl_max, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CASTmod.fit, aes(x=Stigma_length, y=fit.CASTmod)) + 
  # geom_smooth(data=CASTmod.fit, aes(x=Stigma_length, y=fit.CASTmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsCAST_stiglength, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CAST%20stig%20length%20fitted%20values%20plot-1.png)<!-- -->

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CASTmod.fit, aes(x=Date, y=fit.CASTmod)) + 
  # geom_smooth(data=CASTmod.fit, aes(x=Stigma_length, y=fit.CASTmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsCAST_date, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Date", labels=as.character(unique(CASTmod.fit$Date)),
                     breaks=unique(CASTmod.fit$Date)) +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CAST%20date%20fitted%20values%20plot-1.png)<!-- -->

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CASTmod.fit, aes(x=Infl_max, y=fit.CASTmod)) + 
  # geom_smooth(data=CASTmod.fit, aes(x=Stigma_length, y=fit.CASTmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsCAST_inflmax, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Plant height") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CAST%20infl%20max%20fitted%20values%20plot-1.png)<!-- -->

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CASTmod.fit, aes(x=Avg_dist, y=fit.CASTmod)) + 
  # geom_smooth(data=CASTmod.fit, aes(x=Stigma_length, y=fit.CASTmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsCAST_avgdist, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Average distance to 5 nearest pollen-producing neighbours") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CAST%20avg%20dist%20fitted%20values%20plot-1.png)<!-- -->

``` r
ggplot() + 
  geom_point(data = CAST, aes(x=Stigma_length, y=Log_flw_pollen)) + 
  geom_line(data= emmsCAST_stiglength, aes(xvar, yvar), linewidth=0.7) +
  # geom_smooth(method="lm", colour="blue", data = CAST, aes(x=Stigma_length, y=Log_flw_pollen)) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
                     limits=c(0,NA)) +
  scale_x_continuous(name = "Stigma length (cm)") + 

    # scale_x_continuous(name = "Stigma length (cm)",
    #                  breaks = log(c(0.5, 1, 2, 5)),
    #                  labels = c(0.5, 1, 2, 5),
    #                  limits=c(log(0.5),log(5))) +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CAST%20stig%20length%20fitted%20slope%20raw%20values%20plot-1.png)<!-- -->

### Chenopodium album (n=29)

``` r
CHAL <- stigl %>% 
  filter(Species=="Chenopodium album") %>% 
  as.data.frame()
# summary(CHAL) 
```

``` r
# fit full model with all terms
CHALmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist + (1|Plant), data=CHAL)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
summary(CHALmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist +  
    ##     (1 | Plant)
    ##    Data: CHAL
    ## 
    ## REML criterion at convergence: 368.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.2645 -0.8235 -0.0034  0.6779  2.3664 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.0000   0.000   
    ##  Residual             0.6069   0.779   
    ## Number of obs: 149, groups:  Plant, 29
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   305.294823 659.214320 144.000000   0.463    0.644    
    ## Infl_max       -0.003921   0.003795 144.000000  -1.033    0.303    
    ## Stigma_length   2.591149   0.574056 144.000000   4.514 1.31e-05 ***
    ## Date           -0.024178   0.052228 144.000000  -0.463    0.644    
    ## Avg_dist        0.003818   0.004468 144.000000   0.855    0.394    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max    -0.153                     
    ## Stigm_lngth  0.046  0.162              
    ## Date        -1.000  0.153 -0.047       
    ## Avg_dist    -0.088 -0.075  0.096  0.087
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
# anova(CHALmod)
kable(anova(CHALmod))
```

|               |     Sum Sq |    Mean Sq | NumDF |    DenDF |    F value |   Pr(\>F) |
|:--------------|-----------:|-----------:|------:|---------:|-----------:|----------:|
| Infl_max      |  0.5736595 |  0.5736595 |     1 | 144.0000 |  0.9452214 | 0.3325689 |
| Stigma_length | 12.2745943 | 12.2745943 |     1 | 144.0000 | 20.2249068 | 0.0000141 |
| Date          |  0.0002553 |  0.0002553 |     1 | 144.0045 |  0.0004206 | 0.9836654 |
| Avg_dist      |  0.4905378 |  0.4905378 |     1 | 144.0000 |  0.8082614 | 0.3701356 |

boundary (singular) fit: see help(‘isSingular’)…

Significant effect of stigma length - but not sure it’s valid

``` r
CHALmod.fit <- CHAL
fit.CHALmod <- fitted(CHALmod)
CHALmod.fit <- mutate(CHALmod.fit, fit.CHALmod)
emmsCHAL <- emmip(CHALmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = CHALmod.fit, aes(x=Stigma_length, y=fit.CHALmod)) + 
  # geom_smooth(data=CHALmod.fit, aes(x=Stigma_length, y=fit.CHALmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsCHAL, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CHAL%20fitted%20values%20plot-1.png)<!-- -->

``` r
ggplot() + 
  geom_point(data = CHAL, aes(x=Stigma_length, y=Log_flw_pollen)) + 
  geom_line(data= emmsCHAL, aes(xvar, yvar), linewidth=0.7) +
  # geom_smooth(method="lm", colour="blue", data = CHAL, aes(x=Stigma_length, y=Log_flw_pollen)) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
                     limits=c(0,NA)) +
  scale_x_continuous(name = "Stigma length (cm)") + 

    # scale_x_continuous(name = "Stigma length (cm)",
    #                  breaks = log(c(0.5, 1, 2, 5)),
    #                  labels = c(0.5, 1, 2, 5),
    #                  limits=c(log(0.5),log(5))) +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/CHAL%20fitted%20slope%20raw%20values%20plot-1.png)<!-- -->

### Plantago lanceolata (n=30)

``` r
PLLA <- stigl %>% 
  filter(Species=="Plantago lanceolata") %>% 
  as.data.frame()
# summary(PLLA) 
```

``` r
# fit full model with all terms
PLLAmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist + (1|Plant), data=PLLA)
summary(PLLAmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist +  
    ##     (1 | Plant)
    ##    Data: PLLA
    ## 
    ## REML criterion at convergence: 911.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.0580 -0.5262  0.1600  0.6804  1.9119 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.2758   0.5252  
    ##  Residual             1.3783   1.1740  
    ## Number of obs: 272, groups:  Plant, 30
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -1.805e+03  2.586e+02  2.545e+02  -6.978 2.58e-11 ***
    ## Infl_max      -5.891e-03  1.718e-02  2.619e+01  -0.343    0.734    
    ## Stigma_length  5.966e-01  1.321e-01  2.653e+02   4.516 9.49e-06 ***
    ## Date           1.437e-01  2.056e-02  2.545e+02   6.989 2.42e-11 ***
    ## Avg_dist       3.693e-03  3.056e-03  2.701e+01   1.209    0.237    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max    -0.020                     
    ## Stigm_lngth -0.050 -0.187              
    ## Date        -1.000  0.017  0.050       
    ## Avg_dist     0.012 -0.106 -0.113 -0.012

``` r
# anova(PLLAmod)
kable(anova(PLLAmod))
```

|               |     Sum Sq |    Mean Sq | NumDF |     DenDF |    F value |   Pr(\>F) |
|:--------------|-----------:|-----------:|------:|----------:|-----------:|----------:|
| Infl_max      |  0.1620357 |  0.1620357 |     1 |  26.19283 |  0.1175597 | 0.7344329 |
| Stigma_length | 28.1084221 | 28.1084221 |     1 | 265.34150 | 20.3931480 | 0.0000095 |
| Date          | 67.3198183 | 67.3198183 |     1 | 254.51238 | 48.8416964 | 0.0000000 |
| Avg_dist      |  2.0138966 |  2.0138966 |     1 |  27.00592 |  1.4611169 | 0.2372299 |

Significant effect of date and stigma length

``` r
PLLAmod.fit <- PLLA
fit.PLLAmod <- fitted(PLLAmod) 
PLLAmod.fit <- mutate(PLLAmod.fit, fit.PLLAmod)
emmsPLLA_stiglength <- emmip(PLLAmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
emmsPLLA_date <- emmip(PLLAmod, ~ Date, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = PLLAmod.fit, aes(x=Stigma_length, y=fit.PLLAmod)) + 
  # geom_smooth(data=PLLAmod.fit, aes(x=Stigma_length, y=fit.PLLAmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsPLLA_stiglength, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/PLLA%20stig%20length%20fitted%20values%20plot-1.png)<!-- -->

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = PLLAmod.fit, aes(x=Date, y=fit.PLLAmod)) + 
  # geom_smooth(data=PLLAmod.fit, aes(x=Stigma_length, y=fit.PLLAmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsPLLA_date, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Date", labels=as.character(unique(PLLAmod.fit$Date)), 
                     breaks=unique(PLLAmod.fit$Date)) +
  theme_cs() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](anal-stig-full_files/figure-gfm/PLLA%20date%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = PLLA, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsPLLA, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = PLLA, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Rumex acetosella (n=30)

``` r
RUAC <- stigl %>% 
  filter(Species=="Rumex acetosella") %>% 
  as.data.frame()
# summary(RUAC)
RUAC$Inv_avg_dist <-1/RUAC$Avg_dist
```

``` r
# fit full model with all terms
RUACmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=RUAC)
summary(RUACmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: RUAC
    ## 
    ## REML criterion at convergence: 471.8
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6040 -0.6825 -0.0184  0.5677  3.1786 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.06196  0.2489  
    ##  Residual             0.30708  0.5541  
    ## Number of obs: 261, groups:  Plant, 30
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -7.574e+02  9.427e+01  2.553e+02  -8.034 3.49e-14 ***
    ## Infl_max      -1.366e-02  8.582e-03  3.956e+01  -1.592 0.119363    
    ## Stigma_length  1.312e+00  3.940e-01  2.538e+02   3.331 0.000993 ***
    ## Date           6.026e-02  7.504e-03  2.553e+02   8.030 3.58e-14 ***
    ## Inv_avg_dist   1.473e+01  5.413e+00  3.531e+01   2.721 0.010043 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.354                     
    ## Stigm_lngth  0.022  0.040              
    ## Date        -1.000 -0.356 -0.025       
    ## Inv_avg_dst  0.038 -0.346  0.141 -0.037

``` r
anova(RUACmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
    ## Infl_max       0.7782  0.7782     1  39.563  2.5342 0.1193629    
    ## Stigma_length  3.4081  3.4081     1 253.808 11.0985 0.0009925 ***
    ## Date          19.8014 19.8014     1 255.329 64.4831 3.584e-14 ***
    ## Inv_avg_dist   2.2731  2.2731     1  35.307  7.4025 0.0100425 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
kable(anova(RUACmod))
```

|               |     Sum Sq |    Mean Sq | NumDF |     DenDF |   F value |   Pr(\>F) |
|:--------------|-----------:|-----------:|------:|----------:|----------:|----------:|
| Infl_max      |  0.7782011 |  0.7782011 |     1 |  39.56288 |  2.534202 | 0.1193629 |
| Stigma_length |  3.4081304 |  3.4081304 |     1 | 253.80785 | 11.098535 | 0.0009925 |
| Date          | 19.8014327 | 19.8014327 |     1 | 255.32934 | 64.483122 | 0.0000000 |
| Inv_avg_dist  |  2.2731427 |  2.2731427 |     1 |  35.30737 |  7.402461 | 0.0100425 |

Significant effect of date, stigma length, and avg dist.

``` r
RUACmod.fit <- RUAC
fit.RUACmod <- fitted(RUACmod) 
RUACmod.fit <- mutate(RUACmod.fit, fit.RUACmod)
emmsRUAC_stiglength <- emmip(RUACmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
emmsRUAC_date <- emmip(RUACmod, ~ Date, cov.reduce=range,
              type="response", plotit=F)
# emmsRUAC_avgdist <- emmip(RUACmod, ~ Avg_dist, cov.reduce=range,
#               type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = RUACmod.fit, aes(x=Stigma_length, y=fit.RUACmod)) + 
  # geom_smooth(data=RUACmod.fit, aes(x=Stigma_length, y=fit.RUACmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsRUAC_stiglength, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/RUAC%20stig%20length%20fitted%20values%20plot-1.png)<!-- -->

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = RUACmod.fit, aes(x=Date, y=fit.RUACmod)) + 
  # geom_smooth(data=RUACmod.fit, aes(x=Stigma_length, y=fit.RUACmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsRUAC_date, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Date", labels=as.character(unique(RUACmod.fit$Date)), 
                     breaks=unique(RUACmod.fit$Date)) +
  theme_cs() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](anal-stig-full_files/figure-gfm/RUAC%20date%20fitted%20values%20plot-1.png)<!-- -->

``` r
# converted avg_dist to inv_avg_dist -- not significant now
# # emmip_ggplot(emms) #plots predicted slope
# ggplot() + 
#   geom_point(data = RUACmod.fit, aes(x=Avg_dist, y=fit.RUACmod)) + 
#   # geom_smooth(data=RUACmod.fit, aes(x=Stigma_length, y=fit.RUACmod),
#   #             method="lm",
#   #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
#   geom_line(data= emmsRUAC_avgdist, aes(xvar, yvar), linewidth=0.7, colour="blue") +
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
#   scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
#   theme_cs()
```

``` r
# ggplot() + 
#   geom_point(data = RUAC, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsRUAC, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = RUAC, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Rumex crispus (n=25)

``` r
RUCR <- stigl %>% 
  filter(Species=="Rumex crispus") %>% 
  as.data.frame()
# summary(RUCR) 
RUCR$Inv_avg_dist <- 1/(RUCR$Avg_dist)
```

``` r
#Taking a look at the data
hist(RUCR$Flw_pollen) # response variable. right-skewed
```

![](anal-stig-full_files/figure-gfm/RUCR%20take%20a%20look-1.png)<!-- -->

``` r
hist(RUCR$Log_flw_pollen) # more symmetric 
```

![](anal-stig-full_files/figure-gfm/RUCR%20take%20a%20look-2.png)<!-- -->

``` r
hist(RUCR$Infl_max) # looks ok, kind of right skewed
```

![](anal-stig-full_files/figure-gfm/RUCR%20take%20a%20look-3.png)<!-- -->

``` r
hist(RUCR$Stigma_length) # looks ok  
```

![](anal-stig-full_files/figure-gfm/RUCR%20take%20a%20look-4.png)<!-- -->

``` r
hist(RUCR$Avg_dist) # not very symmetric at all!
```

![](anal-stig-full_files/figure-gfm/RUCR%20take%20a%20look-5.png)<!-- -->

``` r
hist(RUCR$Inv_avg_dist) # more symmetric 
```

![](anal-stig-full_files/figure-gfm/RUCR%20take%20a%20look-6.png)<!-- -->

``` r
table(RUCR$Date) # lots of dates
```

    ## 
    ## 2004-06-12 2004-06-13 2004-06-15 2004-06-16 2004-06-18 2004-06-19 2004-06-20 
    ##         14         14         26         44         45         24         12

``` r
# fit full model with all terms
RUCRmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=RUCR)
summary(RUCRmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: RUCR
    ## 
    ## REML criterion at convergence: 482.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.90505 -0.67486  0.04729  0.63199  3.04827 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.1429   0.3780  
    ##  Residual             0.7605   0.8721  
    ## Number of obs: 179, groups:  Plant, 25
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -825.36963  390.77865  173.98888  -2.112   0.0361 *  
    ## Infl_max         0.01575    0.01258   20.13124   1.251   0.2251    
    ## Stigma_length    1.63654    0.38708  173.42036   4.228 3.81e-05 ***
    ## Date             0.06546    0.03105  173.98912   2.108   0.0364 *  
    ## Inv_avg_dist   -12.25763   11.92750   17.09437  -1.028   0.3184    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max    -0.014                     
    ## Stigm_lngth  0.097  0.042              
    ## Date        -1.000  0.010 -0.098       
    ## Inv_avg_dst -0.005 -0.097  0.060  0.005

``` r
# anova(RUCRmod)
kable(anova(RUCRmod))
```

|               |     Sum Sq |    Mean Sq | NumDF |     DenDF |   F value |   Pr(\>F) |
|:--------------|-----------:|-----------:|------:|----------:|----------:|----------:|
| Infl_max      |  1.1910374 |  1.1910374 |     1 |  20.13124 |  1.566094 | 0.2251227 |
| Stigma_length | 15.0030904 | 15.0030904 |     1 | 173.31203 | 19.727545 | 0.0000159 |
| Date          |  1.5124599 |  1.5124599 |     1 |  25.71928 |  1.988732 | 0.1704536 |
| Inv_avg_dist  |  0.8059628 |  0.8059628 |     1 |  17.09404 |  1.059759 | 0.3176169 |

Significant effect of stigma length

``` r
RUCRmod.fit <- RUCR
fit.RUCRmod <- fitted(RUCRmod) 
RUCRmod.fit <- mutate(RUCRmod.fit, fit.RUCRmod)
emmsRUCR <- emmip(RUCRmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = RUCRmod.fit, aes(x=Stigma_length, y=fit.RUCRmod)) + 
  # geom_smooth(data=RUCRmod.fit, aes(x=Stigma_length, y=fit.RUCRmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsRUCR, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/RUCR%20stig%20length%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = RUCR, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsRUCR, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = RUCR, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Schizacne purpurascens (n=30)

``` r
SCPU <- stigl %>% 
  filter(Species=="Schizacne purpurascens") %>% 
  as.data.frame()
# summary(SCPU) 
SCPU$Inv_avg_dist <- 1/(SCPU$Avg_dist)
```

``` r
# fit full model with all terms
SCPUmod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=SCPU)
```

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

``` r
summary(SCPUmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: SCPU
    ## 
    ## REML criterion at convergence: 893.7
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.7219 -0.6224  0.1603  0.7781  1.9408 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.1542   0.3927  
    ##  Residual             2.1071   1.4516  
    ## Number of obs: 246, groups:  Plant, 30
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   3080.67087  630.87672  220.64976   4.883 2.00e-06 ***
    ## Infl_max         0.02251    0.02476   27.76606   0.909    0.371    
    ## Stigma_length    0.10680    0.30505  233.66821   0.350    0.727    
    ## Date            -0.24506    0.05020  220.62569  -4.881 2.02e-06 ***
    ## Inv_avg_dist    66.77495   48.30969   29.48168   1.382    0.177    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max    -0.055                     
    ## Stigm_lngth  0.153 -0.006              
    ## Date        -1.000  0.053 -0.154       
    ## Inv_avg_dst  0.015 -0.086 -0.072 -0.015
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling

``` r
anova(SCPUmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##               Sum Sq Mean Sq NumDF   DenDF F value Pr(>F)
    ## Infl_max      2.8327  2.8327     1  27.657  1.3444 0.2562
    ## Stigma_length 0.3669  0.3669     1 234.340  0.1741 0.6768
    ## Date          2.4599  2.4599     1  33.461  1.1674 0.2877
    ## Inv_avg_dist  3.6042  3.6042     1  29.467  1.7105 0.2010

No significant effect of anything - gives a warning that some predictors
are on very different scales? weird that only this one gives me that
error.

Looks like it deviates from normal qqline at upper extreme, no obvious
deviations from homoskedasticity in residuals vs fitted value plot

``` r
SCPUmod.fit <- SCPU
fit.SCPUmod <- fitted(SCPUmod) 
SCPUmod.fit <- mutate(SCPUmod.fit, fit.SCPUmod)
emmsSCPU <- emmip(SCPUmod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = SCPUmod.fit, aes(x=Stigma_length, y=fit.SCPUmod)) + 
  # geom_smooth(data=SCPUmod.fit, aes(x=Stigma_length, y=fit.SCPUmod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsSCPU, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/SCPU%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = SCPU, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsSCPU, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = SCPU, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Scirpus microcarpus (n=30)

``` r
SCMI <- stigl %>% 
  filter(Species=="Scirpus microcarpus") %>% 
  as.data.frame()
# summary(SCMI) 
SCMI$Inv_avg_dist <- 1/(SCMI$Avg_dist)
```

``` r
# fit full model with all terms
SCMImod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist + (1|Plant), data=SCMI)
summary(SCMImod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Avg_dist +  
    ##     (1 | Plant)
    ##    Data: SCMI
    ## 
    ## REML criterion at convergence: 651.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.95627 -0.72480 -0.07467  0.59454  2.67761 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Plant    (Intercept) 0.06093  0.2468  
    ##  Residual             0.56763  0.7534  
    ## Number of obs: 268, groups:  Plant, 30
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -3.580e+03  3.286e+02  2.429e+02 -10.895   <2e-16 ***
    ## Infl_max       8.011e-03  6.932e-03  2.125e+01   1.156    0.261    
    ## Stigma_length  6.959e-01  4.796e-01  2.582e+02   1.451    0.148    
    ## Date           2.844e-01  2.610e-02  2.429e+02  10.894   <2e-16 ***
    ## Avg_dist       5.560e-04  1.863e-03  2.134e+01   0.298    0.768    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.027                     
    ## Stigm_lngth  0.022  0.037              
    ## Date        -1.000 -0.029 -0.024       
    ## Avg_dist     0.031  0.448  0.008 -0.032

``` r
anova(SCMImod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF   DenDF F value  Pr(>F)  
    ## Infl_max      1.13551 1.13551     1  21.226  2.0004 0.17176  
    ## Stigma_length 1.67213 1.67213     1 258.129  2.9458 0.08730 .
    ## Date          2.74060 2.74060     1  53.303  4.8281 0.03237 *
    ## Avg_dist      0.26481 0.26481     1  21.325  0.4665 0.50194  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Significant effect of Date on pollen loads

``` r
# Check for collinearity between predictors
vif(mod=SCMImod, type="predictor") # variance inflation factors from package "car"
```

    ##      Infl_max Stigma_length          Date      Avg_dist 
    ##      1.252766      1.001972      1.001796      1.251393

``` r
#normality of residuals
hist(residuals(SCMImod))
```

![](anal-stig-full_files/figure-gfm/test%20SCMI%20model%20assumptions-1.png)<!-- -->

``` r
qqnorm(residuals(SCMImod))
qqline(residuals(SCMImod))
```

![](anal-stig-full_files/figure-gfm/test%20SCMI%20model%20assumptions-2.png)<!-- -->

``` r
# homogeneity
plot(SCMImod)
```

![](anal-stig-full_files/figure-gfm/test%20SCMI%20model%20assumptions-3.png)<!-- -->

``` r
#No obvious deviations from normality, some structure in the residuals vs fitted values plot.
```

``` r
SCMImod.fit <- SCMI
fit.SCMImod <- fitted(SCMImod) 
SCMImod.fit <- mutate(SCMImod.fit, fit.SCMImod)
emmsSCMI_date <- emmip(SCMImod, ~ Date, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = SCMImod.fit, aes(x=Date, y=fit.SCMImod)) + 
  # geom_smooth(data=SCMImod.fit, aes(x=Stigma_length, y=fit.SCMImod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsSCMI_date, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Date", labels=as.character(unique(SCMImod.fit$Date)), 
                     breaks=unique(SCMImod.fit$Date)) +
  theme_cs() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

![](anal-stig-full_files/figure-gfm/SCMI%20date%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = SCMI, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsSCMI, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = SCMI, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### Thalictrum dioicum (n=20)

``` r
THDI <- stigl %>% 
  filter(Species=="Thalictrum dioicum") %>% 
  as.data.frame()
# summary(THDI) 
THDI$Inv_avg_dist <- 1/(THDI$Avg_dist)
```

``` r
# fit full model with all terms
THDImod <- lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=THDI)
```

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

    ## Warning: Some predictor variables are on very different scales: consider
    ## rescaling

``` r
summary(THDImod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist +  
    ##     (1 | Plant)
    ##    Data: THDI
    ## 
    ## REML criterion at convergence: 235.8
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
    ##                 Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   -1.810e+03  3.964e+02  8.990e+01  -4.567 1.56e-05 ***
    ## Infl_max       2.177e-02  1.485e-02  1.261e+01   1.466   0.1672    
    ## Stigma_length  3.407e-01  1.355e-01  8.995e+01   2.514   0.0137 *  
    ## Date           1.443e-01  3.158e-02  8.990e+01   4.568 1.56e-05 ***
    ## Inv_avg_dist   6.283e+00  1.870e+01  1.161e+01   0.336   0.7429    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Infl_m Stgm_l Date  
    ## Infl_max     0.015                     
    ## Stigm_lngth  0.029  0.011              
    ## Date        -1.000 -0.016 -0.029       
    ## Inv_avg_dst  0.016 -0.025 -0.060 -0.016
    ## fit warnings:
    ## Some predictor variables are on very different scales: consider rescaling

``` r
anova(THDImod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## Infl_max      0.92605 0.92605     1 12.610  2.2939 0.154537   
    ## Stigma_length 2.86273 2.86273     1 89.843  7.0914 0.009179 **
    ## Date          0.19739 0.19739     1 25.458  0.4890 0.490735   
    ## Inv_avg_dist  0.05311 0.05311     1 11.606  0.1316 0.723329   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Significant effect of stigma length

No obvious deviations from normality or homoskedasticity

``` r
THDImod.fit <- THDI
fit.THDImod <- fitted(THDImod) 
THDImod.fit <- mutate(THDImod.fit, fit.THDImod)
emmsTHDI <- emmip(THDImod, ~ Stigma_length, cov.reduce=range,
              type="response", plotit=F)
```

``` r
# emmip_ggplot(emms) #plots predicted slope
ggplot() + 
  geom_point(data = THDImod.fit, aes(x=Stigma_length, y=fit.THDImod)) + 
  # geom_smooth(data=THDImod.fit, aes(x=Stigma_length, y=fit.THDImod),
  #             method="lm",
  #             formula = y ~ x, linetype=1, colour="black", se=F, linewidth=0.7) + 
  geom_line(data= emmsTHDI, aes(xvar, yvar), linewidth=0.7, colour="blue") +
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs()
```

![](anal-stig-full_files/figure-gfm/THDI%20fitted%20values%20plot-1.png)<!-- -->

``` r
# ggplot() + 
#   geom_point(data = THDI, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   geom_line(data= emmsTHDI, aes(xvar, yvar), linewidth=0.7) +
#   # geom_smooth(method="lm", colour="blue", data = THDI, aes(x=Stigma_length, y=Log_flw_pollen)) + 
#   scale_y_continuous(name="Stigmatic pollen load",
#                      breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#                      labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
#                      limits=c(0,NA)) +
#   scale_x_continuous(name = "Stigma length (cm)") + 
# 
#     # scale_x_continuous(name = "Stigma length (cm)",
#     #                  breaks = log(c(0.5, 1, 2, 5)),
#     #                  labels = c(0.5, 1, 2, 5),
#     #                  limits=c(log(0.5),log(5))) +
#   theme_cs()
```

### post

Run models all at once to get coefficient summaries easily:

``` r
#nest data by species to make it easier to handle
by_spp_stig <- stigl %>% group_by(Species, Sex_sys) %>% 
  nest() %>% 
  arrange(Sex_sys, Species)
```

``` r
stigl %>% 
  ggplot(aes(x=Stigma_length, y=Log_flw_pollen)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  facet_wrap(~Species, scales = "free", ncol=4) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs(font = "sans", fontsize=20) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 18, face="italic"))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-full_files/figure-gfm/all%20stig%20length%20plots-1.png)<!-- -->
