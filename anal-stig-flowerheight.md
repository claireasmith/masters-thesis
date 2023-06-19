Effects of flower height on stigmatic pollen load
================
Claire Smith
2023-06-19

``` r
# Load packages
library(tidyverse)
library(lme4) # for linear mixed effects models
```

    ## Warning: package 'lme4' was built under R version 4.3.1

``` r
library(lmerTest)
```

    ## Warning: package 'lmerTest' was built under R version 4.3.1

``` r
# Source files
source("clean-stig-01.R")
source("clean-stig-02.R")
source("theme_cs.R")
# Load data
stig <- read.csv("processed-data/stig-no-rep-spp.csv", stringsAsFactors = T)
```

``` r
# Take only entries with flower height data and pollen load data
stigfh <- stig %>% 
  filter(!is.na(Flower_height)&!is.na(Flw_pollen))
# head(stigfh)
# summary(stigfh)

# How many species have flower height data? 
# unique(stigfh$Species) # 7 species
# [1] Ambrosia artemisiifolia     Amaranthus retroflexus      Dichanthelium linearifolium Dichanthelium implicatum  
# [5] Phleum pratense             Setaria viridis             Rumex acetosella 
# View(stigfh)

# Flower is nested within plant - each plant has (potentially) multiple flowers measured within it -- need to 
# create a unique plant ID for each plant within a species
stigfh$Ind_ID=paste(stigfh$Date, stigfh$Site, stigfh$Plant, sep = "-")
#How many unique individuals per species? 
stigfh %>% group_by(Species) %>%
  distinct(Ind_ID) %>% 
  summarize(n=n())
```

    ## # A tibble: 7 × 2
    ##   Species                         n
    ##   <fct>                       <int>
    ## 1 Amaranthus retroflexus         30
    ## 2 Ambrosia artemisiifolia        30
    ## 3 Dichanthelium implicatum       40
    ## 4 Dichanthelium linearifolium    13
    ## 5 Phleum pratense                18
    ## 6 Rumex acetosella               50
    ## 7 Setaria viridis                45

``` r
# Species                         n
# <fct>                       <int>
# 1 Amaranthus retroflexus         30
# 2 Ambrosia artemisiifolia        30
# 3 Dichanthelium implicatum       40
# 4 Dichanthelium linearifolium    13
# 5 Phleum pratense                18
# 6 Rumex acetosella               50
# 7 Setaria viridis                45
# Will exclude Dichanthelium linearifolium from this because it only has 13 entries. 

# I'll go through each species and run a model predicting pollen capture by flower height
# Since many likely have log-normal distributions of pollen capture, I'll create a log flw pollen variable now to # save time -- but will still check each species individually 
stigfh$Log_flw_pollen <- log(stigfh$Flw_pollen+1)
```

### Amaranthus retroflexus (n=30)

``` r
AMRE <- stigfh %>% 
  filter(Species=="Amaranthus retroflexus")
# summary(AMRE) # 2 sites, FBF and FBF2

#Taking a look at the data
hist(AMRE$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-1.png)<!-- -->

``` r
hist(AMRE$Log_flw_pollen) # much more symmetric
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-2.png)<!-- -->

``` r
hist(AMRE$Flower_height) # looks a little right-skewed
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-3.png)<!-- -->

``` r
hist(log(AMRE$Flower_height+1)) # now a little left-skewed... I'll just use raw flower height
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-4.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Flw_height 
# looks like no big difference between sites in how much pollen they receive wrt height
plot(AMRE$Flower_height, AMRE$Log_flw_pollen, col=as.factor(AMRE$Site))
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-5.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Flw_height 
# Looks like no big differences between sites in flower height or pollen load
ggplot(data=AMRE, aes(x=Flower_height, y=Log_flw_pollen, color=Date, shape=Site)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Flower height") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Amaranthus%20retroflexus-1.png)<!-- -->

``` r
#model for height
# AMREmod <- lmer(Log_flw_pollen ~ poly(Flower_height, 2) + Site + (1|Ind_ID), data=AMRE)
# no sig quadratic term
AMREmod <- lmer(Log_flw_pollen ~ Flower_height + Site + (1|Ind_ID), data=AMRE)
summary(AMREmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Flower_height + Site + (1 | Ind_ID)
    ##    Data: AMRE
    ## 
    ## REML criterion at convergence: 268.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.34663 -0.68827  0.07175  0.57242  2.28491 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1031   0.3211  
    ##  Residual             1.0247   1.0123  
    ## Number of obs: 88, groups:  Ind_ID, 30
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)    3.015824   0.309523 51.515492   9.743 2.79e-13 ***
    ## Flower_height -0.007328   0.006886 67.432275  -1.064    0.291    
    ## SiteFBF2      -0.179445   0.254799 27.550579  -0.704    0.487    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Flwr_h
    ## Flower_hght -0.827       
    ## SiteFBF2    -0.165 -0.264

``` r
anova(AMREmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
    ## Flower_height 1.16047 1.16047     1 67.432  1.1325 0.2910
    ## Site          0.50824 0.50824     1 27.551  0.4960 0.4872

``` r
# No significant effect of flower height or site on pollen loads in Amaranthus retroflexus. 
```

### Ambrosia artemisiifolia (n=30)

``` r
AMAR <- stigfh %>% 
  filter(Species=="Ambrosia artemisiifolia")
summary(AMAR) # 2 sites, FBF and LP
```

    ##                     Species           Date         Site        Plant   
    ##  Ambrosia artemisiifolia:88   2021-09-14:44   FBF    :44   1      : 6  
    ##  Amaranthus retroflexus : 0   2021-09-20:44   LP     :44   10     : 6  
    ##  Bromus inermis         : 0   1189      : 0   BI     : 0   11     : 6  
    ##  Carex communis         : 0   1190      : 0   CT     : 0   12     : 6  
    ##  Carex hirtifolia       : 0   1191      : 0   FBF2   : 0   13     : 6  
    ##  Carex pedunculata      : 0   1192      : 0   KM     : 0   14     : 6  
    ##  (Other)                : 0   (Other)   : 0   (Other): 0   (Other):52  
    ##      Flower     Flw_pollen     Stigmas_per_flw Flower_height      Infl_max  
    ##  A      :30   Min.   :  0.00   Min.   :2       Min.   :10.00   Min.   : NA  
    ##  B      :29   1st Qu.: 11.00   1st Qu.:2       1st Qu.:30.75   1st Qu.: NA  
    ##  C      :29   Median : 23.00   Median :2       Median :49.00   Median : NA  
    ##  1      : 0   Mean   : 51.03   Mean   :2       Mean   :48.74   Mean   :NaN  
    ##  1-1    : 0   3rd Qu.: 70.25   3rd Qu.:2       3rd Qu.:67.00   3rd Qu.: NA  
    ##  1-10   : 0   Max.   :460.00   Max.   :2       Max.   :91.00   Max.   : NA  
    ##  (Other): 0                                                    NA's   :88   
    ##     Infl_min         D1               D2               D3        
    ##  Min.   : NA   Min.   :  8.00   Min.   : 16.00   Min.   : 17.00  
    ##  1st Qu.: NA   1st Qu.: 19.00   1st Qu.: 36.00   1st Qu.: 37.00  
    ##  Median : NA   Median : 37.00   Median : 67.00   Median : 59.00  
    ##  Mean   :NaN   Mean   : 55.62   Mean   : 67.36   Mean   : 68.11  
    ##  3rd Qu.: NA   3rd Qu.: 80.00   3rd Qu.: 95.00   3rd Qu.: 85.00  
    ##  Max.   : NA   Max.   :215.00   Max.   :165.00   Max.   :200.00  
    ##  NA's   :88                                                      
    ##        D4               D5         Stigma_length    source      Ind_ID         
    ##  Min.   : 15.00   Min.   : 19.00   Min.   : NA   CS2021:88   Length:88         
    ##  1st Qu.: 37.00   1st Qu.: 36.00   1st Qu.: NA   JF2001: 0   Class :character  
    ##  Median : 83.00   Median : 92.00   Median : NA   JF2004: 0   Mode  :character  
    ##  Mean   : 85.22   Mean   : 95.65   Mean   :NaN                                 
    ##  3rd Qu.:119.00   3rd Qu.:130.00   3rd Qu.: NA                                 
    ##  Max.   :240.00   Max.   :216.00   Max.   : NA                                 
    ##                                    NA's   :88                                  
    ##  Log_flw_pollen 
    ##  Min.   :0.000  
    ##  1st Qu.:2.485  
    ##  Median :3.178  
    ##  Mean   :3.233  
    ##  3rd Qu.:4.266  
    ##  Max.   :6.133  
    ## 

``` r
#Taking a look at the data
hist(AMAR$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Ambrosia%20artemisiifolia-1.png)<!-- -->

``` r
hist(AMAR$Log_flw_pollen) # much more symmetric
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Ambrosia%20artemisiifolia-2.png)<!-- -->

``` r
hist(AMAR$Flower_height) # looks ok
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Ambrosia%20artemisiifolia-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Flw_height 
# Seems like FBF received more pollen overall and was shorter. Collected a few days 
# earlier than LP - maybe plays a role? 
ggplot(data=AMAR, aes(x=Flower_height, y=Log_flw_pollen, color=Date, shape=Site)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Flower height") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Ambrosia%20artemisiifolia-1.png)<!-- -->

``` r
#model for height
# AMARmod <- lmer(Log_flw_pollen ~ poly(Flower_height,2) + Site + (1|Ind_ID), data=AMAR)
# no sig quadratic term
AMARmod <- lmer(Log_flw_pollen ~ Flower_height + Site + (1|Ind_ID), data=AMAR)
summary(AMARmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Flower_height + Site + (1 | Ind_ID)
    ##    Data: AMAR
    ## 
    ## REML criterion at convergence: 264.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.47153 -0.60514  0.05969  0.51063  2.38253 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1829   0.4276  
    ##  Residual             0.9217   0.9601  
    ## Number of obs: 88, groups:  Ind_ID, 30
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)    3.713080   0.311355 39.313351  11.926 1.24e-14 ***
    ## Flower_height  0.004939   0.006714 48.704063   0.736    0.465    
    ## SiteLP        -1.441993   0.298786 31.900158  -4.826 3.31e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Flwr_h
    ## Flower_hght -0.811       
    ## SiteLP       0.054 -0.506

``` r
anova(AMARmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Flower_height  0.4989  0.4989     1 48.704  0.5413    0.4654    
    ## Site          21.4687 21.4687     1 31.900 23.2919 3.314e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# anova(AMARmod2)
# anova(AMARmod, AMARmod2) # no sig difference between models
```

Site has a significant effect on pollen load, but not flower height, in
Ambrosia artemisiifolia. No significant interaction of site and flower
height. No difference in explanatory power between model with site x
flower height interaction vs without.

``` r
ggplot(data=AMAR, aes(x=Site, y=Log_flw_pollen, fill=Date)) + 
  geom_violin() + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/Ambrosia%20artemisiifolia%20site%20effects-1.png)<!-- -->

``` r
# FBF was taken on Sep 14 2021, LP was taken on Sep 20 2021 - may also be an effect of date on pollen capture.
# View(AMAR)
```

### Dichanthelium implicatum (n=40)

``` r
DIIM <- stigfh %>% 
  filter(Species=="Dichanthelium implicatum")
summary(DIIM) # 2 sites, KME and KMW, tho 3x the entries in KME than KMW
```

    ##                      Species           Date         Site        Plant   
    ##  Dichanthelium implicatum:58   2021-06-24:34   KME    :43   S1     : 7  
    ##  Amaranthus retroflexus  : 0   2021-06-21:15   KMW    :15   S2     : 5  
    ##  Ambrosia artemisiifolia : 0   2021-06-20: 9   BI     : 0   S3     : 5  
    ##  Bromus inermis          : 0   1189      : 0   CT     : 0   S10    : 3  
    ##  Carex communis          : 0   1190      : 0   FBF    : 0   S12    : 3  
    ##  Carex hirtifolia        : 0   1191      : 0   FBF2   : 0   S13    : 3  
    ##  (Other)                 : 0   (Other)   : 0   (Other): 0   (Other):32  
    ##      Flower     Flw_pollen    Stigmas_per_flw Flower_height      Infl_max    
    ##  A      :33   Min.   : 0.00   Min.   :2       Min.   :13.00   Min.   :13.00  
    ##  B      :16   1st Qu.: 2.25   1st Qu.:2       1st Qu.:22.00   1st Qu.:25.00  
    ##  m1     : 3   Median : 8.00   Median :2       Median :31.00   Median :33.00  
    ##  t1     : 2   Mean   :14.55   Mean   :2       Mean   :31.09   Mean   :33.32  
    ##  b1     : 1   3rd Qu.:23.75   3rd Qu.:2       3rd Qu.:40.88   3rd Qu.:43.88  
    ##  b2     : 1   Max.   :79.00   Max.   :2       Max.   :48.50   Max.   :59.00  
    ##  (Other): 2                                                                  
    ##     Infl_min           D1              D2              D3       
    ##  Min.   :12.00   Min.   : 6.00   Min.   : 7.00   Min.   : 8.00  
    ##  1st Qu.:21.50   1st Qu.: 9.00   1st Qu.:10.75   1st Qu.:13.00  
    ##  Median :30.50   Median :14.00   Median :16.00   Median :19.00  
    ##  Mean   :29.93   Mean   :14.79   Mean   :18.05   Mean   :20.29  
    ##  3rd Qu.:39.50   3rd Qu.:19.00   3rd Qu.:21.00   3rd Qu.:23.25  
    ##  Max.   :47.50   Max.   :41.00   Max.   :47.00   Max.   :57.00  
    ##                  NA's   :2       NA's   :2       NA's   :2      
    ##        D4              D5        Stigma_length    source      Ind_ID         
    ##  Min.   : 8.00   Min.   :10.00   Min.   : NA   CS2021:58   Length:58         
    ##  1st Qu.:11.50   1st Qu.:14.75   1st Qu.: NA   JF2001: 0   Class :character  
    ##  Median :19.00   Median :23.50   Median : NA   JF2004: 0   Mode  :character  
    ##  Mean   :21.45   Mean   :25.79   Mean   :NaN                                 
    ##  3rd Qu.:26.25   3rd Qu.:33.00   3rd Qu.: NA                                 
    ##  Max.   :54.00   Max.   :82.00   Max.   : NA                                 
    ##  NA's   :2       NA's   :2       NA's   :58                                  
    ##  Log_flw_pollen 
    ##  Min.   :0.000  
    ##  1st Qu.:1.171  
    ##  Median :2.197  
    ##  Mean   :2.153  
    ##  3rd Qu.:3.209  
    ##  Max.   :4.382  
    ## 

``` r
#Taking a look at the data
hist(DIIM$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Dichanthelium%20implicatum-1.png)<!-- -->

``` r
hist(DIIM$Log_flw_pollen) # much more symmetric
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Dichanthelium%20implicatum-2.png)<!-- -->

``` r
hist(DIIM$Flower_height) # looks ok
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Dichanthelium%20implicatum-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Flw_height 
# Much fewer points from KMW than KME, but they seemed to have received more pollen and were taller and less variable in height than KMW plants
ggplot(data=DIIM, aes(x=Flower_height, y=Log_flw_pollen, color=Date, shape=Site)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Flower height") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Dichanthelium%20implicatum-1.png)<!-- -->

``` r
#model for height and site
DIIMmod <- lmer(Log_flw_pollen ~ poly(Flower_height,2) + Site + (1|Ind_ID), data=DIIM)
summary(DIIMmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ poly(Flower_height, 2) + Site + (1 | Ind_ID)
    ##    Data: DIIM
    ## 
    ## REML criterion at convergence: 162.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.61828 -0.67683  0.06548  0.50260  1.99780 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.3145   0.5608  
    ##  Residual             0.7945   0.8913  
    ## Number of obs: 58, groups:  Ind_ID, 40
    ## 
    ## Fixed effects:
    ##                         Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)               1.8300     0.1864 29.9655   9.815 7.14e-11 ***
    ## poly(Flower_height, 2)1   0.2218     1.1831 38.2477   0.187  0.85227    
    ## poly(Flower_height, 2)2  -1.9877     1.1742 38.4250  -1.693  0.09858 .  
    ## SiteKMW                   1.0786     0.3531 42.3669   3.055  0.00389 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) p(F_,2)1 p(F_,2)2
    ## ply(Fl_,2)1  0.202                  
    ## ply(Fl_,2)2 -0.168 -0.049           
    ## SiteKMW     -0.591 -0.316    0.232

``` r
anova(DIIMmod) # quad term not sig, site sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                        Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## poly(Flower_height, 2) 2.2855  1.1427     2 38.291  1.4384 0.249844   
    ## Site                   7.4140  7.4140     1 42.367  9.3321 0.003885 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
DIIMmod2 <- lmer(Log_flw_pollen ~ poly(Flower_height,2) + Date + (1|Ind_ID), data=DIIM)
summary(DIIMmod2)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ poly(Flower_height, 2) + Date + (1 | Ind_ID)
    ##    Data: DIIM
    ## 
    ## REML criterion at convergence: 153.2
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.70738 -0.73107  0.01275  0.56363  1.81219 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.06562  0.2562  
    ##  Residual             0.84646  0.9200  
    ## Number of obs: 58, groups:  Ind_ID, 40
    ## 
    ## Fixed effects:
    ##                         Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)              3.02815    0.36715  8.43593   8.248 2.56e-05 ***
    ## poly(Flower_height, 2)1 -0.87828    1.06964 31.02446  -0.821  0.41786    
    ## poly(Flower_height, 2)2 -0.94778    1.07701 35.25512  -0.880  0.38482    
    ## Date2021-06-21          -0.01146    0.42798 13.70045  -0.027  0.97903    
    ## Date2021-06-24          -1.49094    0.43024 10.95354  -3.465  0.00532 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) p(F_,2)1 p(F_,2)2 D2021-06-21
    ## ply(Fl_,2)1 -0.219                              
    ## ply(Fl_,2)2  0.284 -0.161                       
    ## D2021-06-21 -0.796  0.021   -0.119              
    ## D2021-06-24 -0.908  0.313   -0.367    0.694

``` r
anova(DIIMmod2) # quad term not sig, date sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## poly(Flower_height, 2)  1.4609  0.7305     2 32.182  0.8629 0.4314492    
    ## Date                   19.4216  9.7108     2 15.432 11.4723 0.0008855 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
DIIMmod3 <- lmer(Log_flw_pollen ~ Flower_height + Site + (1|Ind_ID), data=DIIM)
summary(DIIMmod3)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Flower_height + Site + (1 | Ind_ID)
    ##    Data: DIIM
    ## 
    ## REML criterion at convergence: 176.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.52169 -0.72195  0.07591  0.55252  1.94650 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.3963   0.6295  
    ##  Residual             0.7654   0.8749  
    ## Number of obs: 58, groups:  Ind_ID, 40
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)    1.721209   0.487436 43.722884   3.531 0.000989 ***
    ## Flower_height  0.001494   0.015739 40.552157   0.095 0.924843    
    ## SiteKMW        1.226797   0.354335 45.304190   3.462 0.001180 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Flwr_h
    ## Flower_hght -0.923       
    ## SiteKMW      0.092 -0.321

``` r
anova(DIIMmod3) # flower height not sig, site sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##               Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)   
    ## Flower_height 0.0069  0.0069     1 40.552   0.009 0.92484   
    ## Site          9.1746  9.1746     1 45.304  11.987 0.00118 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
DIIMmod4 <- lmer(Log_flw_pollen ~ Flower_height + Date + (1|Ind_ID), data=DIIM)
summary(DIIMmod4)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Flower_height + Date + (1 | Ind_ID)
    ##    Data: DIIM
    ## 
    ## REML criterion at convergence: 164.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.69804 -0.73081 -0.02549  0.53307  1.98837 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.07127  0.2670  
    ##  Residual             0.83770  0.9153  
    ## Number of obs: 58, groups:  Ind_ID, 40
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)     3.53002    0.59819 15.33507   5.901 2.66e-05 ***
    ## Flower_height  -0.01321    0.01357 30.96896  -0.974  0.33781    
    ## Date2021-06-21 -0.05567    0.42609 13.77413  -0.131  0.89793    
    ## Date2021-06-24 -1.62952    0.40189  9.07052  -4.055  0.00282 ** 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Flwr_h D2021-06-21
    ## Flower_hght -0.814                   
    ## D2021-06-21 -0.476  0.002            
    ## D2021-06-24 -0.728  0.277  0.707

``` r
anova(DIIMmod4) # flower height not sig, date sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##               Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Flower_height  0.794   0.794     1 30.969  0.9479 0.3378063    
    ## Date          26.308  13.154     2 14.325 15.7026 0.0002451 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant quadratic or linear effect of flower height. Date and site both significant. 
```

``` r
ggplot(data=DIIM, aes(x=Date, y=Log_flw_pollen, fill=Site)) + 
  geom_violin() + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/Dichanthelium%20implicatum%20site%20and%20date%20effects-1.png)<!-- -->

### Phleum pratense (n=18)

``` r
PHPR <- stigfh %>% 
  filter(Species=="Phleum pratense")
summary(PHPR) # 1 site, 1 date
```

    ##                     Species           Date         Site        Plant   
    ##  Phleum pratense        :53   2021-07-05:53   KME    :53   S8     : 5  
    ##  Amaranthus retroflexus : 0   1189      : 0   BI     : 0   S18    : 4  
    ##  Ambrosia artemisiifolia: 0   1190      : 0   CT     : 0   S4     : 4  
    ##  Bromus inermis         : 0   1191      : 0   FBF    : 0   S6     : 4  
    ##  Carex communis         : 0   1192      : 0   FBF2   : 0   S7     : 4  
    ##  Carex hirtifolia       : 0   1194      : 0   KM     : 0   S10    : 3  
    ##  (Other)                : 0   (Other)   : 0   (Other): 0   (Other):29  
    ##      Flower     Flw_pollen     Stigmas_per_flw Flower_height   
    ##  B      :17   Min.   :  0.00   Min.   :2       Min.   : 22.00  
    ##  C      :17   1st Qu.: 14.00   1st Qu.:2       1st Qu.: 55.50  
    ##  A      :13   Median : 24.00   Median :2       Median : 63.00  
    ##  D      : 5   Mean   : 39.81   Mean   :2       Mean   : 66.06  
    ##  E      : 1   3rd Qu.: 50.00   3rd Qu.:2       3rd Qu.: 74.50  
    ##  1      : 0   Max.   :226.00   Max.   :2       Max.   :109.00  
    ##  (Other): 0                                                    
    ##     Infl_max         Infl_min            D1              D2        
    ##  Min.   : 24.00   Min.   : 21.50   Min.   : 10.0   Min.   : 12.00  
    ##  1st Qu.: 56.50   1st Qu.: 54.50   1st Qu.: 15.0   1st Qu.: 16.00  
    ##  Median : 64.50   Median : 62.00   Median : 24.0   Median : 54.00  
    ##  Mean   : 67.74   Mean   : 64.16   Mean   : 31.6   Mean   : 46.72  
    ##  3rd Qu.: 78.00   3rd Qu.: 72.50   3rd Qu.: 41.0   3rd Qu.: 62.00  
    ##  Max.   :111.00   Max.   :106.00   Max.   :111.0   Max.   :111.00  
    ##                                                                    
    ##        D3               D4              D5         Stigma_length    source  
    ##  Min.   : 13.00   Min.   : 26.0   Min.   : 12.00   Min.   : NA   CS2021:53  
    ##  1st Qu.: 17.00   1st Qu.: 55.0   1st Qu.: 62.00   1st Qu.: NA   JF2001: 0  
    ##  Median : 50.00   Median : 75.0   Median : 81.00   Median : NA   JF2004: 0  
    ##  Mean   : 50.68   Mean   : 74.3   Mean   : 78.06   Mean   :NaN              
    ##  3rd Qu.: 71.00   3rd Qu.: 88.0   3rd Qu.: 89.00   3rd Qu.: NA              
    ##  Max.   :117.00   Max.   :138.0   Max.   :125.00   Max.   : NA              
    ##                                                    NA's   :53               
    ##     Ind_ID          Log_flw_pollen 
    ##  Length:53          Min.   :0.000  
    ##  Class :character   1st Qu.:2.708  
    ##  Mode  :character   Median :3.219  
    ##                     Mean   :3.240  
    ##                     3rd Qu.:3.932  
    ##                     Max.   :5.425  
    ## 

``` r
#Taking a look at the data
hist(PHPR$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Phleum%20pratense-1.png)<!-- -->

``` r
hist(PHPR$Log_flw_pollen) # much more symmetric
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Phleum%20pratense-2.png)<!-- -->

``` r
hist(PHPR$Flower_height) # looks ok
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Phleum%20pratense-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Flw_height 
# Much fewer points from KMW than KME, but they seemed to have received more pollen and were taller and less variable in height than KMW plants
ggplot(data=PHPR, aes(x=Flower_height, y=Log_flw_pollen)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Flower height") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Phleum%20pratense-1.png)<!-- -->

``` r
#model for height and site
PHPRmod <- lmer(Log_flw_pollen ~ poly(Flower_height, 2) + (1|Ind_ID), data=PHPR)
# quad term not sig
summary(PHPRmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ poly(Flower_height, 2) + (1 | Ind_ID)
    ##    Data: PHPR
    ## 
    ## REML criterion at convergence: 147.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.76379 -0.59924  0.01586  0.41845  1.83401 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.3706   0.6088  
    ##  Residual             0.7955   0.8919  
    ## Number of obs: 53, groups:  Ind_ID, 18
    ## 
    ## Fixed effects:
    ##                         Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)               3.1963     0.1912 13.6720  16.719 1.71e-10 ***
    ## poly(Flower_height, 2)1  -0.5437     1.3200 16.4450  -0.412    0.686    
    ## poly(Flower_height, 2)2  -0.7419     1.2635 18.5310  -0.587    0.564    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) p(F_,2)1
    ## ply(Fl_,2)1 -0.021         
    ## ply(Fl_,2)2 -0.056  0.001

``` r
anova(PHPRmod) # no sig effect of flower height
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
    ## poly(Flower_height, 2) 0.40894 0.20447     2 17.418   0.257 0.7762

### Rumex acetosella (n=50)

``` r
RUAC <- stigfh %>% 
  filter(Species=="Rumex acetosella")
summary(RUAC) # 3 sites, 4 dates from June 4 to 9 2021
```

    ##                     Species            Date         Site        Plant   
    ##  Rumex acetosella       :110   2021-06-09:48   KMW    :78   F1     : 8  
    ##  Amaranthus retroflexus :  0   2021-06-07:30   WCK    :18   F14    : 7  
    ##  Ambrosia artemisiifolia:  0   2021-06-06:18   KM     :14   F15    : 7  
    ##  Bromus inermis         :  0   2021-06-04:14   BI     : 0   F5     : 7  
    ##  Carex communis         :  0   1189      : 0   CT     : 0   F4     : 6  
    ##  Carex hirtifolia       :  0   1190      : 0   FBF    : 0   F10    : 5  
    ##  (Other)                :  0   (Other)   : 0   (Other): 0   (Other):70  
    ##      Flower     Flw_pollen     Stigmas_per_flw Flower_height      Infl_max  
    ##  1      :50   Min.   : 0.000   Min.   :3       Min.   : 6.70   Min.   : NA  
    ##  2      :50   1st Qu.: 0.000   1st Qu.:3       1st Qu.:13.90   1st Qu.: NA  
    ##  3      : 9   Median : 1.000   Median :3       Median :19.30   Median : NA  
    ##  4      : 1   Mean   : 2.273   Mean   :3       Mean   :18.61   Mean   :NaN  
    ##  1-1    : 0   3rd Qu.: 4.000   3rd Qu.:3       3rd Qu.:22.95   3rd Qu.: NA  
    ##  1-10   : 0   Max.   :11.000   Max.   :3       Max.   :35.00   Max.   : NA  
    ##  (Other): 0                                                    NA's   :110  
    ##     Infl_min         D1               D2               D3        
    ##  Min.   : NA   Min.   :  5.00   Min.   :  4.00   Min.   :  6.50  
    ##  1st Qu.: NA   1st Qu.: 22.25   1st Qu.: 29.00   1st Qu.: 28.00  
    ##  Median : NA   Median : 45.00   Median : 56.00   Median : 52.50  
    ##  Mean   :NaN   Mean   : 56.74   Mean   : 67.49   Mean   : 70.63  
    ##  3rd Qu.: NA   3rd Qu.: 71.00   3rd Qu.: 88.50   3rd Qu.: 87.25  
    ##  Max.   : NA   Max.   :231.00   Max.   :230.00   Max.   :237.00  
    ##  NA's   :110                                                     
    ##        D4              D5         Stigma_length    source       Ind_ID         
    ##  Min.   :  9.0   Min.   :  8.00   Min.   : NA   CS2021:110   Length:110        
    ##  1st Qu.: 36.0   1st Qu.: 39.00   1st Qu.: NA   JF2001:  0   Class :character  
    ##  Median : 60.0   Median : 61.00   Median : NA   JF2004:  0   Mode  :character  
    ##  Mean   : 77.4   Mean   : 77.86   Mean   :NaN                                  
    ##  3rd Qu.:110.0   3rd Qu.:106.00   3rd Qu.: NA                                  
    ##  Max.   :236.0   Max.   :269.00   Max.   : NA                                  
    ##                                   NA's   :110                                  
    ##  Log_flw_pollen  
    ##  Min.   :0.0000  
    ##  1st Qu.:0.0000  
    ##  Median :0.6931  
    ##  Mean   :0.9272  
    ##  3rd Qu.:1.6094  
    ##  Max.   :2.4849  
    ## 

``` r
View(RUAC)

#Taking a look at the data
hist(RUAC$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Rumex%20acetosella-1.png)<!-- -->

``` r
hist(RUAC$Log_flw_pollen) # better
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Rumex%20acetosella-2.png)<!-- -->

``` r
hist(RUAC$Flower_height) # looks ok
```

![](anal-stig-flowerheight_files/figure-gfm/pre-fit%20Rumex%20acetosella-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Flw_height 
# overall doesn't seem to be much of a trend of flower height and pollen load
ggplot(data=RUAC, aes(x=Flower_height, y=Log_flw_pollen, color=Date, shape=Site)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Flower height") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Rumex%20acetosella-1.png)<!-- -->

``` r
#model for height and site
RUACmod <- lmer(Log_flw_pollen ~ poly(Flower_height,2) + Site + (1|Ind_ID), data=RUAC)
summary(RUACmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ poly(Flower_height, 2) + Site + (1 | Ind_ID)
    ##    Data: RUAC
    ## 
    ## REML criterion at convergence: 231.4
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.70233 -0.72381  0.00467  0.66810  1.78050 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.2193   0.4683  
    ##  Residual             0.3309   0.5752  
    ## Number of obs: 110, groups:  Ind_ID, 50
    ## 
    ## Fixed effects:
    ##                          Estimate Std. Error        df t value Pr(>|t|)   
    ## (Intercept)              0.953569   0.271374 36.431212   3.514   0.0012 **
    ## poly(Flower_height, 2)1  0.048982   0.981627 45.141445   0.050   0.9604   
    ## poly(Flower_height, 2)2 -0.629233   0.892121 45.207063  -0.705   0.4842   
    ## SiteKMW                  0.009856   0.296330 37.807076   0.033   0.9736   
    ## SiteWCK                 -0.222550   0.352103 34.881507  -0.632   0.5315   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) p(F_,2)1 p(F_,2)2 SitKMW
    ## ply(Fl_,2)1  0.261                         
    ## ply(Fl_,2)2  0.019  0.015                  
    ## SiteKMW     -0.938 -0.325   -0.034         
    ## SiteWCK     -0.729 -0.045    0.037    0.671

``` r
anova(RUACmod) # quad term not sig, site not sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
    ## poly(Flower_height, 2) 0.16582 0.082909     2 45.155  0.2506 0.7794
    ## Site                   0.25820 0.129098     2 36.535  0.3902 0.6797

``` r
RUACmod2 <- lmer(Log_flw_pollen ~ poly(Flower_height,2) + Date + (1|Ind_ID), data=RUAC)
summary(RUACmod2)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ poly(Flower_height, 2) + Date + (1 | Ind_ID)
    ##    Data: RUAC
    ## 
    ## REML criterion at convergence: 231.7
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.75828 -0.68082 -0.00278  0.62374  1.72532 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.2188   0.4678  
    ##  Residual             0.3308   0.5752  
    ## Number of obs: 110, groups:  Ind_ID, 50
    ## 
    ## Fixed effects:
    ##                         Estimate Std. Error       df t value Pr(>|t|)   
    ## (Intercept)              0.94482    0.27130 35.66907   3.483  0.00133 **
    ## poly(Flower_height, 2)1 -0.07474    0.98817 44.16360  -0.076  0.94005   
    ## poly(Flower_height, 2)2 -0.59710    0.89203 44.24594  -0.669  0.50674   
    ## Date2021-06-06          -0.22006    0.35183 34.14844  -0.625  0.53583   
    ## Date2021-06-07          -0.10956    0.31772 38.05055  -0.345  0.73211   
    ## Date2021-06-09           0.10363    0.30962 37.62096   0.335  0.73972   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) p(F_,2)1 p(F_,2)2 D2021-06-06 D2021-06-07
    ## ply(Fl_,2)1  0.263                                          
    ## ply(Fl_,2)2  0.018  0.011                                   
    ## D2021-06-06 -0.729 -0.046    0.037                          
    ## D2021-06-07 -0.863 -0.257   -0.045    0.623                 
    ## D2021-06-09 -0.906 -0.344   -0.023    0.644       0.785

``` r
anova(RUACmod2) # quad term not sig, date not sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                         Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
    ## poly(Flower_height, 2) 0.14976 0.074881     2 44.191  0.2264 0.7984
    ## Date                   0.61402 0.204673     3 38.492  0.6187 0.6071

``` r
RUACmod3 <- lmer(Log_flw_pollen ~ Flower_height + Site + (1|Ind_ID), data=RUAC)
summary(RUACmod3)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Flower_height + Site + (1 | Ind_ID)
    ##    Data: RUAC
    ## 
    ## REML criterion at convergence: 241.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6853 -0.7212 -0.0255  0.6573  1.8020 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.2155   0.4642  
    ##  Residual             0.3307   0.5751  
    ## Number of obs: 110, groups:  Ind_ID, 50
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error         df t value Pr(>|t|)   
    ## (Intercept)    0.9400336  0.3360451 40.2613626   2.797  0.00787 **
    ## Flower_height  0.0009336  0.0151981 46.1636830   0.061  0.95128   
    ## SiteKMW        0.0024426  0.2944877 38.5806852   0.008  0.99342   
    ## SiteWCK       -0.2136017  0.3498111 35.5859558  -0.611  0.54533   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Flwr_h SitKMW
    ## Flower_hght -0.632              
    ## SiteKMW     -0.480 -0.325       
    ## SiteWCK     -0.548 -0.046  0.673

``` r
anova(RUACmod3) # flower height not sig, site not sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                 Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
    ## Flower_height 0.001248 0.001248     1 46.164  0.0038 0.9513
    ## Site          0.229746 0.114873     2 37.245  0.3473 0.7088

``` r
RUACmod4 <- lmer(Log_flw_pollen ~ Flower_height + Date + (1|Ind_ID), data=RUAC)
summary(RUACmod4)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Flower_height + Date + (1 | Ind_ID)
    ##    Data: RUAC
    ## 
    ## REML criterion at convergence: 242.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.74410 -0.64246 -0.03057  0.63777  1.74456 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.2145   0.4631  
    ##  Residual             0.3307   0.5751  
    ## Number of obs: 110, groups:  Ind_ID, 50
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error        df t value Pr(>|t|)   
    ## (Intercept)     0.967696   0.336562 39.456855   2.875  0.00648 **
    ## Flower_height  -0.001043   0.015290 45.202486  -0.068  0.94594   
    ## Date2021-06-06 -0.211555   0.349287 34.848122  -0.606  0.54865   
    ## Date2021-06-07 -0.119306   0.315417 38.855060  -0.378  0.70731   
    ## Date2021-06-09  0.098659   0.307599 38.427201   0.321  0.75015   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Flwr_h D2021-06-06 D2021-06-07
    ## Flower_hght -0.635                               
    ## D2021-06-06 -0.546 -0.046                        
    ## D2021-06-07 -0.474 -0.257  0.626                 
    ## D2021-06-09 -0.435 -0.344  0.645       0.785

``` r
anova(RUACmod4) # flower height not sig, date not sig
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##                Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
    ## Flower_height 0.00154 0.001538     1 45.202  0.0046 0.9459
    ## Date          0.60657 0.202190     3 39.308  0.6114 0.6116

``` r
# No significant quadratic or linear effect of flower height. No effect of date or site.
```

``` r
ggplot(data=RUAC, aes(x=Date, y=Log_flw_pollen, fill=Site)) + 
  geom_violin() + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-flowerheight_files/figure-gfm/Rumex%20acetosella%20site%20and%20date%20effects-1.png)<!-- -->
