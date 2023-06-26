Effects of local density on pollen capture
================
Claire Smith
2023-06-20

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
# Take only entries with density data and pollen load data
stigd <- stig %>% 
  filter(!is.na(D1)&!is.na(Flw_pollen))
# head(stigd)
# summary(stigd)

# How many species have max infl height data? 
unique(stigd$Species) # 22 species!
```

    ##  [1] Ambrosia artemisiifolia     Amaranthus retroflexus     
    ##  [3] Dichanthelium linearifolium Dichanthelium implicatum   
    ##  [5] Plantago lanceolata         Phleum pratense            
    ##  [7] Setaria viridis             Rumex acetosella           
    ##  [9] Thalictrum dioicum          Plantago major             
    ## [11] Elymus repens               Bromus inermis             
    ## [13] Elymus innovatus            Festuca campestris         
    ## [15] Chenopodium album           Carex communis             
    ## [17] Carex hirtifolia            Carex pedunculata          
    ## [19] Carex plantaginea           Carex stipata              
    ## [21] Rumex crispus               Scirpus microcarpus        
    ## [23] Schizacne purpurascens     
    ## 23 Levels: Amaranthus retroflexus Ambrosia artemisiifolia ... Thalictrum dioicum

``` r
# View(stigd)

# Flower is nested within plant - each plant has (potentially) multiple flowers measured within it -- need to 
# create a unique plant ID for each plant within a species
stigd$Ind_ID=paste(stigd$Date, stigd$Site, stigd$Plant, sep = "-")
#How many unique individuals per species? 
print(n=22, stigd %>% group_by(Species) %>%
  distinct(Ind_ID) %>% 
  summarize(n=n()))
```

    ## # A tibble: 23 × 2
    ##    Species                         n
    ##    <fct>                       <int>
    ##  1 Amaranthus retroflexus         30
    ##  2 Ambrosia artemisiifolia        30
    ##  3 Bromus inermis                151
    ##  4 Carex communis                 81
    ##  5 Carex hirtifolia               89
    ##  6 Carex pedunculata              31
    ##  7 Carex plantaginea              38
    ##  8 Carex stipata                 142
    ##  9 Chenopodium album              95
    ## 10 Dichanthelium implicatum       39
    ## 11 Dichanthelium linearifolium    13
    ## 12 Elymus innovatus              163
    ## 13 Elymus repens                 100
    ## 14 Festuca campestris            319
    ## 15 Phleum pratense                18
    ## 16 Plantago lanceolata            42
    ## 17 Plantago major                 15
    ## 18 Rumex acetosella               50
    ## 19 Rumex crispus                 108
    ## 20 Schizacne purpurascens        154
    ## 21 Scirpus microcarpus           137
    ## 22 Setaria viridis                45
    ## # ℹ 1 more row

``` r
#    Species                      n
#    <fct>                    <int>
#  1 Amaranthus retroflexus         30
#  2 Ambrosia artemisiifolia        30
#  3 Bromus inermis                151
#  4 Carex communis                 81
#  5 Carex hirtifolia               89
#  6 Carex pedunculata              31
#  7 Carex plantaginea              38
#  8 Carex stipata                 142
#  9 Chenopodium album              95
# 10 Dichanthelium implicatum       39
# 11 Dichanthelium linearifolium    13* not enough samples
# 12 Elymus innovatus              163
# 13 Elymus repens                 100
# 14 Festuca campestris            319
# 15 Phleum pratense                18
# 16 Plantago lanceolata            27
# 17 Rumex acetosella               50
# 18 Rumex crispus                 108
# 19 Schizacne purpurascens        154
# 20 Scirpus microcarpus           137
# 21 Setaria viridis                45
# 22 Thalictrum dioicum             59

# I'll go through each species and run a model predicting pollen capture by density
# Since many likely have log-normal distributions of pollen capture, I'll create a log flw pollen variable now to # save time -- but will still check each species individually 
stigd$Log_flw_pollen <- log(stigd$Flw_pollen+1)

# I'll also create an estimate of the average distance to a plant's 5 nearest (staminate) neighbours
stigd <- stigd %>% 
  mutate(Avg_dist = (D1 + D2 + D3 + D4 + D5)/5)
```

### Amaranthus retroflexus (n=30)

``` r
AMRE <- stigd %>% 
  filter(Species=="Amaranthus retroflexus")
# summary(AMRE) 

#Taking a look at the data
hist(AMRE$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-1.png)<!-- -->

``` r
hist(AMRE$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-2.png)<!-- -->

``` r
hist(AMRE$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Amaranthus%20retroflexus-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist 
ggplot(data=AMRE, aes(x=Avg_dist, y=Log_flw_pollen, color=Date, shape=Site)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Amaranthus%20retroflexus-1.png)<!-- -->

``` r
AMREmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=AMRE)
summary(AMREmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: AMRE
    ## 
    ## REML criterion at convergence: 270
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.08476 -0.71881  0.00021  0.65501  2.49921 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.07108  0.2666  
    ##  Residual             1.05627  1.0278  
    ## Number of obs: 88, groups:  Ind_ID, 30
    ## 
    ## Fixed effects:
    ##                 Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)     2.846238   0.469645 25.914545   6.060 2.14e-06 ***
    ## Avg_dist       -0.002645   0.003286 25.713281  -0.805    0.428    
    ## Date2021-08-16  0.015528   0.377007 26.163479   0.041    0.967    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Avg_ds
    ## Avg_dist    -0.932       
    ## D2021-08-16 -0.882  0.771

``` r
anova(AMREmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
    ## Avg_dist 0.68439 0.68439     1 25.713  0.6479 0.4282
    ## Date     0.00179 0.00179     1 26.163  0.0017 0.9675

``` r
# No significant effect of density or date on pollen loads in Amaranthus retroflexus.  
```

### Ambrosia artemisiifolia (n=30)

``` r
AMAR <- stigd %>% 
  filter(Species=="Ambrosia artemisiifolia")
# summary(AMAR) 

#Taking a look at the data
hist(AMAR$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Ambrosia%20artemisiifolia-1.png)<!-- -->

``` r
hist(AMAR$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Ambrosia%20artemisiifolia-2.png)<!-- -->

``` r
hist(AMAR$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Ambrosia%20artemisiifolia-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=AMAR, aes(x=Avg_dist, y=Log_flw_pollen, color=Date, shape=Site)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Ambrosia%20artemisiifolia-1.png)<!-- -->

``` r
AMARmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=AMAR)
summary(AMARmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: AMAR
    ## 
    ## REML criterion at convergence: 266.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.46390 -0.60510  0.03388  0.53561  2.34903 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1827   0.4274  
    ##  Residual             0.9277   0.9632  
    ## Number of obs: 88, groups:  Ind_ID, 30
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)     3.8482769  0.2673803 26.8959330  14.393 3.73e-14 ***
    ## Avg_dist        0.0007097  0.0027425 28.0188484   0.259    0.798    
    ## Date2021-09-20 -1.3355238  0.2588200 27.5820199  -5.160 1.87e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Avg_ds
    ## Avg_dist    -0.731       
    ## D2021-09-20 -0.429 -0.071

``` r
anova(AMARmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Avg_dist  0.0621  0.0621     1 28.019   0.067    0.7977    
    ## Date     24.7002 24.7002     1 27.582  26.626 1.867e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Ambrosia artemisiifolia.  
```

### Bromus inermis (n=151)

``` r
BRIN <- stigd %>% 
  filter(Species=="Bromus inermis")
# summary(BRIN) 

#Taking a look at the data
hist(BRIN$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Bromus%20inermis-1.png)<!-- -->

``` r
hist(BRIN$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Bromus%20inermis-2.png)<!-- -->

``` r
hist(BRIN$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Bromus%20inermis-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=BRIN, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Bromus%20inermis-1.png)<!-- -->

``` r
BRINmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=BRIN)
summary(BRINmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: BRIN
    ## 
    ## REML criterion at convergence: 1273.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.54298 -0.59650  0.01525  0.64093  2.76813 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1918   0.438   
    ##  Residual             1.3265   1.152   
    ## Number of obs: 391, groups:  Ind_ID, 151
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   3.481604   0.462518 187.208829   7.527 2.12e-12 ***
    ## Avg_dist     -0.014792   0.005013 105.313005  -2.951  0.00391 ** 
    ## Date1190     -0.011279   0.480011 176.199107  -0.023  0.98128    
    ## Date1191      0.374577   0.444083 205.312502   0.843  0.39994    
    ## Date1192     -0.155549   0.468839 206.427809  -0.332  0.74039    
    ## Date1194      0.959340   0.483870 189.292816   1.983  0.04885 *  
    ## Date1195      0.008438   0.432014 203.500385   0.020  0.98444    
    ## Date1196      0.209681   0.485767 211.824705   0.432  0.66644    
    ## Date1200     -0.346025   0.579275 160.365178  -0.597  0.55112    
    ## Date1201      0.617384   0.522378 193.026744   1.182  0.23871    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt1190 Dt1191 Dt1192 Dt1194 Dt1195 Dt1196 Dt1200
    ## Avg_dist -0.455                                                        
    ## Date1190 -0.809  0.098                                                 
    ## Date1191 -0.879  0.115  0.807                                          
    ## Date1192 -0.812  0.064  0.760  0.822                                   
    ## Date1194 -0.833  0.165  0.747  0.809  0.759                            
    ## Date1195 -0.902  0.117  0.830  0.898  0.845  0.831                     
    ## Date1196 -0.786  0.068  0.734  0.795  0.750  0.733  0.817              
    ## Date1200 -0.661  0.060  0.616  0.667  0.629  0.615  0.685  0.607       
    ## Date1201 -0.742  0.087  0.685  0.742  0.699  0.686  0.762  0.675  0.566

``` r
anova(BRINmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## Avg_dist 11.548 11.5482     1 105.31  8.7059 0.003909 **
    ## Date     26.993  3.3741     8 135.83  2.5437 0.012985 * 
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Significant effect of density, sig effect of date on pollen loads in Bromus inermis.  
```

### Carex communis (n=81)

``` r
CACO <- stigd %>% 
  filter(Species=="Carex communis")
# summary(CACO) 

#Taking a look at the data
hist(CACO$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20communis-1.png)<!-- -->

``` r
hist(CACO$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20communis-2.png)<!-- -->

``` r
hist(CACO$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20communis-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=CACO, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

    ## Warning: Removed 9 rows containing missing values (`geom_point()`).

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Carex%20communis-1.png)<!-- -->

``` r
CACOmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=CACO)
summary(CACOmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: CACO
    ## 
    ## REML criterion at convergence: 533.9
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.92988 -0.61487 -0.01966  0.48766  3.11499 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.2889   0.5375  
    ##  Residual             0.6781   0.8235  
    ## Number of obs: 192, groups:  Ind_ID, 78
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)  1.336855   0.303197 69.708226   4.409 3.68e-05 ***
    ## Avg_dist    -0.005197   0.008606 67.917015  -0.604 0.547933    
    ## Date4127     0.527528   0.243529 70.218579   2.166 0.033694 *  
    ## Date4129     0.988191   0.244623 70.862971   4.040 0.000134 ***
    ## Date4130     1.444033   0.245584 68.071479   5.880 1.37e-07 ***
    ## Date4132     1.584891   0.416707 99.366516   3.803 0.000247 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4127 Dt4129 Dt4130
    ## Avg_dist -0.837                            
    ## Date4127 -0.396  0.028                     
    ## Date4129 -0.379  0.010  0.462              
    ## Date4130 -0.331 -0.045  0.458  0.457       
    ## Date4132 -0.197 -0.024  0.270  0.269  0.270

``` r
anova(CACOmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
    ## Avg_dist  0.2473  0.2473     1 67.917  0.3647   0.5479    
    ## Date     29.2781  7.3195     4 77.666 10.7943 5.21e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Carex communis.  
```

### Carex hirtifolia (n=89)

``` r
CAHI <- stigd %>% 
  filter(Species=="Carex hirtifolia")
# summary(CAHI) 

#Taking a look at the data
hist(CAHI$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20hirtifolia-1.png)<!-- -->

``` r
hist(CAHI$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20hirtifolia-2.png)<!-- -->

``` r
hist(CAHI$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20hirtifolia-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=CAHI, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Carex%20hirtifolia-1.png)<!-- -->

``` r
CAHImod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=CAHI)
summary(CAHImod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: CAHI
    ## 
    ## REML criterion at convergence: 701.2
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9824 -0.6259 -0.1317  0.5172  2.9889 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1284   0.3584  
    ##  Residual             1.5295   1.2367  
    ## Number of obs: 207, groups:  Ind_ID, 89
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)  1.806550   0.265415 88.991720   6.807 1.13e-09 ***
    ## Avg_dist    -0.003765   0.004350 73.071734  -0.865 0.389655    
    ## Date4133     0.833048   0.271674 91.034500   3.066 0.002853 ** 
    ## Date4134     1.031945   0.273095 91.369165   3.779 0.000281 ***
    ## Date4135     0.453104   0.282235 94.663311   1.605 0.111734    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4133 Dt4134
    ## Avg_dist -0.650                     
    ## Date4133 -0.512 -0.081              
    ## Date4134 -0.477 -0.130  0.559       
    ## Date4135 -0.457 -0.133  0.542  0.546

``` r
anova(CAHImod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## Avg_dist  1.1455  1.1455     1 73.072  0.7489 0.389655   
    ## Date     25.3745  8.4582     3 83.492  5.5300 0.001637 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Carex hirtifolia  
```

### Carex pedunculata (n=31)

``` r
CAPE <- stigd %>% 
  filter(Species=="Carex pedunculata")
# summary(CAPE) 

#Taking a look at the data
hist(CAPE$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20pedunculata-1.png)<!-- -->

``` r
hist(CAPE$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20pedunculata-2.png)<!-- -->

``` r
hist(CAPE$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20pedunculata-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=CAPE, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Carex%20pedunculata-1.png)<!-- -->

``` r
CAPEmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=CAPE)
summary(CAPEmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: CAPE
    ## 
    ## REML criterion at convergence: 206
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.31085 -0.47941 -0.00848  0.57243  2.07328 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.3816   0.6178  
    ##  Residual             0.4665   0.6830  
    ## Number of obs: 80, groups:  Ind_ID, 31
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)  1.02646    0.71027 30.12904   1.445    0.159    
    ## Avg_dist    -0.00880    0.01674 28.75892  -0.526    0.603    
    ## Date4118     1.76491    0.38615 29.47189   4.571 8.09e-05 ***
    ## Date4119     1.74576    0.32838 30.73087   5.316 8.89e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4118
    ## Avg_dist -0.928              
    ## Date4118 -0.316  0.064       
    ## Date4119 -0.239 -0.068  0.551

``` r
anova(CAPEmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Avg_dist  0.1289  0.1289     1 28.759  0.2763    0.6032    
    ## Date     14.9930  7.4965     2 28.605 16.0700 2.099e-05 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Carex pedunculata
```

### Carex plantaginea (n=38)

``` r
CAPL <- stigd %>% 
  filter(Species=="Carex plantaginea")
# summary(CAPL) 

#Taking a look at the data
hist(CAPL$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20plantaginea-1.png)<!-- -->

``` r
hist(CAPL$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20plantaginea-2.png)<!-- -->

``` r
hist(CAPL$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20plantaginea-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=CAPL, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

    ## Warning: Removed 32 rows containing missing values (`geom_point()`).

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Carex%20plantaginea-1.png)<!-- -->

``` r
CAPLmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=CAPL)
summary(CAPLmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: CAPL
    ## 
    ## REML criterion at convergence: 128.1
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.74587 -0.59732 -0.05723  0.63824  2.60979 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1303   0.3610  
    ##  Residual             0.2761   0.5255  
    ## Number of obs: 64, groups:  Ind_ID, 25
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error        df t value Pr(>|t|)   
    ## (Intercept)  0.901889   0.246609 17.996810   3.657   0.0018 **
    ## Avg_dist     0.006413   0.003649 15.281181   1.757   0.0989 . 
    ## Date4119    -0.579936   0.276647 19.202303  -2.096   0.0495 * 
    ## Date4128     0.827972   0.308462 16.174188   2.684   0.0162 * 
    ## Date4130     0.130127   0.333605 16.418743   0.390   0.7015   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4119 Dt4128
    ## Avg_dist -0.485                     
    ## Date4119 -0.582 -0.207              
    ## Date4128 -0.499 -0.233  0.594       
    ## Date4130 -0.426 -0.287  0.564  0.519

``` r
anova(CAPLmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
    ## Avg_dist 0.8527 0.85267     1 15.281  3.0883 0.098869 .  
    ## Date     7.8825 2.62750     3 15.847  9.5165 0.000783 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Carex plantaginea
```

### Carex stipata (n=142)

``` r
CAST <- stigd %>% 
  filter(Species=="Carex stipata")
# summary(CAST) 

#Taking a look at the data
hist(CAST$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20stipata-1.png)<!-- -->

``` r
hist(CAST$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20stipata-2.png)<!-- -->

``` r
hist(CAST$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Carex%20stipata-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=CAST, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Carex%20stipata-1.png)<!-- -->

``` r
CASTmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=CAST)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
# boundary fit is singluar
summary(CASTmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: CAST
    ## 
    ## REML criterion at convergence: 901.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6394 -0.6836 -0.1690  0.3835  3.8522 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.0000   0.0000  
    ##  Residual             0.9165   0.9574  
    ## Number of obs: 321, groups:  Ind_ID, 142
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)  6.631e-01  1.770e-01  3.150e+02   3.747 0.000213 ***
    ## Avg_dist    -2.300e-03  7.006e-04  3.150e+02  -3.282 0.001146 ** 
    ## Date4147     6.999e-01  1.703e-01  3.150e+02   4.109 5.07e-05 ***
    ## Date4148     1.079e+00  1.583e-01  3.150e+02   6.814 4.83e-11 ***
    ## Date4149     8.368e-01  1.722e-01  3.150e+02   4.858 1.87e-06 ***
    ## Date4150     1.016e+00  2.102e-01  3.150e+02   4.832 2.11e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4147 Dt4148 Dt4149
    ## Avg_dist -0.710                            
    ## Date4147 -0.510 -0.008                     
    ## Date4148 -0.513 -0.058  0.576              
    ## Date4149 -0.485 -0.034  0.530  0.572       
    ## Date4150 -0.420  0.003  0.434  0.467  0.429
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
anova(CASTmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF DenDF F value   Pr(>F)    
    ## Avg_dist  9.873  9.8734     1   315  10.773 0.001146 ** 
    ## Date     46.473 11.6183     4   315  12.676 1.41e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Significant effect of density, sig effect of date on pollen loads in Carex stipata
```

### Chenopodium album (n=95)

``` r
CHAL <- stigd %>% 
  filter(Species=="Chenopodium album")
# summary(CHAL) 

#Taking a look at the data
hist(CHAL$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Chenopodium%20album-1.png)<!-- -->

``` r
hist(CHAL$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Chenopodium%20album-2.png)<!-- -->

``` r
hist(CHAL$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Chenopodium%20album-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=CHAL, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Chenopodium%20album-1.png)<!-- -->

``` r
CHALmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=CHAL)
summary(CHALmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: CHAL
    ## 
    ## REML criterion at convergence: 377.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.67796 -0.72644 -0.03995  0.67742  2.71077 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1097   0.3312  
    ##  Residual             0.5805   0.7619  
    ## Number of obs: 149, groups:  Ind_ID, 95
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)   
    ## (Intercept)   0.987970   0.321245 108.984384   3.075  0.00266 **
    ## Avg_dist      0.001638   0.005048  77.292377   0.325  0.74643   
    ## Date4201     -0.279102   0.336824  94.056759  -0.829  0.40941   
    ## Date4202      0.085856   0.288641 131.592345   0.297  0.76659   
    ## Date4203      0.185285   0.270505 110.304415   0.685  0.49481   
    ## Date4204     -0.115501   0.270157 114.221457  -0.428  0.66980   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4201 Dt4202 Dt4203
    ## Avg_dist -0.689                            
    ## Date4201 -0.605  0.150                     
    ## Date4202 -0.698  0.165  0.583              
    ## Date4203 -0.748  0.180  0.622  0.724       
    ## Date4204 -0.739  0.166  0.621  0.723  0.772

``` r
anova(CHALmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
    ## Avg_dist 0.06112 0.06112     1 77.292  0.1053 0.7464
    ## Date     2.60343 0.65086     4 82.990  1.1213 0.3521

``` r
# No significant effect of density or date on pollen loads in Chenopodium album
```

### Dichanthelium implicatum (n=39)

``` r
DIIM <- stigd %>% 
  filter(Species=="Dichanthelium implicatum")
# summary(DIIM) 

#Taking a look at the data
hist(DIIM$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Dichanthelium%20implicatum-1.png)<!-- -->

``` r
hist(DIIM$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Dichanthelium%20implicatum-2.png)<!-- -->

``` r
hist(DIIM$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Dichanthelium%20implicatum-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=DIIM, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Dichanthelium%20implicatum-1.png)<!-- -->

``` r
DIIMmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=DIIM)
summary(DIIMmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: DIIM
    ## 
    ## REML criterion at convergence: 157.1
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.6562 -0.7211  0.0793  0.5567  1.8866 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1001   0.3163  
    ##  Residual             0.7908   0.8893  
    ## Number of obs: 56, groups:  Ind_ID, 39
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)     3.16921    0.46172 10.06508   6.864 4.24e-05 ***
    ## Avg_dist       -0.02335    0.01704 37.24925  -1.370   0.1787    
    ## Date2021-06-21  0.43540    0.53182 13.45077   0.819   0.4272    
    ## Date2021-06-24 -1.18706    0.46038  9.00958  -2.578   0.0297 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Avg_ds D2021-06-21
    ## Avg_dist    -0.451                   
    ## D2021-06-21 -0.495 -0.436            
    ## D2021-06-24 -0.683 -0.257  0.805

``` r
anova(DIIMmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Avg_dist  1.4853  1.4853     1 37.249  1.8782 0.1787454    
    ## Date     24.1065 12.0533     2 15.264 15.2418 0.0002301 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, significant effect of date on pollen loads in Dichanthelium implicatum
```

### Elymus innovatus (n=163)

``` r
ELIN <- stigd %>% 
  filter(Species=="Elymus innovatus")
# summary(ELIN) 

#Taking a look at the data
hist(ELIN$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Elymus%20innovatus-1.png)<!-- -->

``` r
hist(ELIN$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Elymus%20innovatus-2.png)<!-- -->

``` r
hist(ELIN$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Elymus%20innovatus-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=ELIN, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Elymus%20innovatus-1.png)<!-- -->

``` r
ELINmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=ELIN)
summary(ELINmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: ELIN
    ## 
    ## REML criterion at convergence: 1017.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.8421 -0.4824  0.0444  0.5438  3.4210 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1129   0.3359  
    ##  Residual             0.8947   0.9459  
    ## Number of obs: 354, groups:  Ind_ID, 163
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)  4.188e+00  1.553e-01  1.416e+02  26.973  < 2e-16 ***
    ## Avg_dist     7.134e-04  2.184e-03  1.551e+02   0.327   0.7443    
    ## Date2188    -5.606e-01  2.251e-01  1.769e+02  -2.490   0.0137 *  
    ## Date2189    -1.283e+00  3.202e-01  2.109e+02  -4.006 8.55e-05 ***
    ## Date2190     1.998e-02  1.423e-01  1.403e+02   0.140   0.8886    
    ## Date2191    -1.185e-01  1.821e-01  1.385e+02  -0.651   0.5163    
    ## Date2192    -5.777e-03  2.448e-01  1.600e+02  -0.024   0.9812    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt2188 Dt2189 Dt2190 Dt2191
    ## Avg_dist -0.762                                   
    ## Date2188 -0.262 -0.037                            
    ## Date2189 -0.207  0.005  0.140                     
    ## Date2190 -0.422 -0.047  0.318  0.222              
    ## Date2191 -0.230 -0.167  0.253  0.173  0.398       
    ## Date2192 -0.138 -0.169  0.190  0.128  0.298  0.255

``` r
anova(ELINmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Avg_dist  0.0955  0.0955     1 155.12  0.1067 0.7443297    
    ## Date     20.5633  4.1127     5 168.67  4.5964 0.0005881 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Elymus innovatus
```

### Elymus repens (n=100)

``` r
ELRE <- stigd %>% 
  filter(Species=="Elymus repens")
# summary(ELRE) 

#Taking a look at the data
hist(ELRE$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Elymus%20repens-1.png)<!-- -->

``` r
hist(ELRE$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Elymus%20repens-2.png)<!-- -->

``` r
hist(ELRE$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Elymus%20repens-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=ELRE, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Elymus%20repens-1.png)<!-- -->

``` r
ELREmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=ELRE)
summary(ELREmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: ELRE
    ## 
    ## REML criterion at convergence: 623.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.65523 -0.59491  0.05273  0.60278  2.37690 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.09201  0.3033  
    ##  Residual             0.50267  0.7090  
    ## Number of obs: 266, groups:  Ind_ID, 100
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)  3.705614   0.245524 88.364171  15.093  < 2e-16 ***
    ## Avg_dist     0.002942   0.007293 78.250812   0.403   0.6877    
    ## Date1213     1.160128   0.241702 84.516284   4.800 6.78e-06 ***
    ## Date1214     0.304357   0.232713 87.701808   1.308   0.1943    
    ## Date1215     0.175381   0.220204 76.502755   0.796   0.4282    
    ## Date1216     0.424070   0.237614 75.228659   1.785   0.0783 .  
    ## Date1217     0.402669   0.254592 95.825231   1.582   0.1170    
    ## Date1219    -0.456323   0.321836 99.629743  -1.418   0.1593    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt1213 Dt1214 Dt1215 Dt1216 Dt1217
    ## Avg_dist -0.613                                          
    ## Date1213 -0.647  0.020                                   
    ## Date1214 -0.638 -0.035  0.668                            
    ## Date1215 -0.691 -0.009  0.707  0.735                     
    ## Date1216 -0.701  0.091  0.657  0.677  0.718              
    ## Date1217 -0.599 -0.005  0.611  0.635  0.671  0.622       
    ## Date1219 -0.465 -0.019  0.483  0.503  0.531  0.490  0.459

``` r
anova(ELREmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Avg_dist  0.0818  0.0818     1 78.251  0.1628    0.6877    
    ## Date     24.0553  4.0092     6 81.756  7.9758 8.761e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Elymus repens
```

### Festuca campestris (n=319)

``` r
FECA <- stigd %>% 
  filter(Species=="Festuca campestris")
# summary(FECA) 

#Taking a look at the data
hist(FECA$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Festuca%20campestris-1.png)<!-- -->

``` r
hist(FECA$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Festuca%20campestris-2.png)<!-- -->

``` r
hist(FECA$Avg_dist) # looks ok
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Festuca%20campestris-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=FECA, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Festuca%20campestris-1.png)<!-- -->

``` r
FECAmod <- lmer(Log_flw_pollen ~ Avg_dist + Date + (1|Ind_ID), data=FECA)
summary(FECAmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: FECA
    ## 
    ## REML criterion at convergence: 2463.4
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.9792 -0.6945 -0.1117  0.6382  3.0296 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1066   0.3265  
    ##  Residual             1.4417   1.2007  
    ## Number of obs: 748, groups:  Ind_ID, 319
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   1.652674   0.322381 145.759070   5.126 9.24e-07 ***
    ## Avg_dist      0.000318   0.001697 287.865980   0.187   0.8515    
    ## Date2185      0.414322   0.370118 182.760402   1.119   0.2644    
    ## Date2187      0.374561   0.338997 147.492664   1.105   0.2710    
    ## Date2188      0.356978   0.334173 151.895141   1.068   0.2871    
    ## Date2189      0.690629   0.338371 161.076100   2.041   0.0429 *  
    ## Date2190     -0.143140   0.335033 154.165573  -0.427   0.6698    
    ## Date2191      0.781801   0.342298 150.230776   2.284   0.0238 *  
    ## Date2192      0.622097   0.393514 217.673550   1.581   0.1154    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt2185 Dt2187 Dt2188 Dt2189 Dt2190 Dt2191
    ## Avg_dist -0.188                                                 
    ## Date2185 -0.841  0.002                                          
    ## Date2187 -0.919  0.010  0.799                                   
    ## Date2188 -0.930 -0.003  0.811  0.885                            
    ## Date2189 -0.920  0.007  0.801  0.874  0.887                     
    ## Date2190 -0.928  0.000  0.809  0.883  0.896  0.884              
    ## Date2191 -0.909  0.000  0.791  0.864  0.877  0.866  0.874       
    ## Date2192 -0.794  0.019  0.688  0.752  0.762  0.753  0.761  0.744

``` r
anova(FECAmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)    
    ## Avg_dist  0.051  0.0506     1 287.87  0.0351   0.8515    
    ## Date     60.233  8.6047     7 256.31  5.9686 1.89e-06 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Festuca campestris
```

### Phleum pratense (n=18)

``` r
PHPR <- stigd %>% 
  filter(Species=="Phleum pratense")
# summary(PHPR) 

#Taking a look at the data
hist(PHPR$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Phleum%20pratense-1.png)<!-- -->

``` r
hist(PHPR$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Phleum%20pratense-2.png)<!-- -->

``` r
hist(PHPR$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Phleum%20pratense-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=PHPR, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Phleum%20pratense-1.png)<!-- -->

``` r
PHPRmod <- lmer(Log_flw_pollen ~ Avg_dist + (1|Ind_ID), data=PHPR)
summary(PHPRmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + (1 | Ind_ID)
    ##    Data: PHPR
    ## 
    ## REML criterion at convergence: 159.8
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.67847 -0.56108  0.00625  0.48222  1.89965 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.3323   0.5765  
    ##  Residual             0.7977   0.8932  
    ## Number of obs: 53, groups:  Ind_ID, 18
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)  2.901613   0.522885 15.161114   5.549 5.36e-05 ***
    ## Avg_dist     0.005094   0.008612 15.732448   0.591    0.563    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr)
    ## Avg_dist -0.935

``` r
anova(PHPRmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
    ## Avg_dist 0.27908 0.27908     1 15.732  0.3498 0.5626

``` r
# No significant effect of density on pollen loads in Phleum pratense
```

### Plantago lanceolata (n=27)

``` r
PLLA <- stigd %>% 
  filter(Species=="Plantago lanceolata")
# summary(PLLA) 

#Taking a look at the data
hist(PLLA$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Plantago%20lanceolata-1.png)<!-- -->

``` r
hist(PLLA$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Plantago%20lanceolata-2.png)<!-- -->

``` r
hist(PLLA$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Plantago%20lanceolata-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=PLLA, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Plantago%20lanceolata-1.png)<!-- -->

``` r
PLLAmod <- lmer(Log_flw_pollen ~ Avg_dist + Date +  (1|Ind_ID), data=PLLA)
summary(PLLAmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: PLLA
    ## 
    ## REML criterion at convergence: 364
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.9204 -0.6007  0.1097  0.6300  1.6369 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.4802   0.6930  
    ##  Residual             0.8251   0.9084  
    ## Number of obs: 120, groups:  Ind_ID, 42
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)     4.15640    0.30711 36.88897  13.534 7.03e-16 ***
    ## Avg_dist        0.01513    0.01579 37.05481   0.958   0.3444    
    ## Date2021-07-31 -1.26950    0.53488 37.56157  -2.373   0.0228 *  
    ## Date2021-08-18 -0.15805    0.32438 37.42704  -0.487   0.6289    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Avg_ds D2021-07
    ## Avg_dist    -0.674                
    ## D2021-07-31  0.208 -0.773         
    ## D2021-08-18 -0.590  0.109  0.212

``` r
anova(PLLAmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
    ## Avg_dist 0.7569 0.75694     1 37.055  0.9173 0.34438  
    ## Date     4.6483 2.32415     2 37.657  2.8167 0.07245 .
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density or date on pollen loads in Plantago lanceolata
```

### Rumex acetosella (n=50)

``` r
RUAC <- stigd %>% 
  filter(Species=="Rumex acetosella")
# summary(RUAC) 

#Taking a look at the data
hist(RUAC$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Rumex%20acetosella-1.png)<!-- -->

``` r
hist(RUAC$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Rumex%20acetosella-2.png)<!-- -->

``` r
hist(RUAC$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Rumex%20acetosella-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=RUAC, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Rumex%20acetosella-1.png)<!-- -->

``` r
RUACmod <- lmer(Log_flw_pollen ~ Avg_dist + Date +  (1|Ind_ID), data=RUAC)
summary(RUACmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: RUAC
    ## 
    ## REML criterion at convergence: 246.3
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.75214 -0.64294 -0.02885  0.63493  1.73567 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.2143   0.4629  
    ##  Residual             0.3308   0.5751  
    ## Number of obs: 110, groups:  Ind_ID, 50
    ## 
    ## Fixed effects:
    ##                  Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)     0.9575795  0.2644662 36.1704670   3.621 0.000893 ***
    ## Avg_dist       -0.0001696  0.0018549 45.4123610  -0.091 0.927531    
    ## Date2021-06-06 -0.2065281  0.3552324 35.1094348  -0.581 0.564693    
    ## Date2021-06-07 -0.1165427  0.3179923 38.9491939  -0.366 0.715979    
    ## Date2021-06-09  0.1009832  0.3070744 38.4234684   0.329 0.744050    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Avg_ds D2021-06-06 D2021-06-07
    ## Avg_dist    -0.184                               
    ## D2021-06-06 -0.685 -0.189                        
    ## D2021-06-07 -0.751 -0.285  0.652                 
    ## D2021-06-09 -0.770 -0.340  0.684       0.789

``` r
anova(RUACmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq  Mean Sq NumDF  DenDF F value Pr(>F)
    ## Avg_dist 0.00277 0.002767     1 45.412  0.0084 0.9275
    ## Date     0.63362 0.211208     3 39.185  0.6385 0.5947

``` r
# No significant effect of density or date on pollen loads in Rumex acetosella
```

### Rumex crispus (n=108)

``` r
RUCR <- stigd %>% 
  filter(Species=="Rumex crispus")
# summary(RUCR) 

#Taking a look at the data
hist(RUCR$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Rumex%20crispus-1.png)<!-- -->

``` r
hist(RUCR$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Rumex%20crispus-2.png)<!-- -->

``` r
hist(RUCR$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Rumex%20crispus-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=RUCR, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

    ## Warning: Removed 8 rows containing missing values (`geom_point()`).

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Rumex%20crispus-1.png)<!-- -->

``` r
RUCRmod <- lmer(Log_flw_pollen ~ Avg_dist + Date +  (1|Ind_ID), data=RUCR)
summary(RUCRmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: RUCR
    ## 
    ## REML criterion at convergence: 493
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -1.91800 -0.66208 -0.01091  0.65354  3.15354 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1064   0.3263  
    ##  Residual             0.8746   0.9352  
    ## Number of obs: 171, groups:  Ind_ID, 103
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   1.125951   0.291841 149.091567   3.858  0.00017 ***
    ## Avg_dist      0.003047   0.001134  57.201799   2.688  0.00939 ** 
    ## Date4164     -0.341347   0.394640 136.396666  -0.865  0.38858    
    ## Date4166      0.008090   0.340574 150.960201   0.024  0.98108    
    ## Date4167      0.138337   0.318359 136.566552   0.435  0.66459    
    ## Date4169      0.275546   0.320531 121.942561   0.860  0.39167    
    ## Date4170     -0.113510   0.347894 121.869907  -0.326  0.74477    
    ## Date4171      0.267063   0.410504  97.000492   0.651  0.51686    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4164 Dt4166 Dt4167 Dt4169 Dt4170
    ## Avg_dist -0.338                                          
    ## Date4164 -0.640 -0.045                                   
    ## Date4166 -0.756 -0.010  0.562                            
    ## Date4167 -0.818  0.017  0.600  0.696                     
    ## Date4169 -0.801 -0.017  0.597  0.691  0.739              
    ## Date4170 -0.736 -0.022  0.551  0.637  0.681  0.677       
    ## Date4171 -0.601 -0.085  0.470  0.541  0.576  0.575  0.530

``` r
anova(RUCRmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value   Pr(>F)   
    ## Avg_dist 6.3204  6.3204     1 57.202  7.2264 0.009393 **
    ## Date     4.5183  0.7531     6 78.434  0.8610 0.527395   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# Significant effect of density, no sig effect of date on pollen loads in Rumex crispus
```

### Schizacne purpurascens (n=154)

``` r
SCPU <- stigd %>% 
  filter(Species=="Schizacne purpurascens")
# summary(SCPU) 

#Taking a look at the data
hist(SCPU$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Schizacne%20purpurascens-1.png)<!-- -->

``` r
hist(SCPU$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Schizacne%20purpurascens-2.png)<!-- -->

``` r
hist(SCPU$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Schizacne%20purpurascens-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=SCPU, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Schizacne%20purpurascens-1.png)<!-- -->

``` r
SCPUmod <- lmer(Log_flw_pollen ~ Avg_dist + Date +  (1|Ind_ID), data=SCPU)
summary(SCPUmod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: SCPU
    ## 
    ## REML criterion at convergence: 878.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.6378 -0.7240  0.1132  0.7521  1.9514 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.04695  0.2167  
    ##  Residual             1.96323  1.4012  
    ## Number of obs: 246, groups:  Ind_ID, 154
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)   3.116985   0.507794 173.713196   6.138 5.50e-09 ***
    ## Avg_dist     -0.003987   0.003930 145.026241  -1.015  0.31200    
    ## Date4147      1.876202   0.377758 178.579276   4.967 1.58e-06 ***
    ## Date4148      0.995638   0.380244 178.698245   2.618  0.00959 ** 
    ## Date4149      0.644283   0.378063 179.013217   1.704  0.09008 .  
    ## Date4150      0.821383   0.413661 190.584022   1.986  0.04851 *  
    ## Date4153     -0.428716   0.428992 225.173572  -0.999  0.31869    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4147 Dt4148 Dt4149 Dt4150
    ## Avg_dist -0.766                                   
    ## Date4147 -0.575  0.025                            
    ## Date4148 -0.604  0.068  0.744                     
    ## Date4149 -0.591  0.047  0.747  0.745              
    ## Date4150 -0.563  0.073  0.684  0.682  0.685       
    ## Date4153 -0.511  0.029  0.658  0.655  0.658  0.603

``` r
anova(SCPUmod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Avg_dist   2.021  2.0209     1 145.03  1.0294     0.312    
    ## Date     115.584 23.1168     5 144.60 11.7749 1.462e-09 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Schizacne purpurascens
```

### Scirpus microcarpus (n=137)

``` r
SCMI <- stigd %>% 
  filter(Species=="Scirpus microcarpus")
# summary(SCMI) 

#Taking a look at the data
hist(SCMI$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Scirpus%20microcarpus-1.png)<!-- -->

``` r
hist(SCMI$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Scirpus%20microcarpus-2.png)<!-- -->

``` r
hist(SCMI$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Scirpus%20microcarpus-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=SCMI, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Scirpus%20microcarpus-1.png)<!-- -->

``` r
SCMImod <- lmer(Log_flw_pollen ~ Avg_dist + Date +  (1|Ind_ID), data=SCMI)
summary(SCMImod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + Date + (1 | Ind_ID)
    ##    Data: SCMI
    ## 
    ## REML criterion at convergence: 626.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5766 -0.5849 -0.0933  0.5327  2.8329 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.1482   0.3850  
    ##  Residual             0.4391   0.6627  
    ## Number of obs: 268, groups:  Ind_ID, 137
    ## 
    ## Fixed effects:
    ##               Estimate Std. Error         df t value Pr(>|t|)    
    ## (Intercept)  5.298e-01  1.330e-01  1.277e+02   3.984 0.000113 ***
    ## Avg_dist    -2.751e-04  1.376e-03  1.165e+02  -0.200 0.841925    
    ## Date4170     6.002e-01  1.622e-01  1.323e+02   3.699 0.000316 ***
    ## Date4171     1.188e+00  1.630e-01  1.296e+02   7.290 2.72e-11 ***
    ## Date4173     1.370e+00  1.619e-01  1.264e+02   8.463 5.54e-14 ***
    ## Date4174     1.451e+00  1.781e-01  1.323e+02   8.149 2.43e-13 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) Avg_ds Dt4170 Dt4171 Dt4173
    ## Avg_dist -0.483                            
    ## Date4170 -0.612 -0.034                     
    ## Date4171 -0.604 -0.045  0.514              
    ## Date4173 -0.606 -0.050  0.518  0.516       
    ## Date4174 -0.569 -0.007  0.469  0.467  0.471

``` r
anova(SCMImod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## Avg_dist  0.018  0.0175     1 116.52   0.040    0.8419    
    ## Date     47.148 11.7870     4 129.43  26.843 3.193e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Scirpus microcarpus
```

### Setaria viridis (n=45)

``` r
SEVI <- stigd %>% 
  filter(Species=="Setaria viridis")
# summary(SEVI) 

#Taking a look at the data
hist(SEVI$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Setaria%20viridis-1.png)<!-- -->

``` r
hist(SEVI$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Setaria%20viridis-2.png)<!-- -->

``` r
hist(SEVI$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Setaria%20viridis-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=SEVI, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Setaria%20viridis-1.png)<!-- -->

``` r
SEVImod <- lmer(Log_flw_pollen ~ Avg_dist +  (1|Ind_ID), data=SEVI)
summary(SEVImod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + (1 | Ind_ID)
    ##    Data: SEVI
    ## 
    ## REML criterion at convergence: 301.9
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -3.6380 -0.5518  0.1252  0.6380  1.6837 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.0422   0.2054  
    ##  Residual             0.5126   0.7159  
    ## Number of obs: 131, groups:  Ind_ID, 45
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error       df t value Pr(>|t|)    
    ## (Intercept)  4.13267    0.26453 40.36028  15.623   <2e-16 ***
    ## Avg_dist    -0.01766    0.01542 40.97475  -1.145    0.259    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr)
    ## Avg_dist -0.965

``` r
anova(SEVImod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
    ## Avg_dist 0.67172 0.67172     1 40.975  1.3105 0.2589

``` r
# No significant effect of density, sig effect of date on pollen loads in Setaria viridis
```

### Thalictrum dioicum (n=29)

- 59 total but 29 when just considering latest date. Looking at just
  2021-05-17 rather than 05-5 or 05-6 because I think those were too
  early, not representative of actual pollen capture.

``` r
THDI <- stigd %>% 
  filter(Species=="Thalictrum dioicum" & Date=="2021-05-17")
# summary(THDI) 
# str(THDI)
length(unique(THDI$Plant))
```

    ## [1] 29

``` r
#Taking a look at the data
hist(THDI$Flw_pollen) # response variable. looks very right skewed!
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Thalictrum%20dioicum-1.png)<!-- -->

``` r
hist(THDI$Log_flw_pollen) # much more symmetric
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Thalictrum%20dioicum-2.png)<!-- -->

``` r
hist(THDI$Avg_dist) 
```

![](anal-stig-density_files/figure-gfm/pre-fit%20Thalictrum%20dioicum-3.png)<!-- -->

``` r
# Look at plot of Log_flw_pollen and Avg_dist
ggplot(data=THDI, aes(x=Avg_dist, y=Log_flw_pollen, color=Date)) + 
  geom_point(size=3) + 
  scale_x_continuous(name="Average distance to 5 nearest neighbours") + 
  scale_y_continuous(name="Log stigmatic pollen load") + 
  theme_cs()
```

![](anal-stig-density_files/figure-gfm/flower%20height%20vs%20pollen%20load%20Thalictrum%20dioicum-1.png)<!-- -->

``` r
THDImod <- lmer(Log_flw_pollen ~ Avg_dist + (1|Ind_ID), data=THDI)
summary(THDImod)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Log_flw_pollen ~ Avg_dist + (1 | Ind_ID)
    ##    Data: THDI
    ## 
    ## REML criterion at convergence: 634.6
    ## 
    ## Scaled residuals: 
    ##      Min       1Q   Median       3Q      Max 
    ## -2.61229 -0.41784  0.09304  0.63613  2.64650 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Ind_ID   (Intercept) 0.6686   0.8177  
    ##  Residual             0.2429   0.4929  
    ## Number of obs: 369, groups:  Ind_ID, 29
    ## 
    ## Fixed effects:
    ##              Estimate Std. Error        df t value Pr(>|t|)    
    ## (Intercept)  2.514073   0.271025 26.961231   9.276 7.03e-10 ***
    ## Avg_dist    -0.013310   0.005364 26.549024  -2.482   0.0197 *  
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr)
    ## Avg_dist -0.822

``` r
anova(THDImod)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
    ## Avg_dist 1.4959  1.4959     1 26.549   6.158 0.01972 *
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
# No significant effect of density, sig effect of date on pollen loads in Thalictrum dioicum
```
