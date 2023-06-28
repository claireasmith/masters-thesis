Testing for intraspecific size number tradeoff
================
Claire Smith
2023-06-28

Investigating the relationship between pollen size and number within
wind-pollinated species using linear regression: pollen number ~ pollen
volume

``` r
library(tidyverse)
library(ggplot2)
library(lme4)
library(lmerTest)
source("theme_cs.R")
```

``` r
prodfull <- read.csv("processed-data/size-prod-norep.csv", stringsAsFactors = T) # has individual-level pollen size and number data
# head(prodfull)
# str(prodfull)
# summary(prodfull)
```

``` r
prod_summary <- prodfull  %>% group_by(Species, source, Sex_sys) %>% 
  summarize(mean_diam = mean(Avg_diam, na.rm=T), 
            sd_diam = sd(Avg_diam, na.rm=T),
            cv_diam = sd_diam/mean_diam,
            mean_polanth = mean(Avg_pol_anth, na.rm=T),
            sd_polanth = sd(Avg_pol_anth, na.rm=T),
            cv_polanth = sd_polanth/mean_polanth,
            n=n()) %>%  
  filter(n>=5) %>% 
  droplevels() %>% 
  arrange(Sex_sys, Species) %>% 
  as.data.frame()
```

    ## `summarise()` has grouped output by 'Species', 'source'. You can override using
    ## the `.groups` argument.

``` r
prod_summary$anth_msd <- paste0(round(prod_summary$mean_polanth,0)," \U00B1 ",round(prod_summary$sd_polanth,0))
prod_summary$diam_msd <- paste0(round(prod_summary$mean_diam,0)," \U00B1 ",round(prod_summary$sd_diam,0))

prod_summary_filt <- prod_summary 
#write it to a file I will make a table with 
write.csv(prod_summary_filt, "processed-data/sizeprod-summary.csv", row.names = F, fileEncoding = "UTF-8")

# how variable is pollen vol/num across species? 
# Pollen size seems much less variable - I'll run a t test of mean cv for diameter vs mean cv for pollen per anther to determine if pollen size is significantly less variable than pollen production
mean(prod_summary$cv_diam)
```

    ## [1] 0.04286847

``` r
sd(prod_summary$cv_diam)
```

    ## [1] 0.01997849

``` r
mean(prod_summary$cv_polanth)
```

    ## [1] 0.3684211

``` r
sd(prod_summary$cv_polanth)
```

    ## [1] 0.1613425

``` r
cvtest <- prod_summary %>% select(cv_diam, cv_polanth)
t.test(x=cvtest$cv_polanth, y=cvtest$cv_diam, paired=T) # t test across all species
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  cvtest$cv_polanth and cvtest$cv_diam
    ## t = 10.683, df = 26, p-value = 5.267e-11
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  0.2629104 0.3881948
    ## sample estimates:
    ## mean difference 
    ##       0.3255526

Pollen diameter is less variable than pollen size. The CV for pollen
diameter is significantly less than the CV for pollen production per
anther.

``` r
head(prodfull)
```

    ##   source        Sex_sys       Species Site Ind Label Avg_diam  Sd_diam N_diam
    ## 1 JF2001 hermaphroditic     Acristata <NA>   1  <NA> 21.18655 6.260917  15230
    ## 2 JF2001 hermaphroditic     Acristata <NA>   2  <NA> 22.58325 5.632595  14830
    ## 3 JF2001 hermaphroditic        Adasyd <NA>   1  <NA> 24.94907 7.320320   3648
    ## 4 JF2001 hermaphroditic Elymus repens <NA> A01  <NA> 38.30692 2.009918  14823
    ## 5 JF2001 hermaphroditic Elymus repens <NA> A02  <NA> 37.85869 2.367872  15265
    ## 6 JF2001 hermaphroditic Elymus repens <NA> A03  <NA> 38.05908 2.905632  14821
    ##   Avg_area Avg_pol_anth Sd_pol_anth Anth_per_flw     count_method
    ## 1       NA     5076.667          NA            3 particle counter
    ## 2       NA     4943.333          NA            3 particle counter
    ## 3       NA     1824.000          NA            3 particle counter
    ## 4       NA     4941.000          NA            3 particle counter
    ## 5       NA     5088.333          NA            3 particle counter
    ## 6       NA     4940.333          NA            3 particle counter
    ##   Anth_per_sample Pol_flw Sd_pol_flw
    ## 1               3   15230         NA
    ## 2               3   14830         NA
    ## 3               2    5472         NA
    ## 4               3   14823         NA
    ## 5               3   15265         NA
    ## 6               3   14821         NA

``` r
summary(prodfull)
```

    ##     source              Sex_sys                  Species      Site    
    ##  CS2021: 16   dioecious     : 44   Hierochloe odorata: 67   FBF : 10  
    ##  JF2001:358   hermaphroditic:416   Festuca campestris: 60   FBF2:  6  
    ##  JF2004:252   monoecious    :166   Carex stipata     : 36   NA's:610  
    ##                                    Elymus innovatus  : 35             
    ##                                    Elymus repens     : 33             
    ##                                    Bromus inermis    : 30             
    ##                                    (Other)           :365             
    ##       Ind          Label        Avg_diam        Sd_diam           N_diam      
    ##  1      : 10   A1     :  2   Min.   :15.88   Min.   :0.1445   Min.   :     2  
    ##  2      :  9   A2     :  2   1st Qu.:22.89   1st Qu.:1.5343   1st Qu.:  3947  
    ##  a      :  9   A3     :  2   Median :25.17   Median :1.9706   Median :  9074  
    ##  3      :  8   A4     :  2   Mean   :26.55   Mean   :2.3518   Mean   : 14025  
    ##  4      :  7   A5     :  2   3rd Qu.:30.19   3rd Qu.:2.6167   3rd Qu.: 18551  
    ##  (Other):571   (Other):  6   Max.   :40.41   Max.   :9.4961   Max.   :141705  
    ##  NA's   : 12   NA's   :610   NA's   :17      NA's   :17                       
    ##     Avg_area      Avg_pol_anth      Sd_pol_anth      Anth_per_flw   
    ##  Min.   :379.7   Min.   :   46.0   Min.   : 46.19   Min.   : 3.000  
    ##  1st Qu.:442.3   1st Qu.:  627.7   1st Qu.:114.75   1st Qu.: 3.000  
    ##  Median :458.1   Median : 1859.1   Median :312.65   Median : 3.000  
    ##  Mean   :483.8   Mean   : 2774.8   Mean   :344.80   Mean   : 4.683  
    ##  3rd Qu.:530.2   3rd Qu.: 4123.3   3rd Qu.:495.09   3rd Qu.: 3.000  
    ##  Max.   :622.3   Max.   :15583.5   Max.   :884.84   Max.   :31.000  
    ##  NA's   :610                       NA's   :610      NA's   :4       
    ##            count_method Anth_per_sample      Pol_flw         Sd_pol_flw    
    ##  manual          : 16   Min.   :  2.000   Min.   :   138   Min.   : 230.9  
    ##  particle counter:610   1st Qu.:  3.000   1st Qu.:  2321   1st Qu.: 573.7  
    ##                         Median :  3.000   Median :  6915   Median :1563.3  
    ##                         Mean   :  9.551   Mean   : 11925   Mean   :1724.0  
    ##                         3rd Qu.:  8.000   3rd Qu.: 16447   3rd Qu.:2475.5  
    ##                         Max.   :109.000   Max.   :141705   Max.   :4424.2  
    ##                         NA's   :16        NA's   :4        NA's   :610

``` r
prod <- prodfull %>% 
  filter(!is.na(Avg_diam) & !is.na(Avg_pol_anth)) %>% # want only entries with both size and prod values
  # remove outliers
  filter(!(Species=="Phleum pratense"&Ind=="C11")) %>% 
  filter(!(Species=="Elymus innovatus"&Ind=="EI01"))
# summary(prod)
```

Now I want to build a model for each species. First I’ll go through each
species and look at the distributions of pollen size and number, as well
as if they have any other explanatory variables (like site) that might
be useful

``` r
#nest data by species to make it easier to handle
by_spp <- prod %>% group_by(Species, Sex_sys) %>% 
  nest() %>% 
  arrange(Sex_sys, Species) %>% 
  # Keep only species with at least 10 entries %>% 
  mutate(nind = purrr::map(data, nrow)) %>% 
  filter(nind>=10)
# by_spp # 23 species remaining

#inspect distributions of mean pollen diam and pollen production for each species - are they symmetric? any obvious outliers? 
#define function for this
checkfun <- function(df,i){
  spec = as.character(df[i,]$Species)
  n=length(df[i,]$data[[1]]$Avg_diam)
  hist(df[i,]$data[[1]]$Avg_diam, main=paste(spec, "n=",n), xlab="Mean pollen diameter")
  hist(log(df[i,]$data[[1]]$Avg_diam), main=paste(spec, "n=",n), xlab="Log mean pollen diameter")
  hist(df[i,]$data[[1]]$Avg_pol_anth, main=paste(spec, "n=",n), xlab="Mean pollen per anther")
  hist(log(df[i,]$data[[1]]$Avg_pol_anth), main=paste(spec, "n=",n), xlab="Log mean pollen per anther")
}

checkfun(by_spp,2)
```

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/prepare%20data%20for%20models-1.png)<!-- -->![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/prepare%20data%20for%20models-2.png)<!-- -->![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/prepare%20data%20for%20models-3.png)<!-- -->![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/prepare%20data%20for%20models-4.png)<!-- -->

``` r
#run through species and record in an excel file "within species size number tradeoff pre model fit sum.xlsx"
```

### LINEAR REGRESSION MODELS AND FIGURES

Run linear models and plot them

``` r
#define function to run model 
spp_model <- function(df) {
  lm(Avg_pol_anth ~ Avg_diam, data = df)
}
#run models and add model to by_spp dataframe
by_spp <- by_spp %>% 
  mutate(model = purrr::map(data, spp_model)) %>% 
  arrange(Sex_sys, Species) # arrange in order of sex sys, then alphabetical
#get model summary using glance() from "broom"
glance <- by_spp %>% 
  mutate(glance = purrr::map(model, broom::glance)) %>% 
  unnest(glance, .drop=T)
```

    ## Warning: The `.drop` argument of `unnest()` is deprecated as of tidyr 1.0.0.
    ## ℹ All list-columns are now preserved.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
# put this in a viewable format (remove nested data and model info)
glancetab <- glance %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species)
#save this in a new table without the ugly list elements, which I can print to a file
glancetab_write <- glance %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species) %>% apply(2, as.character) # coercing all columns to characters so that I can write this to a file
write.csv(glancetab_write, "processed-data/tradeoff-reg-model-sum.csv", row.names = F)
#take a look
print(glancetab, n=Inf)
```

    ## # A tibble: 23 × 15
    ## # Groups:   Species, Sex_sys [23]
    ##    Sex_sys  Species nind  r.squared adj.r.squared  sigma statistic p.value    df
    ##    <fct>    <fct>   <lis>     <dbl>         <dbl>  <dbl>     <dbl>   <dbl> <dbl>
    ##  1 dioecio… Rumex … <int>  0.000147      -0.0525   579.    0.00280 9.58e-1     1
    ##  2 dioecio… Thalic… <int>  0.0550         0.00999 1288.    1.22    2.81e-1     1
    ##  3 hermaph… Agropy… <int>  0.186          0.0840   887.    1.83    2.14e-1     1
    ##  4 hermaph… Bromus… <int>  0.221          0.193   2202.    7.95    8.73e-3     1
    ##  5 hermaph… Chenop… <int>  0.105          0.0601   133.    2.34    1.41e-1     1
    ##  6 hermaph… Elymus… <int>  0.214          0.189   2191.    8.72    5.87e-3     1
    ##  7 hermaph… Elymus… <int>  0.00859       -0.0234  1870.    0.269   6.08e-1     1
    ##  8 hermaph… Festuc… <int>  0.176          0.162   1783.   12.4     8.49e-4     1
    ##  9 hermaph… Festuc… <int>  0.439          0.392    726.    9.39    9.81e-3     1
    ## 10 hermaph… Hieroc… <int>  0.0102        -0.0104   788.    0.494   4.86e-1     1
    ## 11 hermaph… Koeler… <int>  0.249          0.216    925.    7.62    1.11e-2     1
    ## 12 hermaph… Phleum… <int>  0.00294       -0.0327   605.    0.0824  7.76e-1     1
    ## 13 hermaph… Planta… <int>  0.202          0.166    702.    5.58    2.74e-2     1
    ## 14 hermaph… Schiza… <int>  0.00359       -0.0961   104.    0.0360  8.53e-1     1
    ## 15 hermaph… Stipa … <int>  0.352          0.271    330.    4.35    7.05e-2     1
    ## 16 monoeci… Amaran… <int>  0.314          0.265    530.    6.42    2.39e-2     1
    ## 17 monoeci… Ambros… <int>  0.0126        -0.0394   395.    0.242   6.28e-1     1
    ## 18 monoeci… Carex … <int>  0.214          0.164     91.6   4.34    5.35e-2     1
    ## 19 monoeci… Carex … <int>  0.0531        -0.00613  153.    0.896   3.58e-1     1
    ## 20 monoeci… Carex … <int>  0.129          0.0774   121.    2.51    1.31e-1     1
    ## 21 monoeci… Carex … <int>  0.0389         0.0106    72.5   1.38    2.49e-1     1
    ## 22 monoeci… Rumex … <int>  0.0356        -0.0180   168.    0.665   4.26e-1     1
    ## 23 monoeci… Scirpu… <int>  0.00334       -0.0797    66.1   0.0403  8.44e-1     1
    ## # ℹ 6 more variables: logLik <dbl>, AIC <dbl>, BIC <dbl>, deviance <dbl>,
    ## #   df.residual <int>, nobs <int>

``` r
#similar - but using broom::tidy(), gives me model coefficients
tidy <- by_spp %>% 
  mutate(tidy = purrr::map(model, broom::tidy)) %>% 
  unnest(tidy, .drop=T)
```

    ## Warning: The `.drop` argument of `unnest()` is deprecated as of tidyr 1.0.0.
    ## ℹ All list-columns are now preserved.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
tidytab <- tidy %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species)
print(tidytab, n=46)
```

    ## # A tibble: 46 × 8
    ## # Groups:   Species, Sex_sys [23]
    ##    Sex_sys        Species       nind  term  estimate std.error statistic p.value
    ##    <fct>          <fct>         <lis> <chr>    <dbl>     <dbl>     <dbl>   <dbl>
    ##  1 dioecious      Rumex acetos… <int> (Int…   2.30e3    4340.     0.530  6.02e-1
    ##  2 dioecious      Rumex acetos… <int> Avg_…  -1.15e1     217.    -0.0529 9.58e-1
    ##  3 dioecious      Thalictrum d… <int> (Int…   1.16e4    7467.     1.55   1.36e-1
    ##  4 dioecious      Thalictrum d… <int> Avg_…  -4.97e2     449.    -1.11   2.81e-1
    ##  5 hermaphroditic Agropyron tr… <int> (Int…   1.51e4    7157.     2.10   6.85e-2
    ##  6 hermaphroditic Agropyron tr… <int> Avg_…  -2.56e2     190.    -1.35   2.14e-1
    ##  7 hermaphroditic Bromus inerm… <int> (Int…   2.07e4    5339.     3.87   5.92e-4
    ##  8 hermaphroditic Bromus inerm… <int> Avg_…  -4.49e2     159.    -2.82   8.73e-3
    ##  9 hermaphroditic Chenopodium … <int> (Int…   2.45e3    1259.     1.95   6.56e-2
    ## 10 hermaphroditic Chenopodium … <int> Avg_…  -7.96e1      52.0   -1.53   1.41e-1
    ## 11 hermaphroditic Elymus innov… <int> (Int…  -1.63e4    8087.    -2.01   5.26e-2
    ## 12 hermaphroditic Elymus innov… <int> Avg_…   6.78e2     230.     2.95   5.87e-3
    ## 13 hermaphroditic Elymus repens <int> (Int…   9.49e2    7813.     0.121  9.04e-1
    ## 14 hermaphroditic Elymus repens <int> Avg_…   1.09e2     210.     0.518  6.08e-1
    ## 15 hermaphroditic Festuca camp… <int> (Int…   2.16e4    4295.     5.03   5.01e-6
    ## 16 hermaphroditic Festuca camp… <int> Avg_…  -5.22e2     148.    -3.52   8.49e-4
    ## 17 hermaphroditic Festuca prat… <int> (Int…   1.87e4    5302.     3.53   4.15e-3
    ## 18 hermaphroditic Festuca prat… <int> Avg_…  -5.48e2     179.    -3.06   9.81e-3
    ## 19 hermaphroditic Hierochloe o… <int> (Int…   7.71e2    1346.     0.573  5.70e-1
    ## 20 hermaphroditic Hierochloe o… <int> Avg_…   4.17e1      59.4    0.703  4.86e-1
    ## 21 hermaphroditic Koeleria cri… <int> (Int…   1.23e4    3595.     3.42   2.35e-3
    ## 22 hermaphroditic Koeleria cri… <int> Avg_…  -4.50e2     163.    -2.76   1.11e-2
    ## 23 hermaphroditic Phleum prate… <int> (Int…   2.24e3    1891.     1.18   2.47e-1
    ## 24 hermaphroditic Phleum prate… <int> Avg_…  -1.81e1      63.2   -0.287  7.76e-1
    ## 25 hermaphroditic Plantago lan… <int> (Int…   8.92e3    2618.     3.41   2.52e-3
    ## 26 hermaphroditic Plantago lan… <int> Avg_…  -2.91e2     123.    -2.36   2.74e-2
    ## 27 hermaphroditic Schizachne p… <int> (Int…  -5.53e1    1804.    -0.0307 9.76e-1
    ## 28 hermaphroditic Schizachne p… <int> Avg_…   1.17e1      61.6    0.190  8.53e-1
    ## 29 hermaphroditic Stipa columb… <int> (Int…  -4.24e3    2360.    -1.80   1.10e-1
    ## 30 hermaphroditic Stipa columb… <int> Avg_…   1.72e2      82.6    2.09   7.05e-2
    ## 31 monoecious     Amaranthus r… <int> (Int…   7.78e3    2030.     3.83   1.83e-3
    ## 32 monoecious     Amaranthus r… <int> Avg_…  -2.08e2      81.9   -2.53   2.39e-2
    ## 33 monoecious     Ambrosia art… <int> (Int…   2.61e3    3746.     0.697  4.94e-1
    ## 34 monoecious     Ambrosia art… <int> Avg_…  -1.03e2     209.    -0.492  6.28e-1
    ## 35 monoecious     Carex commun… <int> (Int…   4.95e1     268.     0.184  8.56e-1
    ## 36 monoecious     Carex commun… <int> Avg_…   2.28e1      11.0    2.08   5.35e-2
    ## 37 monoecious     Carex hirtif… <int> (Int…   2.92e3    2386.     1.23   2.38e-1
    ## 38 monoecious     Carex hirtif… <int> Avg_…  -8.96e1      94.6   -0.947  3.58e-1
    ## 39 monoecious     Carex pedunc… <int> (Int…   2.09e3     991.     2.11   5.03e-2
    ## 40 monoecious     Carex pedunc… <int> Avg_…  -6.37e1      40.2   -1.58   1.31e-1
    ## 41 monoecious     Carex stipata <int> (Int…   8.76e2     395.     2.22   3.32e-2
    ## 42 monoecious     Carex stipata <int> Avg_…  -1.88e1      16.0   -1.17   2.49e-1
    ## 43 monoecious     Rumex crispus <int> (Int…   2.11e3    1657.     1.27   2.19e-1
    ## 44 monoecious     Rumex crispus <int> Avg_…  -5.62e1      68.9   -0.815  4.26e-1
    ## 45 monoecious     Scirpus micr… <int> (Int…   3.62e2     965.     0.376  7.14e-1
    ## 46 monoecious     Scirpus micr… <int> Avg_…  -8.49e0      42.3   -0.201  8.44e-1

``` r
tidytab_write <- tidy %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species) %>% apply(2, as.character)
write.csv(tidytab_write, "processed-data/tradeoff-reg-model-coefs.csv", row.names = F)
#take a look
# print(tidytab, n=Inf)

#take a look at r.squared values
# glance %>% 
#   ggplot(aes(spp, r.squared)) + 
#   geom_point() + 
#   theme(axis.text = element_text(angle=90))
# take a look at all data
# prod %>% filter(spp %in% keepspp) %>% 
#   ggplot(aes(x=polvol, y=polanth)) +
#   geom_point() +
#   facet_wrap(.~spp, scales="free")
```

``` r
#take only species with significant relationships to plot with linear model overtop
sigsub <- glance %>% filter(p.value<0.05) #7 species
#Plot relationship between size and number for the 7 species with a significant relationship
# sigsub$Species
```

``` r
# Bromus inermis
sigsub[which(sigsub$Species=="Bromus inermis"),]$data %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/Bromus%20inermis%20tradeoff-1.png)<!-- -->

``` r
# Elymus innovatus
sigsub[which(sigsub$Species=="Elymus innovatus"),]$data %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/Elymus%20innovatus%20tradeoff-1.png)<!-- -->

``` r
#Festuca campestris
sigsub[which(sigsub$Species=="Festuca campestris"),]$data %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/Festuca%20campestris%20tradeoff-1.png)<!-- -->

``` r
# Festuca pratensis
sigsub[which(sigsub$Species=="Festuca pratensis"),]$data %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/Festuca%20pratensis%20tradeoff-1.png)<!-- -->

``` r
# Koeleria cristata
sigsub[which(sigsub$Species=="Koeleria cristata"),]$data %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/Koeleria%20cristata%20tradeoff-1.png)<!-- -->

``` r
# Plantago lanceolata
sigsub[which(sigsub$Species=="Plantago lanceolata"),]$data %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/Plantago%20lanceolata%20tradeoff-1.png)<!-- -->

``` r
#Amaranthus retroflexus
sigsub[which(sigsub$Species=="Amaranthus retroflexus"),]$data %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Avg_diam, y=Avg_pol_anth)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  scale_y_continuous(name="Pollen production per anther") + 
  scale_x_continuous(name=expression(paste("Pollen diameter (", mu, "m",")"))) + 
  theme_cs(font = "sans")
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-sizenum-tradeoff-within-spp_files/figure-gfm/Amaranthus%20retroflexus%20tradeoff-1.png)<!-- -->

Found a significantly negative relationship between pollen volume and
number (i.e. size-number trade-off) in 4/15 species (only taking species
with at least 15 entries for pollen size/number). Found a significantly
positive relationship between size and number in 1 species (Leymus
innovatus).
