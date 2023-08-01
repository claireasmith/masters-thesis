Pollen production and size – summary
================
Claire Smith
2023-07-28

### Goal

Examine patterns in pollen size and number across wind-polllinated
flowering plant species.

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
prod <- prodfull %>% 
  filter(!is.na(Avg_diam) & !is.na(Avg_pol_anth)) # want only entries with both size and prod values
```

First I’ll take a look at a summary of the pollen size and number data -
means, standard deviations, and CVs so that I can compare how variable
species area relative to each other. I’ll remove any species with less
than 5 points.

``` r
prod_summary <- prod  %>% group_by(Species, source, Sex_sys) %>% 
  summarize(mean_diam = mean(Avg_diam, na.rm=T), 
            sd_diam = sd(Avg_diam, na.rm=T),
            cv_diam = sd_diam/mean_diam,
            min_diam = min(Avg_diam, na.rm=T),
            max_diam = max(Avg_diam, na.rm=T),
            mean_polanth = mean(Avg_pol_anth, na.rm=T),
            sd_polanth = sd(Avg_pol_anth, na.rm=T),
            cv_polanth = sd_polanth/mean_polanth,
            min_polanth = min(Avg_pol_anth, na.rm=T),
            max_polanth = max(Avg_pol_anth, na.rm=T),
            n=n()) %>%  
  filter(n>=3) %>% 
  droplevels() %>% 
  arrange(Sex_sys, Species) %>% 
  as.data.frame()
```

    ## `summarise()` has grouped output by 'Species', 'source'. You can override using
    ## the `.groups` argument.

``` r
# This will be made into a table - I'll add columns with mean +/- sd for that future table
prod_summary$anth_msd <- paste0(round(prod_summary$mean_polanth,0)," \U00B1 ",round(prod_summary$sd_polanth,0))
prod_summary$diam_msd <- paste0(round(prod_summary$mean_diam,1)," \U00B1 ",round(prod_summary$sd_diam,1))
```

I’ll save this summary data to a file.

``` r
prod_summary_filt <- prod_summary 
#write it to a file I will make a table with 
write.csv(prod_summary_filt, "summary-data/sizeprod-summary.csv", row.names = F, fileEncoding = "UTF-8")
```

Pollen size seems much less variable than pollen number. I’ll test
whether this difference is significant using a paired t test.

``` r
# how variable is pollen size/num across species? 
# Pollen size seems much less variable - I'll run a t test of mean cv for diameter vs mean cv for pollen per anther to determine if pollen size is significantly less variable than pollen production
mean(prod_summary$cv_diam)
```

    ## [1] 0.03678952

``` r
sd(prod_summary$cv_diam)
```

    ## [1] 0.02049357

``` r
mean(prod_summary$cv_polanth)
```

    ## [1] 0.3423091

``` r
sd(prod_summary$cv_polanth)
```

    ## [1] 0.1612152

``` r
cvtest <- prod_summary %>% select(cv_diam, cv_polanth)
t.test(x=cvtest$cv_polanth, y=cvtest$cv_diam, paired=T) # t test across all species
```

    ## 
    ##  Paired t-test
    ## 
    ## data:  cvtest$cv_polanth and cvtest$cv_diam
    ## t = 11.018, df = 30, p-value = 4.562e-12
    ## alternative hypothesis: true mean difference is not equal to 0
    ## 95 percent confidence interval:
    ##  0.2488911 0.3621480
    ## sample estimates:
    ## mean difference 
    ##       0.3055195

Pollen diameter is less variable than pollen size. The CV for pollen
diameter is significantly less than the CV for pollen production per
anther (t = 10.974, df = 30, p-value = 5.033e-12)
