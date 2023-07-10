Untitled
================
Claire Smith
2023-07-09

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
# summary(stig)
```

### Prepare data

``` r
stigl <- stig %>% 
  mutate(Log_flw_pollen=log(Flw_pollen + 1)) %>% 
  rowwise %>% 
  mutate(Avg_dist = mean(c_across(c('D1', 'D2', 'D3', 'D4', 'D5')), na.rm=TRUE),
         Inv_avg_dist=1/Avg_dist) %>% 
  # Take only entries stigma length data
  filter(!is.na(Flw_pollen) & !is.na(Infl_max) & !is.na(Avg_dist) &!is.na(Stigma_length) & !is.na(Date))
# remove outlier identified in plots in anal-stig-full.Rmd...
# CAHI[which(CAHI$Stigma_length<1),] # small stigma has length 0.2 (most are between 2-3)
stigl2 <- stigl %>% filter(Species != "Carex hirtifolia" | 
                            (Species == "Carex hirtifolia" & Stigma_length>0.5))
```

``` r
stigl %>% 
  ggplot(aes(x=Stigma_length, y=Log_flw_pollen)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  facet_wrap(~Species, scales = "free", ncol=3) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Stigma length (cm)") +
  theme_cs(font = "sans", fontsize=16) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 14, face="italic"))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-full-allmods_files/figure-gfm/all%20stig%20length%20plots-1.png)<!-- -->

``` r
stigl %>% 
  ggplot(aes(x=Date, y=Log_flw_pollen)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  facet_wrap(~Species, scales = "free", ncol=3) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_date(name = "Date")+ 
  theme_cs(font = "sans", fontsize=16) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 14, face="italic"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-full-allmods_files/figure-gfm/all%20date%20plots-1.png)<!-- -->

``` r
stigl %>% 
  ggplot(aes(x=Infl_max, y=Log_flw_pollen)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  facet_wrap(~Species, scales = "free", ncol=3) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Max plant height (cm)")+ 
  theme_cs(font = "sans", fontsize=16) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 14, face="italic"))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-full-allmods_files/figure-gfm/all%20infl%20max%20plots-1.png)<!-- -->

``` r
stigl %>% 
  ggplot(aes(x=Inv_avg_dist, y=Log_flw_pollen)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  facet_wrap(~Species, scales = "free", ncol=3) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_continuous(name = "Inverse average distance to nearest pollen-producing neighbour (1/cm)")+ 
  theme_cs(font = "sans", fontsize=16) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 14, face="italic"))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-full-allmods_files/figure-gfm/all%20inv%20avg%20dist%20plots-1.png)<!-- -->

``` r
# stigl %>% 
#   filter(Species=="Rumex acetosella") %>% 
#   ggplot(aes(x=Avg_dist, y=Flw_pollen)) + 
#   geom_point() + 
#   geom_smooth(method="lm", se=F, formula = log(y+1) ~ 1/x) + 
#   facet_wrap(~Species, scales = "free", ncol=3) + 
#   # scale_y_continuous(name="Stigmatic pollen load",
#   #                    breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
#   #                    labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
#   scale_x_continuous(name = "Inverse average distance to nearest pollen-producing neighbour (1/cm)")+ 
#   theme_cs(font = "sans", fontsize=16) + 
#   theme(strip.background = element_rect(linewidth = NULL,
#                                         linetype = NULL,
#                                         colour = "white"),
#         strip.text = element_text(size = 14, face="italic"))
```

Run models all at once to get coefficient summaries easily:

``` r
#nest data by species to make it easier to handle
by_spp <- stigl %>% group_by(Species, Sex_sys) %>% 
  nest() %>% 
  arrange(Sex_sys, Species)
```

Now I’ll run linear models predicting pollen production by pollen size,
and pull out the ones with a significant relationship.

``` r
# , cols.print=14
library(broom.mixed)
#define function to run model 
spp_model <- function(df) {
  lmer(Log_flw_pollen ~ Infl_max + Stigma_length + Date + Inv_avg_dist + (1|Plant), data=df)
}

#run models and add model to by_spp dataframe
by_spp <- by_spp %>% 
  mutate(model = purrr::map(data, spp_model)) %>% 
  arrange(Sex_sys, Species) # arrange in order of sex sys, then alphabetical
```

    ## boundary (singular) fit: see help('isSingular')
    ## boundary (singular) fit: see help('isSingular')

    ## Warning: There were 10 warnings in `mutate()`.
    ## The first warning was:
    ## ℹ In argument: `model = purrr::map(data, spp_model)`.
    ## ℹ In group 1: `Species = Ambrosia artemisiifolia`, `Sex_sys = monoecious`.
    ## Caused by warning:
    ## ! Some predictor variables are on very different scales: consider rescaling
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 9 remaining warnings.

``` r
#get model summary using glance() from "broom"
glance <- by_spp %>% 
  mutate(glance = purrr::map(model, broom.mixed::glance)) %>% #need to use broom.mixed::glance for mixed model
  unnest(glance)
# put this in a viewable format (remove nested data and model info)
glancetab <- glance %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species)

#save this in a new table without the ugly list elements, which I can print to a file
glancetab_write <- glance %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species) %>% apply(2, as.character) # coercing all columns to characters so that I can write this to a file
write.csv(glancetab_write, "processed-data/stigfull-reg-model-sum.csv", row.names = F)
#take a look
# print(glancetab, n=Inf)

#similar - but using broom::tidy(), gives me model coefficients
tidy <- by_spp %>% 
  mutate(tidy = purrr::map(model, broom.mixed::tidy, effects="fixed")) %>% 
  unnest(tidy)
tidytab <- tidy %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species)
print(tidytab, n=46)
```

    ## # A tibble: 65 × 9
    ## # Groups:   Species, Sex_sys [13]
    ##    Species     Sex_sys effect term  estimate std.error statistic     df  p.value
    ##    <fct>       <fct>   <chr>  <chr>    <dbl>     <dbl>     <dbl>  <dbl>    <dbl>
    ##  1 Rumex acet… dioeci… fixed  (Int… -7.57e+2   9.43e+1    -8.03  255.   3.49e-14
    ##  2 Rumex acet… dioeci… fixed  Infl… -1.37e-2   8.58e-3    -1.59   39.6  1.19e- 1
    ##  3 Rumex acet… dioeci… fixed  Stig…  1.31e+0   3.94e-1     3.33  254.   9.93e- 4
    ##  4 Rumex acet… dioeci… fixed  Date   6.03e-2   7.50e-3     8.03  255.   3.58e-14
    ##  5 Rumex acet… dioeci… fixed  Inv_…  1.47e+1   5.41e+0     2.72   35.3  1.00e- 2
    ##  6 Thalictrum… dioeci… fixed  (Int… -1.81e+3   3.96e+2    -4.57   89.9  1.56e- 5
    ##  7 Thalictrum… dioeci… fixed  Infl…  2.18e-2   1.49e-2     1.47   12.6  1.67e- 1
    ##  8 Thalictrum… dioeci… fixed  Stig…  3.41e-1   1.36e-1     2.51   89.9  1.37e- 2
    ##  9 Thalictrum… dioeci… fixed  Date   1.44e-1   3.16e-2     4.57   89.9  1.56e- 5
    ## 10 Thalictrum… dioeci… fixed  Inv_…  6.28e+0   1.87e+1     0.336  11.6  7.43e- 1
    ## 11 Chenopodiu… hermap… fixed  (Int…  3.18e+2   6.60e+2     0.482 144.   6.30e- 1
    ## 12 Chenopodiu… hermap… fixed  Infl… -3.79e-3   3.79e-3    -1.00  144.   3.19e- 1
    ## 13 Chenopodiu… hermap… fixed  Stig…  2.58e+0   5.74e-1     4.49  144.   1.46e- 5
    ## 14 Chenopodiu… hermap… fixed  Date  -2.52e-2   5.23e-2    -0.482 144.   6.31e- 1
    ## 15 Chenopodiu… hermap… fixed  Inv_… -2.88e+0   4.44e+0    -0.648 144.   5.18e- 1
    ## 16 Plantago l… hermap… fixed  (Int… -1.81e+3   2.59e+2    -7.00  254.   2.21e-11
    ## 17 Plantago l… hermap… fixed  Infl… -8.83e-3   1.84e-2    -0.479  26.7  6.36e- 1
    ## 18 Plantago l… hermap… fixed  Stig…  6.06e-1   1.32e-1     4.59  265.   6.96e- 6
    ## 19 Plantago l… hermap… fixed  Date   1.44e-1   2.06e-2     7.02  254.   2.04e-11
    ## 20 Plantago l… hermap… fixed  Inv_… -1.26e+1   1.51e+1    -0.835  26.8  4.11e- 1
    ## 21 Schizacne … hermap… fixed  (Int…  3.08e+3   6.31e+2     4.88  221.   2.00e- 6
    ## 22 Schizacne … hermap… fixed  Infl…  2.25e-2   2.48e-2     0.909  27.8  3.71e- 1
    ## 23 Schizacne … hermap… fixed  Stig…  1.07e-1   3.05e-1     0.350 234.   7.27e- 1
    ## 24 Schizacne … hermap… fixed  Date  -2.45e-1   5.02e-2    -4.88  221.   2.02e- 6
    ## 25 Schizacne … hermap… fixed  Inv_…  6.68e+1   4.83e+1     1.38   29.5  1.77e- 1
    ## 26 Ambrosia a… monoec… fixed  (Int… -1.49e+2   2.84e+2    -0.525 212.   6.00e- 1
    ## 27 Ambrosia a… monoec… fixed  Infl… -7.06e-4   3.97e-3    -0.178  23.9  8.60e- 1
    ## 28 Ambrosia a… monoec… fixed  Stig…  2.18e-1   1.38e-1     1.58  209.   1.15e- 1
    ## 29 Ambrosia a… monoec… fixed  Date   1.21e-2   2.25e-2     0.538 212.   5.91e- 1
    ## 30 Ambrosia a… monoec… fixed  Inv_…  4.67e+0   4.34e+0     1.07   24.3  2.93e- 1
    ## 31 Carex comm… monoec… fixed  (Int… -3.97e+3   4.54e+2    -8.75  185.   1.27e-15
    ## 32 Carex comm… monoec… fixed  Infl… -1.87e-2   3.72e-2    -0.502  19.5  6.21e- 1
    ## 33 Carex comm… monoec… fixed  Stig…  4.59e-1   1.49e-1     3.08  194.   2.35e- 3
    ## 34 Carex comm… monoec… fixed  Date   3.17e-1   3.62e-2     8.75  185.   1.27e-15
    ## 35 Carex comm… monoec… fixed  Inv_…  6.33e+0   1.22e+1     0.517  18.9  6.11e- 1
    ## 36 Carex hirt… monoec… fixed  (Int… -1.22e+3   1.11e+3    -1.10  195.   2.72e- 1
    ## 37 Carex hirt… monoec… fixed  Infl…  5.67e-2   4.02e-2     1.41   21.7  1.73e- 1
    ## 38 Carex hirt… monoec… fixed  Stig…  4.45e-1   2.63e-1     1.69  196.   9.30e- 2
    ## 39 Carex hirt… monoec… fixed  Date   9.74e-2   8.84e-2     1.10  195.   2.72e- 1
    ## 40 Carex hirt… monoec… fixed  Inv_…  7.76e+0   7.76e+0     1.00   20.3  3.29e- 1
    ## 41 Carex pedu… monoec… fixed  (Int… -5.85e+3   8.98e+2    -6.52   68.5  9.93e- 9
    ## 42 Carex pedu… monoec… fixed  Infl… -2.89e-2   8.13e-2    -0.356   9.25 7.30e- 1
    ## 43 Carex pedu… monoec… fixed  Stig…  1.33e-1   2.00e-1     0.666  70.8  5.08e- 1
    ## 44 Carex pedu… monoec… fixed  Date   4.67e-1   7.16e-2     6.52   68.5  1.00e- 8
    ## 45 Carex pedu… monoec… fixed  Inv_…  1.43e+1   3.46e+1     0.412  10.8  6.88e- 1
    ## 46 Carex plan… monoec… fixed  (Int… -6.45e+2   2.14e+2    -3.02   85.4  3.37e- 3
    ## # ℹ 19 more rows

``` r
tidytab_write <- tidy %>% select(-c(data, model)) %>% 
  arrange(Sex_sys, Species) 
# %>% apply(2, as.character)
write.csv(tidytab_write, "processed-data/stigfull-reg-model-coefs.csv", row.names = F)
#take a look
# print(tidytab, n=Inf)
```
