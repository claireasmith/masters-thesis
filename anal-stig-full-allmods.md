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

### post

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
  theme_cs(font = "sans", fontsize=18) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 18, face="italic"))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-full-allmods_files/figure-gfm/all%20stig%20length%20plots-1.png)<!-- -->

``` r
stigl %>% 
  ggplot(aes(x=Date, y=Log_flw_pollen)) + 
  geom_point() + 
  geom_smooth(method="lm", se=F) + 
  facet_wrap(~Species, scales = "free", ncol=4) + 
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500)) +
  scale_x_date(name = "Date")+ 
  theme_cs(font = "sans", fontsize=18) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 18, face="italic"),
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
  theme_cs(font = "sans", fontsize=18) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 18, face="italic"))
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
  scale_x_continuous(name = "1/average distance to nearest pollen-producing neighbour (1/cm)")+ 
  theme_cs(font = "sans", fontsize=18) + 
  theme(strip.background = element_rect(linewidth = NULL,
                                        linetype = NULL,
                                        colour = "white"),
        strip.text = element_text(size = 18, face="italic"))
```

    ## `geom_smooth()` using formula = 'y ~ x'

![](anal-stig-full-allmods_files/figure-gfm/all%20inv%20avg%20dist%20plots-1.png)<!-- -->

Run models all at once to get coefficient summaries easily:

``` r
#nest data by species to make it easier to handle
by_spp_stig <- stigl %>% group_by(Species, Sex_sys) %>% 
  nest() %>% 
  arrange(Sex_sys, Species)
```
