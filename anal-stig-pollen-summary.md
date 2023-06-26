Overall stigmatic pollen capture
================
Claire Smith
2023-06-26

``` r
# Load packages
library(tidyverse)
library(ggplot2)
library(ggridges)
```

    ## Warning: package 'ggridges' was built under R version 4.3.1

``` r
# Source files
source("clean-stig-01.R")
source("clean-stig-02.R")
source("theme_cs.R")
```

``` r
# Load data
stig <- read.csv("processed-data/stig-no-rep-spp.csv", stringsAsFactors = T)
```

``` r
# head(stig)
# summary(stig)
stig[which(is.na(stig$Flw_pollen)),] # some NAs from when the data was joined - some individuals had 
```

    ##                       Species       Date   Site Plant Flower Flw_pollen
    ## 229  Dichanthelium implicatum 2021-06-24    KME   S12      C         NA
    ## 264       Plantago lanceolata 2021-07-31     KM     6   <NA>         NA
    ## 276       Plantago lanceolata 2021-07-31     KM    10      3         NA
    ## 288       Plantago lanceolata 2021-07-31     KM    14      3         NA
    ## 290       Plantago lanceolata 2021-07-31     KM    14      2         NA
    ## 291       Plantago lanceolata 2021-07-31     KM    14      3         NA
    ## 303       Plantago lanceolata 2021-08-18    FBF     1      3         NA
    ## 305       Plantago lanceolata 2021-08-18    FBF     2      2         NA
    ## 306       Plantago lanceolata 2021-08-18    FBF     2      3         NA
    ## 309       Plantago lanceolata 2021-08-18    FBF     3      3         NA
    ## 360           Phleum pratense 2021-07-05    KME   S11      A         NA
    ## 369           Phleum pratense 2021-07-05    KME   S14      B         NA
    ## 423           Setaria viridis 2021-07-16    FBF    13      B         NA
    ## 470           Setaria viridis 2021-07-16    FBF    28      A         NA
    ## 1092       Thalictrum dioicum 2021-05-07     BI    B2    2-8         NA
    ## 1598      Plantago lanceolata  20-Aug-21    QMP     6      3         NA
    ## 1604      Plantago lanceolata  20-Aug-21    QMP     8      3         NA
    ## 1671           Plantago major       <NA> Biosci     3      A         NA
    ## 1672           Plantago major       <NA> Biosci    10      A         NA
    ##      Stigmas_per_flw Flower_height Infl_max Infl_min    D1   D2    D3    D4
    ## 229                2          17.5    18.50     16.0  7.00 10.0 18.00   9.0
    ## 264                1            NA    18.10     17.0 12.00 42.0 29.00  42.0
    ## 276                1            NA    19.40     17.8 13.00 26.0 24.00  22.0
    ## 288                1            NA    16.80     15.7 22.00 27.0 36.00  45.0
    ## 290                1            NA    16.80     15.7 22.00 27.0 36.00  45.0
    ## 291                1            NA    16.80     15.7 22.00 27.0 36.00  45.0
    ## 303                1            NA    33.50     33.0  8.00  6.0  6.00  10.0
    ## 305                1            NA    22.00     21.5  2.00  9.0 18.00  14.0
    ## 306                1            NA    22.00     21.5  2.00  9.0 18.00  14.0
    ## 309                1            NA    32.80     31.4  9.00 14.0 32.00  36.0
    ## 360                2          71.0    75.00     70.5 62.00 68.0 76.00  55.0
    ## 369                2          76.5    78.00     72.5 24.00 48.0 64.00  60.0
    ## 423                2          41.5    44.50     41.5 10.00 14.0 15.00  11.0
    ## 470                2          33.0    38.50     32.0 13.00 22.0 23.00  31.0
    ## 1092               1            NA    57.15       NA 40.64 76.2 58.42 114.3
    ## 1598               1            NA    26.20     24.2  7.00  8.0 12.00  16.0
    ## 1604               1            NA    19.10     18.0  3.00 15.0 18.00  26.0
    ## 1671               1           4.0     8.00      3.5  8.00 11.0 11.00   7.0
    ## 1672               1           5.0     7.50      5.0 28.00 32.0 32.00  23.0
    ##         D5 Stigma_length source        Sex_sys
    ## 229   15.0            NA CS2021 hermaphroditic
    ## 264   66.0            NA CS2021 hermaphroditic
    ## 276   30.0            NA CS2021 hermaphroditic
    ## 288   37.0            NA CS2021 hermaphroditic
    ## 290   37.0            NA CS2021 hermaphroditic
    ## 291   37.0            NA CS2021 hermaphroditic
    ## 303   11.0            NA CS2021 hermaphroditic
    ## 305   22.0            NA CS2021 hermaphroditic
    ## 306   22.0            NA CS2021 hermaphroditic
    ## 309   37.0            NA CS2021 hermaphroditic
    ## 360   75.0            NA CS2021 hermaphroditic
    ## 369   58.0            NA CS2021 hermaphroditic
    ## 423   14.0            NA CS2021 hermaphroditic
    ## 470   32.0            NA CS2021 hermaphroditic
    ## 1092 584.2            NA CS2021      dioecious
    ## 1598  18.0            NA CS2021 hermaphroditic
    ## 1604  17.0            NA CS2021 hermaphroditic
    ## 1671  14.0            NA CS2021 hermaphroditic
    ## 1672  28.0            NA CS2021 hermaphroditic

Visualize distributions of pollen loads across species

``` r
stig_ridges <- stig %>% 
  # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(desc(Sex_sys), desc(Species)) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot(aes(x=log(Flw_pollen+1), y=Species, fill=Sex_sys)) +
  geom_density_ridges2() + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#FC8D62"),
                    labels=c("Dioecious", "Monoecious", "Hermaphroditic")) +
  scale_x_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,25,50,100,200,500,1000, 2500)+1),
                     labels = c(0,1,2,5,10,25,50,100,200,500,1000, 2500),
                     limits=c(0,NA)) +
  scale_y_discrete() +
  ylab(label="") + 
  guides(fill="none") + 
  theme_cs(font="sans", fontsize=18) + 
  theme(axis.text.y = element_text(face = "italic"),
        axis.ticks.y = element_line())

stig_ridges
```

    ## Picking joint bandwidth of 0.358

    ## Warning: Removed 19 rows containing non-finite values
    ## (`stat_density_ridges()`).

![](anal-stig-pollen-summary_files/figure-gfm/stigma%20ridgeplots-1.png)<!-- -->

``` r
stig %>% 
    # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(Sex_sys, Species) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot() + 
  geom_boxplot(aes(x=Species, y=log(Flw_pollen + 1), fill=Sex_sys)) + 
  
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,50,100,200,500,1000)+1),
                     labels = c(0,1,2,5,10,50,100,200,500,1000),
                     limits=c(0,NA)) + 
  scale_x_discrete(name = "") + 
  
  guides(fill="none") + 
  scale_fill_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                    values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  theme_cs(fontsize=18, font="sans") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic"))
```

    ## Warning: Removed 19 rows containing non-finite values (`stat_boxplot()`).

![](anal-stig-pollen-summary_files/figure-gfm/stigma%20box%20plots-1.png)<!-- -->

``` r
stig_box_ss <- stig %>% 
  # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(Sex_sys, Species) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot() + 
  geom_boxplot(aes(x=Sex_sys, y=log(Flw_pollen + 1), fill=Sex_sys),
               width=0.4) + 
  
  scale_y_continuous(name="Stigmatic pollen load",
                     breaks = log(c(0,1,2,5,10,50,100,200,500,1000, 1500)+1),
                     labels = c(0,1,2,5,10,50,100,200,500,1000, 1500),
                     limits=c(0,NA)) + 
  scale_x_discrete(name = "", labels = NULL) + 
  
  guides(fill=guide_legend(title="Sex system")) + 
  scale_fill_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                    values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  theme_cs(fontsize=18, font="sans") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "right")

stig_box_ss
```

    ## Warning: Removed 19 rows containing non-finite values (`stat_boxplot()`).

![](anal-stig-pollen-summary_files/figure-gfm/sex%20sys%20stigma%20box%20plots-1.png)<!-- -->

``` r
pollen_load_tab <- stig %>% group_by(source, Sex_sys, Species, Plant) %>% 
  group_by(Sex_sys, Species, source) %>% 
  summarize(mean_poll = mean(Flw_pollen, na.rm=T),
            sd_poll = sd(Flw_pollen, na.rm=T),
            n_poll_count = n(),
            se_poll = sd_poll/sqrt(n_poll_count),
            CV_poll = sd_poll/mean_poll,
            max_poll = max(Flw_pollen, na.rm=T),
            min_poll = min(Flw_pollen, na.rm=T))
```

    ## `summarise()` has grouped output by 'Sex_sys', 'Species'. You can override
    ## using the `.groups` argument.

``` r
n_plants_stig <- stig %>% group_by(Sex_sys, Species, source, Plant) %>% 
  summarize(g=n()) %>%
  group_by(Sex_sys, Species, source) %>% 
  summarize(nplants = n())
```

    ## `summarise()` has grouped output by 'Sex_sys', 'Species', 'source'. You can
    ## override using the `.groups` argument.
    ## `summarise()` has grouped output by 'Sex_sys', 'Species'. You can override
    ## using the `.groups` argument.

``` r
pollen_load_tab$stig_nplants = n_plants_stig$nplants
# names(pollen_load_tab)
pollen_load_tab
```

    ## # A tibble: 23 × 11
    ## # Groups:   Sex_sys, Species [23]
    ##    Sex_sys        Species  source mean_poll sd_poll n_poll_count se_poll CV_poll
    ##    <fct>          <fct>    <fct>      <dbl>   <dbl>        <int>   <dbl>   <dbl>
    ##  1 dioecious      Rumex a… CS2021      2.26    2.43          117   0.225   1.08 
    ##  2 dioecious      Thalict… CS2021      3.69    6.48          921   0.213   1.75 
    ##  3 hermaphroditic Bromus … JF2001     50.6    75.7           391   3.83    1.50 
    ##  4 hermaphroditic Chenopo… JF2004      3.09    3.73          150   0.304   1.21 
    ##  5 hermaphroditic Dichant… CS2021     14.6    16.3            64   2.04    1.12 
    ##  6 hermaphroditic Dichant… CS2021     38.0    61.8            20  13.8     1.62 
    ##  7 hermaphroditic Elymus … JF2001    107.    195.            354  10.4     1.82 
    ##  8 hermaphroditic Elymus … JF2001     84.3    69.9           266   4.29    0.830
    ##  9 hermaphroditic Festuca… JF2001     17.9    35.7           748   1.31    2.00 
    ## 10 hermaphroditic Phleum … CS2021     38.8    43.6            63   5.49    1.12 
    ## # ℹ 13 more rows
    ## # ℹ 3 more variables: max_poll <int>, min_poll <int>, stig_nplants <int>

``` r
write.csv(pollen_load_tab, "processed-data/pollen-loads-table.csv", row.names=F)
```
