Summary pollen size
================
Claire Smith
2023-06-26

``` r
## Load packages
library(tidyverse)
library(ggplot2)
library(ggridges)
```

    ## Warning: package 'ggridges' was built under R version 4.3.1

``` r
## Source files
# custom
source("theme_cs.R")
# data cleaning
source("clean-stig-01.R")
source("clean-stig-02.R")
source("clean-anthers-04.R") # calls 1, 2, 3 within itself
```

``` r
# Load data - size-prod-all.csv has per individual average size
sizefull <- read.csv("processed-data/size-prod-all.csv", stringsAsFactors = T)
# head(sizefull)
# str(sizefull)
# size-CS2021-within-inds.csv has within-individual variation in pollen size as well
sizewithin <- read.csv("processed-data/size-CS2021-within-inds.csv", stringsAsFactors = T)
# head(sizewithin)
# summary(sizewithin)
```

``` r
# Remove NAs
sizefull <- sizefull %>% 
  filter(!is.na(Avg_diam))

# Keep only species with at least 5 individuals
sf_sum <- sizefull %>% group_by(Species) %>% dplyr::summarize(n=n())
keep_vec <- sf_sum$Species[which(sf_sum$n>=5)]
sf_all_5ormore <- sizefull[sizefull$Species %in% keep_vec,]

sizefull <- sf_all_5ormore %>% droplevels()
```

``` r
size_ridges <- sizefull %>% 
  # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(desc(Sex_sys), desc(Species)) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot(aes(x=Avg_diam, y=Species, fill=Sex_sys)) +
  geom_density_ridges2() + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#FC8D62"),
                    labels=c("Dioecious", "Monoecious", "Hermaphroditic")) +
  scale_x_continuous(name = expression("Pollen diameter ("*mu*"m)")) + 

  scale_y_discrete() +
  ylab(label="") + 
  guides(fill="none") + 
  theme_cs(font="sans", fontsize=15) + 
  theme(axis.text.y = element_text(face = "italic"),
        axis.ticks.y = element_line())

size_ridges
```

    ## Picking joint bandwidth of 0.512

![](anal-size-summary_files/figure-gfm/plot%20size%20distributions-1.png)<!-- -->

``` r
size_box <- sizefull %>% 
    # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(Sex_sys, Species) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot() + 
  geom_boxplot(aes(x=Species, y=Avg_diam, fill=Sex_sys)) + 

  scale_y_continuous(name = expression("Pollen diameter ("*mu*"m)")) + 
  scale_x_discrete(name = "") + 
  
  guides(fill="none") + 
  scale_fill_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                    values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  theme_cs(fontsize=15, font="sans") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic"))

size_box
```

![](anal-size-summary_files/figure-gfm/size%20box%20plots-1.png)<!-- -->

``` r
size_box_ss <- sizefull %>% 
  # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(Sex_sys, Species) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot() + 
  geom_boxplot(aes(x=Sex_sys, y=Avg_diam, fill=Sex_sys),
               width=0.4) + 
  
  scale_y_continuous(name = expression("Pollen diameter ("*mu*"m)")) + 
  scale_x_discrete(name = "", labels = NULL) + 
  
  guides(fill=guide_legend(title="Sex system")) + 
  scale_fill_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                    values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  theme_cs(fontsize=15, font="sans") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "right")

size_box_ss
```

![](anal-size-summary_files/figure-gfm/sex%20sys%20size%20box%20plots-1.png)<!-- -->

``` r
size_sum <- sizefull %>% group_by(Sex_sys, Species) %>% 
  summarize(df = N_diam-1,
            s2 = Sd_diam^2,
            ss = s2 * df,
            m = Avg_diam,
            nn = N_diam) %>% 
  group_by(Sex_sys, Species) %>% 
  summarize(Mean_diam = sum(m*nn)/sum(nn),
            Sd_diam_pooled = sum(ss)/sum(df),
            Cv = Sd_diam_pooled/Mean_diam,
            Max_diam = max(m),
            Min_diam = min(m),
            N_grains = sum(nn),
            N_ind = n()) %>% 
  filter(N_ind>=5) # Keep only species with at least 5 individuals
```

    ## Warning: Returning more (or less) than 1 row per `summarise()` group was deprecated in
    ## dplyr 1.1.0.
    ## ℹ Please use `reframe()` instead.
    ## ℹ When switching from `summarise()` to `reframe()`, remember that `reframe()`
    ##   always returns an ungrouped data frame and adjust accordingly.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

    ## `summarise()` has grouped output by 'Sex_sys', 'Species'. You can override
    ## using the `.groups` argument.
    ## `summarise()` has grouped output by 'Sex_sys'. You can override using the
    ## `.groups` argument.

``` r
print(size_sum, n=27)
```

    ## # A tibble: 27 × 9
    ## # Groups:   Sex_sys [3]
    ##    Sex_sys    Species Mean_diam Sd_diam_pooled     Cv Max_diam Min_diam N_grains
    ##    <fct>      <fct>       <dbl>          <dbl>  <dbl>    <dbl>    <dbl>    <dbl>
    ##  1 dioecious  Rumex …      19.9           3.44 0.173      22.4     15.0  316076 
    ##  2 dioecious  Thalic…      16.5           2.65 0.160      17.7     14.8 1824160 
    ##  3 hermaphro… Agropy…      37.6           4.69 0.125      39.2     34.9  161853 
    ##  4 hermaphro… Avenul…      28.6           6.25 0.219      29.7     25.9   25787 
    ##  5 hermaphro… Bromus…      32.9          13.2  0.400      40.4     30.1  508033 
    ##  6 hermaphro… Chenop…      24.2           2.01 0.0832     29.0     23.3   40246 
    ##  7 hermaphro… Elymus…      35.4           7.02 0.198      38.2     31.9  778882 
    ##  8 hermaphro… Elymus…      37.2           6.65 0.179      39.2     33.8  481773 
    ##  9 hermaphro… Festuc…      28.7           5.01 0.175      32.6     23.9 1153669 
    ## 10 hermaphro… Festuc…      29.4           6.53 0.222      31.1     27.8  103862 
    ## 11 hermaphro… Festuc…      28.7           6.19 0.216      30.5     26.6   77397.
    ## 12 hermaphro… Hieroc…      22.9          26.0  1.13       25.6     18.1  220134.
    ## 13 hermaphro… Koeler…      21.8           3.22 0.148      24.5     20.7  181789 
    ## 14 hermaphro… Phalar…      32.0           6.70 0.209      32.5     31.1   28476.
    ## 15 hermaphro… Phleum…      30.0           5.52 0.184      33.5     26.8  167430.
    ## 16 hermaphro… Planta…      21.1           5.58 0.264      27.3     18.6  265949 
    ## 17 hermaphro… Poa ju…      21.5          10.2  0.474      23.0     20.0   39470 
    ## 18 hermaphro… Schiza…      29.4           3.11 0.106      30.3     28.6   12436 
    ## 19 hermaphro… Stipa …      28.6           4.16 0.145      30.1     26.5   24386.
    ## 20 monoecious Amaran…      24.0           2.78 0.115      28.1     21.9     166 
    ## 21 monoecious Ambros…      17.9           2.73 0.152      19.0     17.0  280272 
    ## 22 monoecious Carex …      24.1           9.99 0.414      26.1     18.6  327616 
    ## 23 monoecious Carex …      25.2           1.93 0.0766     26.1     24.7  372071 
    ## 24 monoecious Carex …      24.6           2.29 0.0931     26.0     23.2  507422 
    ## 25 monoecious Carex …      24.5           2.22 0.0908     26.7     23.8  147145 
    ## 26 monoecious Rumex …      24.0           2.57 0.107      24.9     23.0   87014 
    ## 27 monoecious Scirpu…      22.9           4.04 0.177      23.4     22.0   77145 
    ## # ℹ 1 more variable: N_ind <int>

``` r
# write it to a file
write.csv(size_sum, "processed-data/pollen-size-table.csv", row.names=F)
```

``` r
size_sum_ss <- size_sum %>% 
  group_by(Sex_sys) %>% 
  summarize(msize = mean(Mean_diam),
            sdsize = sd(Mean_diam),
            maxsize = max(Max_diam),
            minsize = min(Min_diam),
            nspp=n(),
            nplantstot = sum(N_ind))
size_sum_ss
```

    ## # A tibble: 3 × 7
    ##   Sex_sys        msize sdsize maxsize minsize  nspp nplantstot
    ##   <fct>          <dbl>  <dbl>   <dbl>   <dbl> <int>      <int>
    ## 1 dioecious       18.2   2.37    22.4    14.8     2         60
    ## 2 hermaphroditic  28.8   5.23    40.4    18.1    17        429
    ## 3 monoecious      23.4   2.32    28.1    17.0     8        172

``` r
# write it to a file
write.csv(size_sum_ss, "processed-data/pollen-size-table-ss.csv", row.names=F)
```

write.csv(size_sum, “Data/clean/size_summary.csv”, row.names = F)
write.csv(size_sum_ss, “Tables/size_summary_sexsys_table.csv”,
row.names=F)
