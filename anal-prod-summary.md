Pollen production summary
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
source("clean-anthers-04.R") # calls 1, 2, 3 within itself
```

``` r
# Pollen production (and size)
sizenum.dat <- read.csv("processed-data/size-prod-all.csv", stringsAsFactors = T)
# head(sizenum.dat)
# str(sizenum.dat)
```

``` r
prod_spp <- sizenum.dat
# Some species were collected by both me and Jannice. She used automated pollen counting 
# with a particle counter. It counts all particles that go through it vs small subsample I 
# counted in mine, makes less assumptions about how uniformly pollen is suspended in the sample -- I'll keep her data where there's overlap. 
# A very long line of code to find out if there are any species in my production data NOT also included in the JF data: 
(unique(prod_spp$Species[which(prod_spp$source == "CS2021")]))[!(unique(prod_spp$Species[which(prod_spp$source == "CS2021")])) %in% unique(prod_spp$Species[which(prod_spp$source == "JF2001" | prod_spp$source == "JF2004" )])]
```

    ## [1] Amaranthus retroflexus
    ## 36 Levels: Acristata Adasyd Agropyron trachycaulum ... Thalictrum dioicum

``` r
# Amaranthus retroflexus is the only one
# Exclude all pollen production data from CS2021 except Amaranthus retroflexus
prod_spp_filt <- prod_spp %>% filter(source != "CS2021" |  Species == "Amaranthus retroflexus")

# Re-order levels in sex system column 
prod_spp_filt$Sex_sys <- as.character(prod_spp_filt$Sex_sys)
prod_spp_filt$Sex_sys <- factor(prod_spp_filt$Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic"))

# Arrange data so that species are grouped by sex system and ordered alphabetically
prod_spp_filt <- arrange(prod_spp_filt, Sex_sys, Species)
prod_spp_filt$Species <- factor(prod_spp_filt$Species, levels = unique(prod_spp_filt$Species), ordered = T)

# Remove NAs
prod_spp_filt <- prod_spp_filt %>% 
  filter(!is.na(Avg_pol_anth))

# Keep only species with at least 5 individuals
p_sum <- prod_spp_filt %>% group_by(Species) %>% dplyr::summarize(n=n())
keep_vec <- p_sum$Species[which(p_sum$n>=5)]
p_all_5ormore <- prod_spp_filt[prod_spp_filt$Species %in% keep_vec,]

prodfull <- p_all_5ormore %>% droplevels()

# After looking at ridgeplot distributions of pollen production it looks like there are some clear outliers in some species - like Phleum pratense
# prodfull %>% filter(Species == "Phleum pratense")
# looks like Ind C11 is the one, with 16161 grains/anther! The rest are all <10,000. I'll
# remove it
prodfull_filt <- prodfull %>% filter(!(Species=="Phleum pratense" & Ind == "C11"))
# prodfull_filt %>% filter(Species == "Phleum pratense") # make sure it's gone
```

``` r
prod_sum <- prodfull_filt %>% 
  group_by(Species, Sex_sys, source) %>% 
  summarize(Avg_prod_anth_sp = mean(Avg_pol_anth, na.rm=T),
            Sd_prod_anth_sp = sd(Avg_pol_anth, na.rm=T),
            n_prod_anth = n(),
            SE_prod_anth_sp = Sd_prod_anth_sp/n_prod_anth) %>% 
  filter(n_prod_anth>=5) %>% 
  droplevels()
```

    ## `summarise()` has grouped output by 'Species', 'Sex_sys'. You can override
    ## using the `.groups` argument.

``` r
print(prod_sum, n=27)
```

    ## # A tibble: 27 × 7
    ## # Groups:   Species, Sex_sys [27]
    ##    Species           Sex_sys source Avg_prod_anth_sp Sd_prod_anth_sp n_prod_anth
    ##    <ord>             <fct>   <fct>             <dbl>           <dbl>       <int>
    ##  1 Rumex acetosella  dioeci… JF2004            2070.           564.           21
    ##  2 Thalictrum dioic… dioeci… JF2004            3320.          1295.           23
    ##  3 Amaranthus retro… monoec… CS2021            2648.           618.           16
    ##  4 Ambrosia artemis… monoec… JF2004             766.           388.           21
    ##  5 Carex communis    monoec… JF2004             607.           100.           18
    ##  6 Carex hirtifolia  monoec… JF2004             665.           152.           18
    ##  7 Carex pedunculata monoec… JF2004             518.           126.           19
    ##  8 Carex stipata     monoec… JF2004             413.            72.9          36
    ##  9 Rumex crispus     monoec… JF2004             758.           166.           20
    ## 10 Scirpus microcar… monoec… JF2004             169.            63.7          14
    ## 11 Agropyron trachy… hermap… JF2001            5395.           927.           10
    ## 12 Avenula hookeri   hermap… JF2001            1719.          1078.            5
    ## 13 Bromus inermis    hermap… JF2001            5659.          2452.           30
    ## 14 Chenopodium album hermap… JF2004             525.           137.           22
    ## 15 Elymus innovatus  hermap… JF2001            7802.          2754.           35
    ## 16 Elymus repens     hermap… JF2001            4995.          1848.           33
    ## 17 Festuca campestr… hermap… JF2001            6518.          1947.           60
    ## 18 Festuca pratensis hermap… JF2001            2473.           931.           14
    ## 19 Festuca rubra     hermap… JF2001            3037.          2565.            8
    ## 20 Hierochloe odora… hermap… JF2001            1382.           900.           67
    ## 21 Koeleria cristata hermap… JF2001            2374.          1045.           25
    ## 22 Phalaris arundin… hermap… JF2001            1496.           493.            7
    ## 23 Phleum pratense   hermap… JF2001            1697.           595.           30
    ## 24 Plantago lanceol… hermap… JF2004            2749.           768.           24
    ## 25 Poa juncifolia    hermap… JF2001            2426.          1037.            6
    ## 26 Schizachne purpu… hermap… JF2004             287.            99.4          12
    ## 27 Stipa columbiana  hermap… JF2001             673.           387.           10
    ## # ℹ 1 more variable: SE_prod_anth_sp <dbl>

``` r
prod_ridges <- prodfull_filt %>% 
  # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(desc(Sex_sys), desc(Species)) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot(aes(x=Avg_pol_anth, y=Species, fill=Sex_sys)) +
  geom_density_ridges2() + 
  scale_fill_manual(values=c("#66C2A5","#8DA0CB","#FC8D62"),
                    labels=c("Dioecious", "Monoecious", "Hermaphroditic")) +
  scale_x_continuous(name = "Pollen per anther") + 

  scale_y_discrete() +
  ylab(label="") + 
  guides(fill="none") + 
  theme_cs(font="sans", fontsize=15) + 
  theme(axis.text.y = element_text(face = "italic"),
        axis.ticks.y = element_line())

prod_ridges
```

    ## Picking joint bandwidth of 377

![](anal-prod-summary_files/figure-gfm/plot%20prod%20distributions-1.png)<!-- -->

``` r
prod_box <- prodfull_filt %>% 
    # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(Sex_sys, Species) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot() + 
  geom_boxplot(aes(x=Species, y=Avg_pol_anth, fill=Sex_sys)) + 

  scale_y_continuous(name = "Pollen per anther") + 
  scale_x_discrete(name = "") + 
  
  guides(fill="none") + 
  scale_fill_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                    values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  theme_cs(fontsize=15, font="sans") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face="italic"))

prod_box
```

![](anal-prod-summary_files/figure-gfm/prod%20boxplots-1.png)<!-- -->

``` r
prod_box_ss <- prodfull_filt %>% 
  # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  arrange(Sex_sys, Species) %>% 
  mutate(Species =  factor(Species, levels = unique(Species), ordered = T)) %>% 
  
  ggplot() + 
  geom_boxplot(aes(x=Sex_sys, y=Avg_pol_anth, fill=Sex_sys),
               width=0.4) + 
  
  scale_y_continuous(name = "Pollen per anther") + 
  scale_x_discrete(name = "", labels = NULL) + 
  
  guides(fill=guide_legend(title="Sex system")) + 
  scale_fill_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                    values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  theme_cs(fontsize=15, font="sans") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "right")

prod_box_ss
```

![](anal-prod-summary_files/figure-gfm/sex%20sys%20prod%20box%20plots-1.png)<!-- -->
