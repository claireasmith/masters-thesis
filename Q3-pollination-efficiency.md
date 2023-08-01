Pollination efficiency
================
Claire Smith
2023-06-23

### Question

Do plants that produce more pollen also capture more pollen?

### Goal

Calculating and plotting pollen transfer efficiency (PTE) for 19
wind-pollinated species. PTE = percentage of grains captured per stigma
over grains produced per anther.

### Loading and preparing data

``` r
## Load packages
library(tidyverse)
library(ggplot2)

## Source files
# custom
source("theme_cs.R")
# data cleaning
source("clean-dat.R")
```

``` r
## Stigmatic pollen capture
pc.dat <- read.csv("processed-data/stig-no-rep-spp.csv", stringsAsFactors = T) 
# head(pc.dat)
# str(pc.dat)
#Fix a typo in stigma data
pc.dat$Species <- gsub("Schizacne purpurascens", "Schizachne purpurascens", pc.dat$Species)

## Pollen production (and size)
sizenum.dat <- read.csv("processed-data/size-prod-norep.csv", stringsAsFactors = T)
# head(sizenum.dat)
# str(sizenum.dat)
```

``` r
# The stigma and anther data is at the per-individual level right now - summarize it so that it's at the per-species level
stig_spp <- pc.dat %>% 
  group_by(Species, Sex_sys, source) %>%
  summarise(Avg_stig_sp = mean(Flw_pollen, na.rm=T), 
            Sd_stig_sp = sd(Flw_pollen, na.rm=T),
            n_stig = n(),
            SE_stig_sp = Sd_stig_sp/sqrt(n_stig))
```

    ## `summarise()` has grouped output by 'Species', 'Sex_sys'. You can override
    ## using the `.groups` argument.

``` r
# View(stig_spp)
prod_spp <- sizenum.dat %>% 
  group_by(Species, Sex_sys, source) %>% 
  summarize(Avg_prod_anth_sp = mean(Avg_pol_anth, na.rm=T),
            Sd_prod_anth_sp = sd(Avg_pol_anth, na.rm=T),
            n_prod_anth = n(),
            SE_prod_anth_sp = Sd_prod_anth_sp/n_prod_anth)
```

    ## `summarise()` has grouped output by 'Species', 'Sex_sys'. You can override
    ## using the `.groups` argument.

``` r
# View(prod_spp)

## Join data
# I'll use inner_join() because not all species have both capture and prod data.
pdat <- inner_join(stig_spp, prod_spp, by = c("Species", "Sex_sys"))
# which species didn't make it from the stigma data?
anti_pdat <- anti_join(stig_spp, prod_spp, by = c("Species"))
# which species didn't make it from the prod data?
anti_pdat2 <- anti_join(prod_spp, stig_spp, by = c("Species"))
# head(pdat)
# summary(pdat)
# View(pdat)

# 19 total species: 2 dioecious species, 8 monoecious species, 9 hermaphroditic species. 

# Re-order levels in sex system column 
pdat$Sex_sys <- as.character(pdat$Sex_sys)
pdat$Sex_sys <- factor(pdat$Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic"))

# Arrange data so that species are grouped by sex system and ordered alphabetically
pdat <- arrange(pdat, Sex_sys, Species)
pdat$Species <- factor(pdat$Species, levels = unique(pdat$Species), ordered = T)
```

### Plots and tables

``` r
## Plot pollen capture vs pollen production:
pdat %>% 
  ggplot(aes(x=Avg_prod_anth_sp, y=Avg_stig_sp, shape=Sex_sys, color=Sex_sys) ) + 
  geom_point(size=5, alpha=0.8) + 
  geom_errorbar(aes(ymin = pmax(0,Avg_stig_sp - Sd_stig_sp), ymax = Avg_stig_sp + Sd_stig_sp), color="black") + 
  geom_errorbarh(aes(xmin = Avg_prod_anth_sp - Sd_prod_anth_sp, xmax = Avg_prod_anth_sp + Sd_prod_anth_sp), color="black") + 
  scale_x_continuous(name = expression("Pollen production per anther"), na.value=0) + 
  scale_y_continuous(name = "Pollen capture per flower", na.value = 0) + 
  guides(shape=guide_legend(title="Sex system"), color=guide_legend(title="Sex system")) + 
  scale_shape_manual(values = c(15, 16, 17),
                     labels=c("Dioecious", "Monoecious", "Hermaphroditic")) + 
  scale_color_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                     values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  # geom_text(aes(label=Species), size=2, color="black") + 
  theme_cs(font="sans", fontsize=20)
```

![](Q3-pollination-efficiency_files/figure-gfm/plot%20efficiency-1.png)<!-- -->

``` r
## A plot with labels, for reference (zoom into saved file to see labels):
pdat %>% 
  ggplot(aes(x=Avg_prod_anth_sp, y=Avg_stig_sp, shape=Sex_sys, color=Sex_sys, label=Species) ) + 
  geom_point(size=3, alpha=0.8) + 
  geom_errorbar(aes(ymin = pmax(0,Avg_stig_sp - Sd_stig_sp), ymax = Avg_stig_sp + Sd_stig_sp), color="black") + 
  geom_errorbarh(aes(xmin = Avg_prod_anth_sp - Sd_prod_anth_sp, xmax = Avg_prod_anth_sp + Sd_prod_anth_sp), color="black") + 
  geom_label(size=2) + 
  
  scale_x_continuous(name = expression("Pollen production per anther"), na.value=0) + 
  scale_y_continuous(name = "Pollen capture per flower", na.value = 0) + 
  guides(shape=guide_legend(title="Sex system"), color=guide_legend(title="Sex system")) + 
  scale_shape_manual(values = c(15, 16, 17),
                     labels=c("Dioecious", "Monoecious", "Hermaphroditic")) + 
  scale_color_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                     values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  # geom_text(aes(label=Species), size=2, color="black") + 
  theme_cs(font="sans", fontsize=20)
```

![](Q3-pollination-efficiency_files/figure-gfm/plot%20efficiency%20w%20labels-1.png)<!-- -->

``` r
# Create table of pollination efficiency

eff_tab <- pdat %>% dplyr::group_by(Sex_sys, Species) %>% 
  dplyr::summarize(Avg_stig_sp = Avg_stig_sp,
            Sd_stig_sp = Sd_stig_sp,
            Avg_prod_anth_sp = Avg_prod_anth_sp,
            Sd_prod_anth_sp = Sd_prod_anth_sp,
            polleff = 100 * Avg_stig_sp/Avg_prod_anth_sp)
```

    ## `summarise()` has grouped output by 'Sex_sys'. You can override using the
    ## `.groups` argument.

``` r
print(eff_tab, width = 90)
```

    ## # A tibble: 18 × 7
    ## # Groups:   Sex_sys [3]
    ##    Sex_sys        Species                 Avg_stig_sp Sd_stig_sp
    ##    <fct>          <ord>                         <dbl>      <dbl>
    ##  1 dioecious      Rumex acetosella               2.26       2.43
    ##  2 dioecious      Thalictrum dioicum             3.69       6.48
    ##  3 monoecious     Amaranthus retroflexus        22.4       26.9 
    ##  4 monoecious     Ambrosia artemisiifolia       51.0       76.9 
    ##  5 monoecious     Carex communis                11.3       16.1 
    ##  6 monoecious     Carex hirtifolia              27.3       54.2 
    ##  7 monoecious     Carex pedunculata             10.6       13.2 
    ##  8 monoecious     Carex plantaginea              9.92      22.8 
    ##  9 monoecious     Carex stipata                  4.66      11.6 
    ## 10 monoecious     Rumex crispus                  6.63      12.5 
    ## 11 monoecious     Scirpus microcarpus            5.16       5.69
    ## 12 hermaphroditic Bromus inermis                50.6       75.7 
    ## 13 hermaphroditic Chenopodium album              3.09       3.73
    ## 14 hermaphroditic Elymus repens                 84.3       69.9 
    ## 15 hermaphroditic Festuca campestris            17.9       35.7 
    ## 16 hermaphroditic Phleum pratense               38.8       43.6 
    ## 17 hermaphroditic Plantago lanceolata           95.0      107.  
    ## 18 hermaphroditic Schizachne purpurascens       92.3      121.  
    ##    Avg_prod_anth_sp Sd_prod_anth_sp polleff
    ##               <dbl>           <dbl>   <dbl>
    ##  1            2070.           564.    0.109
    ##  2            3320.          1295.    0.111
    ##  3            2648.           618.    0.846
    ##  4             741.           391.    6.88 
    ##  5             618.            88.5   1.83 
    ##  6             665.           152.    4.11 
    ##  7             525.           126.    2.03 
    ##  8            1006.           115.    0.985
    ##  9             415.            73.2   1.12 
    ## 10             758.           166.    0.875
    ## 11             169.            63.7   3.05 
    ## 12            5659.          2452.    0.894
    ## 13             525.           137.    0.587
    ## 14            4995.          1848.    1.69 
    ## 15            6518.          1947.    0.274
    ## 16            1697.           595.    2.29 
    ## 17            2749.           768.    3.46 
    ## 18             287.            99.4  32.1

``` r
write.csv(eff_tab, "processed-data/efficiency-table.csv", row.names = F)
```
