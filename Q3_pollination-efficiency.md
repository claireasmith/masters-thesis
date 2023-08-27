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
# Get estimate of pollen production per flower
# first make sure all species have values for anthers per flower
anthsum <- sizenum.dat %>% group_by(Species) %>% 
  summarize(m_anth_flw = mean(Anth_per_flw, na.rm=T))
sizenum.dat <- sizenum.dat %>%
  mutate(Avg_pol_flw = Avg_pol_anth*Anth_per_flw)


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
            SE_prod_anth_sp = Sd_prod_anth_sp/n_prod_anth,
            Avg_prod_flw_sp = mean(Avg_pol_flw, na.rm=T),
            Sd_prod_flw_sp = sd(Avg_pol_flw, na.rm=T),
            n_prod_flw = n(),
            SE_prod_flw_sp = Sd_prod_flw_sp/n_prod_flw)
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
  # filter(Species != "Thalictrum dioicum") %>% 
  ggplot(aes(x=Avg_prod_flw_sp, y=Avg_stig_sp, shape=Sex_sys) ) + 
  geom_point(size=4, alpha=0.7, aes(color=Sex_sys)) + 
  geom_errorbar(aes(ymin = pmax(0,Avg_stig_sp - Sd_stig_sp), ymax = Avg_stig_sp + Sd_stig_sp), color="black") + 
  geom_errorbarh(aes(xmin = Avg_prod_flw_sp - Sd_prod_flw_sp, xmax = Avg_prod_flw_sp + Sd_prod_flw_sp), color="black") + 
  
  # Axes
  scale_x_continuous(name = expression("Pollen production per flower"), na.value=0, expand = c(0, 0)) + 
  scale_y_continuous(name = "Pollen capture per flower", na.value = 0, expand = c(0, 0.2)) + 
  guides(shape=guide_legend(title="Sex system"),
         color=guide_legend(title="Sex system")) + 
  scale_shape_manual(values = c(15, 16, 17),
                     labels=c("Dioecious", "Monoecious", "Hermaphroditic")) + 
  scale_color_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                     values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  # Guide lines
  geom_segment(aes(x=0, xend = 90000+2000, y=0, yend = 0.0001*90000), colour="black", linetype = "dashed") + 
  geom_segment(aes(x=0, xend = 90000+2000, y=0, yend = 0.001*90000), colour="black", linetype = "dashed") + 
  geom_segment(aes(x=0, xend = 20000, y=0, yend = 0.01*20000), colour="black", linetype = "dashed") + 
  geom_segment(aes(x=0, xend = 2000, y=0, yend = 0.1*2000), colour="black", linetype = "dashed") +
  geom_text(x=90000+5000, y=0.0001*90000, label = "0.01%", colour="black") + 
  geom_text(x=90000+4000, y=0.001*90000+5, label = "0.1%", colour="black") +
  geom_text(x=20000+1000, y=0.01*20000+5, label = "1%", colour="black") +
  geom_text(x=2000+400, y=0.1*2000+5, label = "10%", colour="black") +
  
  theme_cs(font="sans", fontsize=18) + 
  theme(legend.position = c(0.85,0.85))
```

![](Q3_pollination-efficiency_files/figure-gfm/plot%20efficiency-1.png)<!-- -->

``` r
## Plot pollen capture vs pollen production:
pdat %>% 
  filter(Species != "Thalictrum dioicum") %>% 
  ggplot(aes(x=Avg_prod_flw_sp, y=Avg_stig_sp, shape=Sex_sys) ) + 
  geom_point(size=4, alpha=0.7, aes(color=Sex_sys)) + 
  geom_errorbar(aes(ymin = pmax(0,Avg_stig_sp - Sd_stig_sp), ymax = Avg_stig_sp + Sd_stig_sp), color="black") + 
  geom_errorbarh(aes(xmin = Avg_prod_flw_sp - Sd_prod_flw_sp, xmax = Avg_prod_flw_sp + Sd_prod_flw_sp), color="black") + 
  
  # Axes
  scale_x_continuous(name = expression("Pollen production per flower"), na.value=0, expand = c(0, 0)) + 
  scale_y_continuous(name = "Pollen capture per flower", na.value = 0, expand = c(0, 0.2)) + 
  guides(shape=guide_legend(title="Sex system"),
         color=guide_legend(title="Sex system")) + 
  scale_shape_manual(values = c(15, 16, 17),
                     labels=c("Dioecious", "Monoecious", "Hermaphroditic")) + 
  scale_color_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                     values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  # Guide lines
  geom_segment(aes(x=0, xend = 20000+2000, y=0, yend = 0.0001*20000), colour="black", linetype = "dashed") + 
  geom_segment(aes(x=0, xend = 20000+2000, y=0, yend = 0.001*20000), colour="black", linetype = "dashed") + 
  geom_segment(aes(x=0, xend = 15000, y=0, yend = 0.01*15000), colour="black", linetype = "dashed") + 
  geom_segment(aes(x=0, xend = 2000, y=0, yend = 0.1*2000), colour="black", linetype = "dashed") +
  geom_text(x=20000+4000, y=0.0001*20000+2, label = "0.01%", colour="black") + 
  geom_text(x=20000+3000, y=0.001*20000+5, label = "0.1%", colour="black") +
  geom_text(x=15000+1000, y=0.01*15000+5, label = "1%", colour="black") +
  geom_text(x=2000+400, y=0.1*2000+5, label = "10%", colour="black") +
  
  theme_cs(font="sans", fontsize=18) + 
    theme(legend.position = c(0.85,0.7))
```

![](Q3_pollination-efficiency_files/figure-gfm/plot%20efficiency%20wo%20T%20dioicum-1.png)<!-- -->

``` r
## A plot with labels, for reference (zoom into saved file to see labels):
pdat %>% 
  ggplot(aes(x=Avg_prod_flw_sp, y=Avg_stig_sp, shape=Sex_sys, color=Sex_sys, label=Species) ) + 
  geom_point(size=3, alpha=0.8) + 
  geom_errorbar(aes(ymin = pmax(0,Avg_stig_sp - Sd_stig_sp), ymax = Avg_stig_sp + Sd_stig_sp), color="black") + 
  geom_errorbarh(aes(xmin = Avg_prod_flw_sp - Sd_prod_flw_sp, xmax = Avg_prod_flw_sp + Sd_prod_flw_sp), color="black") + 
  geom_label(size=2) + 
  # Axes
  scale_x_continuous(name = expression("Pollen production per flower"), na.value=0) + 
  scale_y_continuous(name = "Pollen capture per flower", na.value = 0) + 
  guides(shape=guide_legend(title="Sex system"), color=guide_legend(title="Sex system")) + 
  scale_shape_manual(values = c(15, 16, 17),
                     labels=c("Dioecious", "Monoecious", "Hermaphroditic")) + 
  scale_color_manual(labels=c("Dioecious", "Monoecious", "Hermaphroditic"),
                     values=c("#66C2A5","#8DA0CB","#FC8D62")) + 
  # geom_text(aes(label=Species), size=2, color="black") + 
  theme_cs(font="sans", fontsize=20)
```

![](Q3_pollination-efficiency_files/figure-gfm/plot%20efficiency%20w%20labels-1.png)<!-- -->

``` r
# Create table of pollination efficiency

eff_tab <- pdat %>% dplyr::group_by(Sex_sys, Species) %>% 
  dplyr::summarize(Avg_stig_sp = Avg_stig_sp,
            Sd_stig_sp = Sd_stig_sp,
            Avg_prod_flw_sp = Avg_prod_flw_sp,
            Sd_prod_flw_sp = Sd_prod_flw_sp,
            polleff = Avg_stig_sp/Avg_prod_flw_sp)
```

    ## `summarise()` has grouped output by 'Sex_sys'. You can override using the
    ## `.groups` argument.

``` r
print(eff_tab, width = 90)
```

    ## # A tibble: 18 × 7
    ## # Groups:   Sex_sys [3]
    ##    Sex_sys        Species                 Avg_stig_sp Sd_stig_sp Avg_prod_flw_sp
    ##    <fct>          <ord>                         <dbl>      <dbl>           <dbl>
    ##  1 dioecious      Rumex acetosella               2.26       2.43          12423.
    ##  2 dioecious      Thalictrum dioicum             3.69       6.48          79162.
    ##  3 monoecious     Amaranthus retroflexus        22.4       26.9           13242.
    ##  4 monoecious     Ambrosia artemisiifolia       51.0       76.9            3707.
    ##  5 monoecious     Carex communis                11.3       16.1            1854.
    ##  6 monoecious     Carex hirtifolia              27.3       54.2            1994.
    ##  7 monoecious     Carex pedunculata             10.6       13.2            1574.
    ##  8 monoecious     Carex plantaginea              9.92      22.8            3019.
    ##  9 monoecious     Carex stipata                  4.66      11.6            1246.
    ## 10 monoecious     Rumex crispus                  6.63      12.5            4545.
    ## 11 monoecious     Scirpus microcarpus            5.16       5.69            507.
    ## 12 hermaphroditic Bromus inermis                50.6       75.7           16976.
    ## 13 hermaphroditic Chenopodium album              3.09       3.73           2627.
    ## 14 hermaphroditic Elymus repens                 84.3       69.9           14986.
    ## 15 hermaphroditic Festuca campestris            17.9       35.7           19554.
    ## 16 hermaphroditic Phleum pratense               38.8       43.6            5090.
    ## 17 hermaphroditic Plantago lanceolata           95.0      107.            10995.
    ## 18 hermaphroditic Schizachne purpurascens       92.3      121.              861.
    ##    Sd_prod_flw_sp   polleff
    ##             <dbl>     <dbl>
    ##  1          3386. 0.000182 
    ##  2         34338. 0.0000467
    ##  3          3092. 0.00169  
    ##  4          1955. 0.0138   
    ##  5           265. 0.00609  
    ##  6           456. 0.0137   
    ##  7           378. 0.00675  
    ##  8           345. 0.00328  
    ##  9           220. 0.00374  
    ## 10           997. 0.00146  
    ## 11           191. 0.0102   
    ## 12          7355. 0.00298  
    ## 13           687. 0.00117  
    ## 14          5545. 0.00562  
    ## 15          5841. 0.000915 
    ## 16          1786. 0.00762  
    ## 17          3073. 0.00864  
    ## 18           298. 0.107

``` r
write.csv(eff_tab, "processed-data/efficiency-table.csv", row.names = F)
```
