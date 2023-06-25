Pollination efficiency
================
Claire Smith
2023-06-23

Calculating pollen transfer efficiency - percentage of grains captured
per stigma to grains produced per anther.

Do plants that produce more pollen also capture more pollen? And does
this relationship change with sex system (do dioecious species have
lower efficiencies because they cannot capture self pollen?)

``` r
## Load packages
library(tidyverse)
library(ggplot2)

## Source files
# custom
source("theme_cs.R")
# data cleaning
source("clean-stig-01.R")
source("clean-stig-02.R")
source("clean-anthers-04.R") # calls 1, 2, 3 within itself
```

``` r
# Stigmatic pollen capture
pc.dat <- read.csv("processed-data/stig-no-rep-spp.csv", stringsAsFactors = T) 
# head(pc.dat)
# str(pc.dat)
#Fix a typo in stigma data
pc.dat$Species <- gsub("Schizacne purpurascens", "Schizachne purpurascens", pc.dat$Species)
# Pollen production (and size)
sizenum.dat <- read.csv("processed-data/size-prod-all.csv", stringsAsFactors = T)
# head(sizenum.dat)
# str(sizenum.dat)
```

``` r
# The stigma and anther data is at the per-individual level right now - summarize it so that it's at the per-species level
# # ** go back and do a satterthwaite pooled sd instead of just lumping it all together - for now just combine all by species
# stig_sum_ind <- pc.dat %>% 
#   group_by(Species, Date, Site, Plant) %>% 
#   summarize(Avg_pollen_ind = mean(Flw_pollen, na.rm=T), # within-individual avg pollen receipt
#             Sd_pollen_ind = sd(Flw_pollen, na.rm=T), # within-individual sd pollen receipt
#             N_flw = n()) %>% 
#   # Now ungroup to find total number of individuals sampled in each sp
#   group_by(Species) %>% 
#   mutate(N_ind=n()) %>% 
#   ungroup()
# 
# # Pool together to get per-species mean pollen receipt and sd pollen receipt
# stig_sum_spp <- stig_sum_ind %>% 
#   group_by(Species) %>% 
#   summarise(Avg_pollen_sp = mean(Avg_pollen_ind, na.rm=T),
#             sd_num = sum((N_flw-1)*Sd_pollen_ind^2),
#             sd_denom = sum(N_flw) - N_ind,
#             Sd_pollen_sp_pooled = sd_num/sd_denom)

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
# Join data for plotting - use inner_join() because not all species have both capture and prod data
# Keep JF data for pollen production when species overlap - I trust the automated counting with the
# particle counter more than my subsampling counting. It counts all the particles that go through it
# vs I only count a fraction of the total pollen and assume it's uniformly distributed within the sample.
prod_spp$Species[which(prod_spp$source == "CS2021")]
```

    ## [1] Amaranthus retroflexus  Ambrosia artemisiifolia Chenopodium album      
    ## [4] Plantago lanceolata     Rumex acetosella        Thalictrum dioicum     
    ## 36 Levels: Acristata Adasyd Agropyron trachycaulum ... Thalictrum dioicum

``` r
prod_spp_filt <- prod_spp %>% filter(source != "CS2021" |  Species == "Amaranthus retroflexus")

pdat <- inner_join(stig_spp, prod_spp_filt, by = c("Species", "Sex_sys"))
anti_pdat <- anti_join(stig_spp, prod_spp_filt, by = c("Species"))
anti_pdat2 <- anti_join(prod_spp_filt, stig_spp, by = c("Species"))
# head(pdat)
# summary(pdat)
# View(pdat)
```

2 dioecious species, 8 monoecious species, 9 hermaphroditic species

``` r
# Re-order levels in sex system column 
pdat$Sex_sys <- as.character(pdat$Sex_sys)
pdat$Sex_sys <- factor(pdat$Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic"))

# Arrange data so that species are grouped by sex system and ordered alphabetically
pdat <- arrange(pdat, Sex_sys, Species)
pdat$Species <- factor(pdat$Species, levels = unique(pdat$Species), ordered = T)
# names(pdat)
# [1] "Species"          "Sex_sys"          "Avg_stig_sp"      "Sd_stig_sp"       "n_stig"           "SE_stig_sp"      
#  [7] "Avg_prod_anth_sp" "Sd_prod_anth_sp"  "n_prod_anth"      "SE_prod_anth_sp" 
```

Plots

``` r
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

![](anal_efficiency_files/figure-gfm/plot%20efficiency-1.png)<!-- -->

Tables

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

    ## # A tibble: 19 × 7
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
    ## 14 hermaphroditic Elymus innovatus             107.       195.  
    ## 15 hermaphroditic Elymus repens                 84.3       69.9 
    ## 16 hermaphroditic Festuca campestris            17.9       35.7 
    ## 17 hermaphroditic Phleum pratense               38.8       43.6 
    ## 18 hermaphroditic Plantago lanceolata           87.0      119.  
    ## 19 hermaphroditic Schizachne purpurascens       92.3      121.  
    ##    Avg_prod_anth_sp Sd_prod_anth_sp polleff
    ##               <dbl>           <dbl>   <dbl>
    ##  1            2070.           564.    0.109
    ##  2            3320.          1295.    0.111
    ##  3            2648.           618.    0.846
    ##  4             766.           388.    6.66 
    ##  5             607.           100.    1.86 
    ##  6             665.           152.    4.11 
    ##  7             518.           126.    2.05 
    ##  8            1006.           115.    0.985
    ##  9             413.            72.9   1.13 
    ## 10             758.           166.    0.875
    ## 11             169.            63.7   3.05 
    ## 12            5659.          2452.    0.894
    ## 13             525.           137.    0.587
    ## 14            7802.          2754.    1.38 
    ## 15            4995.          1848.    1.69 
    ## 16            6518.          1947.    0.274
    ## 17            2163.          2663.    1.79 
    ## 18            2749.           768.    3.17 
    ## 19             287.            99.4  32.1

``` r
write.csv(eff_tab, "processed-data/efficiency-table.csv", row.names = F)
```
