Overall stigmatic pollen capture
================
Claire Smith
2023-06-26

``` r
# Load packages
library(tidyverse)
library(ggplot2)
library(ggridges)

# Source files
source("clean-dat.R") # get cleaned data
source("theme_cs.R")
```

``` r
# Load data
stig <- read.csv("processed-data/stig-no-rep-spp.csv", stringsAsFactors = T)
```

``` r
# head(stig)
# summary(stig)
# check for NAs
# stig[which(is.na(stig$Flw_pollen)),] # some NAs from when the data was joined - some individuals had 
```

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
# pollen_load_tab %>% filter(stig_nplants>=10)

# pollen_load_tab %>% group_by(Sex_sys) %>% summarize(N=n())
# pollen_load_tab %>% group_by(Species) %>% summarize(N=n()) %>% print(n=Inf)

write.csv(pollen_load_tab, "processed-data/pollen-loads-table.csv", row.names=F)
```

Visualize distributions of pollen loads across species

``` r
stig %>% 
  # Arrange species in order of sex system and alphabetically by species -- edit: just arrange in order of species alphabetically to match tables for easy comparison
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>%
  # arrange(desc(Sex_sys), desc(Species)) %>% 
  arrange(desc(Species)) %>% 
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
```

    ## Picking joint bandwidth of 0.358

    ## Warning: Removed 19 rows containing non-finite values
    ## (`stat_density_ridges()`).

![](Q2-1_pollen-capture-summary_files/figure-gfm/stigma%20ridgeplots-1.png)<!-- -->

Make this more readable - add sex system to label

``` r
stig %>% 
  mutate(SS = case_when(Sex_sys == "hermaphroditic" ~ "H",
                        Sex_sys == "monoecious" ~ "M",
                        Sex_sys == "dioecious" ~ "D"),
         Species_SS = as.factor(paste0(Species, " (", SS, ")"))) %>% 
  
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic"),
                          ordered=T) ) %>% 
  # arrange(desc(Sex_sys), desc(Species_SS)) %>% 

  arrange(desc(Species_SS)) %>% 

  mutate(Species_SS =  factor(Species_SS, levels = unique(Species_SS), ordered = T)) %>% 

  ggplot(aes(x=log(Flw_pollen+1), y=Species_SS, fill=Sex_sys)) +
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
```

    ## Picking joint bandwidth of 0.358

    ## Warning: Removed 19 rows containing non-finite values
    ## (`stat_density_ridges()`).

![](Q2-1_pollen-capture-summary_files/figure-gfm/stigma%20ss%20labeled%20ridgeplots-1.png)<!-- -->

``` r
stig %>% 
    # Arrange species in order of sex system and alphabetically by species
  mutate(Sex_sys = as.character(Sex_sys),
         Sex_sys = factor(Sex_sys, levels=c("dioecious", "monoecious", "hermaphroditic")) ) %>% 
  # arrange(Sex_sys, Species) %>% 
  arrange(Species) %>% 
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

![](Q2-1_pollen-capture-summary_files/figure-gfm/stigma%20box%20plots-1.png)<!-- -->

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

![](Q2-1_pollen-capture-summary_files/figure-gfm/sex%20sys%20stigma%20box%20plots-1.png)<!-- -->
