## Combine CS stigmatic pollen load data
## Claire Smith
## 16 June 2023

# Create stig-CS2021-all.csv, a file with all the CS2021 stigma data compiled together

##########################################################################################
# Load libraries
library(tidyverse)

# Load raw stigmatic pollen count data
AR <- read.csv("raw-data/stig-CS2021-Amaranthus-retroflexus.csv")
AA <- read.csv("raw-data/stig-CS2021-Ambrosia-artemisiifolia.csv")
CA <- read.csv("raw-data/stig-CS2021-Chenopodium-album.csv")
DL <- read.csv("raw-data/stig-CS2021-Dicanthelium_linearifolium.csv")
DI <- read.csv("raw-data/stig-CS2021-Dicanthelium-implicatum.csv")
PP <- read.csv("raw-data/stig-CS2021-Phleum-pratense.csv")
PL <- read.csv("raw-data/stig-CS2021-Plantago-lanceolata.csv")
SV <- read.csv("raw-data/stig-CS2021-Setaria-viridis.csv")
RA <- read.csv("raw-data/stig-CS2021-Rumex-acetosella.csv")
TD <- read.csv("raw-data/stig-CS2021-Thalictrum-dioicum.csv")

# Define vector of variables I want in final full dataset
sel_vec <- c("Species", "Date", "Site", "Plant", "Flower", "Flw_pollen", "Stigmas_per_flw", "Flower_height", "Infl_max", "Infl_min", "D1", "D2", "D3", "D4", "D5")

# Ambrosia artenmisiifolia
# glimpse(AA)
# summary(AA)
# View(AA)
# 2 sites: FBF and LP
# 2 dates: 2021-09-14 and 2021-09-20
# max 15 plants per site, 3 flowers per plant
# each flower has two stigmas/ovule but pollen counts done at the per flower level, not per stigma
# hence Stig_pollen will be NA, Flw_pollen will be populated. 
# has Flower_height but not Infl_max or Infl_min
AA_pre <- AA %>% 
  mutate(Flw_pollen = Pollen, Infl_max=NA, Infl_min=NA, source="CS2021") %>% 
  select(all_of(sel_vec))

# Amaranthus retroflexus
# glimpse(AR)
# summary(AR)
# View(AR)
# Two sites: FBF1 and FBF2 - both taken from Fruition Berry Farm two days apart from different subpopulations within the farm
AR_pre <- AR %>% 
  filter(Flower_extra==1) %>% # Some flowers have "extras" - other flowers taken from same height, remove these
  mutate(Flw_pollen = Pollen, Infl_max=NA, Infl_min=NA) %>% 
  select(all_of(sel_vec))
# rbind(AA_pre, AR_pre) # test

# Chenopodium album 
# glimpse(CA)
# summary(CA)
# View(CA)
# 1 stigma per flower 
CA_pre <- CA %>% 
  mutate(Flw_pollen = Pollen, Infl_min=NA, Infl_max=Height_max) %>% 
  select(all_of(sel_vec))
# rbind(AR_pre, CA_pre) # test

# Dicanthelium linearifolium 
# glimpse(DL)
# summary(DL)
# View(DL)
# This species has two stigmas/ovule and they were counted separately - I'll add them together to get flower-level data
DL_pre <- DL %>% group_by(Species, Date, Site, Plant, Flower, Flw_height,
                            D1, D2, D3, D4, D5) %>% #group by everything but the variable I'll sum, Pollen, to keep it all
  summarise(Flw_pollen = sum(Pollen)) %>% 
  mutate(Stigmas_per_flw=2,Infl_max=NA, Infl_min=NA,
         Flower_height=Flw_height) %>% 
  ungroup() %>% 
  select(all_of(sel_vec))
# rbind(CA_pre, DL_pre) # test

# Dichanthelium implicatum 
# glimpse(DI)
# summary(DI)
# View(DI)
# This species also has two stigmas/ovule - sum these together to get per flower pollen counts
# Also has extra info I won't include but should note - flowers per stem, inflorescence max width, flower pos
DI_pre <- DI %>% group_by(Species, Date, Site, Plant, Flower, Flower_height, Infl_max, Infl_min,
                          D1, D2, D3, D4, D5) %>% 
  summarise(Flw_pollen = sum(Pollen)) %>% 
  mutate(Stigmas_per_flw=2) %>% 
  ungroup() %>% 
  select(all_of(sel_vec))
# rbind(DI_pre, DL_pre) # test

# Phleum pratense
# glimpse(PP)
# summary(PP)
# View(PP)
# Also two stigmas per ovule - will sum to get per flower pollen count
# Also has flowers per stem estimate
PP_pre <- PP %>% group_by(Species, Date, Site, Plant, Flower, Flower_height, Infl_max, Infl_min,
                          D1, D2, D3, D4, D5) %>% 
  summarise(Flw_pollen = sum(Pollen)) %>% 
  mutate(Stigmas_per_flw=2) %>% 
  ungroup() %>% 
  select(all_of(sel_vec))
# rbind(DI_pre, PP_pre) # test

# Plantago lanceolata
# View(PL)
# This one has distance to TEN nearest neighbours - i'll still do 5 for consistency. Also have heights of male and female 
# flowers in inflorescence as estimate of sex phase. 
PL_pre <- PL %>% 
  mutate(Flw_pollen = Pollen, Infl_min=Min, Infl_max=Max, Site=Population, Stigmas_per_flw=1, Flower_height=NA) %>% 
  select(all_of(sel_vec))
# rbind(PP_pre, PL_pre) # test

# Setaria viridis
# View(SV)
# Also has count of stems
# 2 stigmas per ovule, stigmas counted separately - i'll add them together
SV_pre <- SV %>% mutate(Site = Population, Infl_max=Flower_max, Infl_min=Flower_min) %>%  
  group_by(Species, Date, Site, Plant, Flower, Flower_height, Infl_max, Infl_min,
                          D1, D2, D3, D4, D5) %>% 
  summarise(Flw_pollen = sum(Pollen)) %>% 
  mutate(Stigmas_per_flw=2) %>% 
  ungroup() %>% 
  select(all_of(sel_vec))
# rbind(PL_pre, SV_pre) # test

# Rumex acetosella
# Also has number of flowers per plant estimate
# View(RA)
# This species has three stigmas/ovule and they were counted separately - I'll add them together to get flower-level data
RA_pre <- RA %>% group_by(Species, Date, Site, Plant, Flower, Flower_height,
                          D1, D2, D3, D4, D5) %>% #group by everything but the variable I'll sum, Pollen, to keep it all
  summarise(Flw_pollen = sum(Pollen, na.rm=T)) %>% 
  mutate(Stigmas_per_flw=3,Infl_max=NA, Infl_min=NA) %>% 
  ungroup() %>% 
  select(all_of(sel_vec))
# rbind(SV_pre, RA_pre) # test 

# Thalictrum dioicum 
# View(TD)
# glimpse(TD)
# has flowering stage estimate 
# infl max is plant height (height of tallest flower)
# I mistakenly called each inflorescence a flower - so each "stigma" is really its own flower with its own ovule
# I'll combine flower and stigma and call that combo "flower"
TD_pre <- TD %>% 
  mutate(Flower=paste(Flower,Stigma, sep="-")) %>% 
  mutate(Flw_pollen = Pollen, Infl_min=NA, Infl_max=Height, Flower_height=NA, Stigmas_per_flw=1) %>% 
  select(all_of(sel_vec))
# rbind(RA_pre, TD_pre) # test 

# combine all species data together
stig_CS2021 <- rbind(AA_pre, AR_pre, CA_pre, DL_pre, DI_pre, PL_pre,PP_pre, SV_pre, RA_pre, TD_pre)

## Check data
stig <- stig_CS2021
# head(stig)
# summary(stig)
stig[which(is.na(stig$Flw_pollen)),] # looks like all the population data from the QMP pop of P. lanceolata
# is there but pollen counts are missing! I'll pull up the original stigma data file for P. lanceolata and add that in now
PL_stig <- read.csv("raw-data/Plantago_lanceolata_stigmas.csv", stringsAsFactors = T)
head(PL_stig)
summary(PL_stig) # QMP is the site I'm looking for
PL_site <- read.csv("raw-data/Plantago_lanceolata_site.csv", stringsAsFactors = T)
head(PL_site)
summary(PL_site) # QMP is called "QU" here - I'll replace it with QMP
PL_site$Population <- as.factor(gsub("QU", "QMP", PL_site$Population))
# I just want the QMP sites...
PL_stig_QMP <- PL_stig %>% filter(Site=="QMP")
PL_site_QMP <- PL_site %>% filter(Population=="QMP")
# the individuals in "stig" have an S in front of them, which will get in the way when I join them together
PL_stig_QMP$Plant <- gsub("S", "", PL_stig_QMP$Plant)
PL_site_QMP$Plant <- as.character(PL_site_QMP$Plant)
# join them together
PL_QMP <- full_join(PL_stig_QMP, PL_site_QMP, by=c("Species", "Date", "Site"="Population", "Plant"))
PL_QMP <- PL_QMP %>% 
  mutate(Flw_pollen = Pollen, Stigmas_per_flw=1, Flower_height=NA, Infl_min=Min, Infl_max=Max) %>% 
  select(all_of(sel_vec))
# now take out old PL QMP data from "stig"
stig_noqmp <- stig %>% filter(Site != "QMP")
stig_newqmp <- rbind(stig_noqmp, PL_QMP)

# Add in Plantago major
PM_stig <- read.csv("raw-data/Plantago_major_stigmas.csv", stringsAsFactors = T)
head(PM_stig)
summary(PM_stig) 
PM_site <- read.csv("raw-data/Plantago_major_site.csv", stringsAsFactors = T)
head(PM_site)
summary(PM_site) # QMP is called "QU" here - I'll replace it with QMP
# join them together
PM_site$Plant <- as.character(PM_site$Plant)
PM_stig$Plant <- as.character(PM_stig$Plant)
PM_stig$Site <- "Biosci"
PM <- full_join(PM_stig, PM_site, by=c("Species", "Site"="Population", "Plant", "Flower"))
PM_pre <- PM %>% 
  mutate(Flw_pollen = Pollen, Stigmas_per_flw=1, Infl_min=Flower_min, Infl_max=Flower_max, Date=Date.x) %>% 
  select(all_of(sel_vec))
# join this onto stig_newqmp 
stig_pre <- rbind(stig_newqmp, PM_pre)

# write to a csv
write.csv(stig_pre, "processed-data/stig-CS2021.csv", row.names = F)

