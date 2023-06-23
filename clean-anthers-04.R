## Combine CS pollen size and number data with JF pollen size and number data
## Claire Smith
## Last updated: 22 June 2023

# Create size-prod-all.csv, which combines size-prod-CS2021.csv with pollen size and number data
# from JF2001_pollen and JF_2004_anthers_all_spp

############################################################################################
# Load libraries
library(tidyverse) # need dplyr and tidyr specifically
library(readxl)

# Source files 
source("clean-anthers-01.R")
source("clean-anthers-02.R")
source("clean-anthers-03.R")

############################################################################################
# Read in data
spCS2021 <- read.csv("processed-data/size-prod-CS2021.csv", stringsAsFactors = T)
# head(spCS2021)
spJF2001raw <- read_xls("raw-data/JF2001_pollen.xls")
# head(spJF2001raw)
spJF2004raw <- read_xls("raw-data/JF_2004_anthers_all_spp.xls")
# head(spJF2004raw)

# 2001 data
spJF2001 <- spJF2001raw %>% 
  # keep only relevant columns
  dplyr::select(spp, ind, nanthers, polcnt, pollmsize, sdsize, nanthers) %>% 
  # keep track of which dataset they came from 
  mutate(source = "JF2001") %>% 
  filter(!is.na(polcnt))
# head(spJF2001)
# names(spJF2001)

# 2004 data
spJF2004 <- spJF2004raw %>% 
  dplyr::select(spp, ind, nanth, polcnt, meansize, sdsize, nanth) %>% 
  mutate(source = "JF2004") %>% 
  filter(!is.na(polcnt))
# head(spJF2004)
# Keep col names consistent with JF2001
names(spJF2004) <- c("spp", "ind", "nanthers", "polcnt", "pollmsize", "sdsize", "source")

# Join JF datasets together
JF2001_2004 <- rbind(spJF2001, spJF2004)
# View(JF2001_2004)
JF2001_2004$count_method = "particle counter" # keep track of how pollen was counted

# names(spCS2021)
# [1] "source"       "Sex_sys"      "Species"      "Site"         "Ind"          "Avg_area"     "Avg_diam"   "Sd_diam"   "Date"         
# "Label"        "Observer"     "Anth_per_flw" "count_method" "Avg_pol_anth" "Sd_pol_anth" 

# names(JF2001_2004)
# [1] "spp"          "ind"          "nanthers"     "polcnt"       "pollmsize"    "sdsize"  "source"       "count_method"

sel_vec <- c("source", "Sex_sys", "Species", "Site", "Ind", "Label", "Avg_diam", "Sd_diam", "Avg_area", 
             "Avg_pol_anth", "Sd_pol_anth", "Anth_per_flw", "count_method", "Anth_per_sample")

JF2001_2004_pre <- JF2001_2004 %>% 
  mutate(Sex_sys=NA, #fix this later
         Species=spp,
         Ind=ind,
         Avg_pol_anth=polcnt/nanthers,
         Sd_pol_anth=NA,
         Anth_per_flw=NA,
         Site = NA, 
         Label = NA, 
         Avg_diam = pollmsize, 
         Sd_diam = sdsize, 
         Avg_area = NA,
         Anth_per_sample = nanthers) %>% # could add this in per species later 
  select(all_of(sel_vec))

CS2021_pre <- spCS2021 %>% 
  mutate(Anth_per_sample = NA) %>% 
  select(all_of(sel_vec))

# Bind them together
sp_all <- rbind(JF2001_2004_pre, CS2021_pre)

### Cleaning up joined pollen number data

# Keep only species with at least 5 individuals
sp_sum <- sp_all %>% group_by(Species) %>% dplyr::summarize(n=n())
keep_vec <- sp_sum$Species[which(sp_sum$n>=5)]

sp_all_5ormore <- sp_all[sp_all$Species %in% keep_vec,]
sp_dat <- sp_all_5ormore

# Update names to full species names to full versions
sp_dat <- sp_dat %>% 
  mutate(Species=stringr::str_replace_all(Species, c("A.artemisi"="Ambrosia artemisiifolia",
                                             "Binermis"="Bromus inermis",
                                             "C.album"="Chenopodium album",
                                             "C.communis"="Carex communis",
                                             "C.hirtifol"="Carex hirtifolia",
                                             "C.peduncul"="Carex pedunculata",
                                             "C.stipata"="Carex stipata",
                                             "C.plantagi"="Carex plantaginea",
                                             "Einnovatus"="Leymus innovatus",
                                             "Fcampestris"="Festuca campestris",
                                             "P.lanceola"="Plantago lanceolata",
                                             "Ppratense"="Phleum pratense",
                                             "Ppratensis"="Poa pratensis",
                                             "Purple stigma grass"="Setaria viridis",
                                             "R.acetosel"="Rumex acetosella",
                                             "R.crispus"="Rumex crispus",
                                             "S.microcar"="Scirpus microcarpus",
                                             "S.purp"="Schizachne purpurascens",
                                             "T.dioicum"="Thalictrum dioicum",
                                             "Arepens" = "Elymus repens",
                                             "Atrachg" = "Agropyron trachycaulum",
                                             "Frubra" = "Festuca rubra",
                                             "Fpratensis2" = "Festuca pratensis",
                                             "Hhookeri" = "Avenula hookeri",
                                             "Hodorata" = "Hierochloe odorata",
                                             "Kcristata" = "Koeleria cristata",
                                             "Parudinaceae" = "Phalaris arundinacea",
                                             "Pjuncifolia" = "Poa juncifolia",
                                             "purp$" = "Poa secunda subsp. secunda",
                                             "Scolumbiana" = "Stipa columbiana")))

# Add sex system data
sp_dat_ss <- sp_dat %>% 
  mutate(Sex_sys = case_when(source == "JF2001" ~ "hermaphroditic",
                             Species %in% c("Rumex acetosella","Thalictrum dioicum") ~ "dioecious",
                             Species %in% c("Ambrosia artemisiifolia", "Amaranthus retroflexus", "Carex communis",
                                        "Carex hirtifolia","Carex plantaginea",
                                        "Carex pedunculata",
                                        "Carex stipata",
                                        "Rumex crispus",
                                        "Scirpus microcarpus") ~ "monoecious",
                             Species %in% c("Chenopodium album", "Dicanthelium sp 1", 
                                        "Dicanthelium sp 2", "Phleum pratense", 
                                        "Setaria viridis", "Schizachne purpurascens",
                                        "Elymus repens", "Agropyron trachycaulum",
                                        "Festuca rubra", "Festuca pratensis", 
                                        "Festuca campestris",
                                        "Avenula hookeri", "Hierochloe odorata",
                                        "Koeleria cristata","Phalaris arundinacea",
                                        "Plantago lanceolata", 
                                        "Poa juncifolia", "Poa secunda subsp. secunda",
                                        "Stipa columbiana") ~ "hermaphroditic"))

# View(sp_dat_ss)
# View(sp_dat_ss[which(is.na(sp_dat_ss$Sex_sys),)]) # none NA

# Remove any repeat lines if they exist
sp_dat_distinct <- sp_dat_ss %>% distinct() 

## Add in anthers per flower...
# If it comes from JF 2001 it will be a grass, so Anth_per_flw = 3. Sedges also have 3 stamens. 
sp_dat2 <- sp_dat_distinct %>% 
  mutate(Anth_per_flw = case_when(source == "JF2001" ~ 3,
                                  Species == "Ambrosia artemisiifolia" ~ Anth_per_sample,
                                  Species == "Chenopodium album" ~ 5,
                                  source == "JF2004" & Species == "Carex communis" ~ 3,
                                  source == "JF2004" & Species == "Carex hirtifolia" ~ 3,
                                  source == "JF2004" & Species == "Carex pedunculata" ~ 3,
                                  source == "JF2004" & Species == "Carex stipata" ~ 3,
                                  Species == "Plantago lanceolata" ~ 4,
                                  Species == "Rumex acetosella" ~ Anth_per_sample,
                                  source == "JF2004" & Species == "Rumex crispus" ~ Anth_per_sample,
                                  source == "JF2004" & Species == "Scirpus microcarpus" ~ 3,
                                  source == "JF2004" & Species == "Schizachne purpurascens" ~ 3,
                                  source == "JF2004" & Species == "Thalictrum dioicum" ~ Anth_per_sample,
                                  Species == "Amaranthus retroflexus" ~ 5,
                                  TRUE ~ as.numeric(as.character(Anth_per_flw))))


sp_dat3 <- sp_dat2 %>% mutate(Pol_flw = Avg_pol_anth*Anth_per_flw,
                            Sd_pol_flw = Sd_pol_anth*Anth_per_flw)

write.csv(sp_dat3, "processed-data/size-prod-all.csv", row.names = F)
