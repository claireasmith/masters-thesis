## Combine CS pollen size and number data with JF pollen size and number data
## Claire Smith
## Last updated: 26 Jul 2023

# Create size-prod-all.csv, which combines size-prod-CS2021.csv with pollen size and number data
# from JF2001_pollen and JF_2004_anthers_all_spp

############################################################################################
# Load libraries
library(tidyverse) # need dplyr and tidyr specifically
library(readxl)

############################################################################################
# Read in data
spCS2021 <- read.csv("processed-data/size-prod-CS2021.csv", stringsAsFactors = T)
# head(spCS2021)
spJF2001raw <- read_xls("raw-data/JF2001_pollen.xls")
# head(spJF2001raw)
spJF2004raw <- read_xls("raw-data/JF_2004_anthers_all_spp.xls")
# head(spJF2004raw)

sel_vec1 <- c("spp", "ind", "nanthers","polcnt","pollmsize","sdsize", "source", "Date")

## 2001 data prep
spJF2001 <- spJF2001raw %>% 
  # keep only relevant columns
  dplyr::select(spp, ind, nanthers, polcnt, pollmsize, sdsize, nanthers, date) %>% 
  # keep track of which dataset they came from 
  mutate(source = "JF2001",
         Date = as.Date(date)) %>% 
  filter(!is.na(polcnt)) %>% 
  select(all_of(sel_vec1))
# head(spJF2001)
# names(spJF2001)


## 2004 data prep
# Convert dates from jday to date
spJF2004pre <- spJF2004raw
d2004 <- as.character(spJF2004raw$jday)
d2004_jday <- as.numeric(gsub("^4", "", d2004)) # they are julian dates with a prefix "4"
d2004_date <- as.Date(d2004_jday,    # Convert Julian day to date
                      origin = as.Date("2004-01-01"))
spJF2004pre$Date <- as.Date(d2004_date) # add back to dataframe

spJF2004 <- spJF2004pre %>% 
  mutate(source = "JF2004",
         # Keep col names consistent with JF2001
         nanthers = nanth,
         pollmsize=meansize) %>% 
  filter(!is.na(polcnt)) %>% 
  filter(nsample>1) %>% 
  select(all_of(sel_vec1))

# Join JF datasets together
JF2001_2004 <- rbind(spJF2001, spJF2004)
# View(JF2001_2004)
JF2001_2004$count_method = "particle counter" # keep track of how pollen was counted

# names(spCS2021)
# [1] "source"       "Sex_sys"      "Species"      "Site"         "Ind"          "Avg_area"     "Avg_diam"   "Sd_diam"   "Date"         
# "Label"        "Observer"     "Anth_per_flw" "count_method" "Avg_pol_anth" "Sd_pol_anth" 

# names(JF2001_2004)
# [1] "spp"          "ind"          "nanthers"     "polcnt"       "pollmsize"    "sdsize"  "source"       "count_method"

sel_vec <- c("source", "Sex_sys", "Species", "Site", "Ind", "Label", "Avg_diam", "Sd_diam", "N_diam",
             "Avg_area", "Avg_pol_anth", "Sd_pol_anth", "Anth_per_flw", "count_method", "Anth_per_sample",
             "Date")

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
         N_diam = polcnt, # since all grains counted and measured, the # of grains diam is measured from 
         # is just the number of grains in the sample
         Avg_area = NA,
         Anth_per_sample = nanthers) %>% # could add this in per species later 
  select(all_of(sel_vec))

CS2021_pre <- spCS2021 %>% 
  mutate(Anth_per_sample = Anth_per_flw,
         Date = dmy(spCS2021$Date)) %>% 
  select(all_of(sel_vec))


# Bind them together
sp_all <- rbind(JF2001_2004_pre, CS2021_pre)

### Cleaning up joined pollen number data

## Keep only species with at least 3 individuals
sp_sum <- sp_all %>% group_by(Species) %>% dplyr::summarize(n=n())
keep_vec <- sp_sum$Species[which(sp_sum$n>=3)]
sp_3ormore <- sp_all[sp_all$Species %in% keep_vec,]

# Update names to full species names to full versions
sp_dat <- sp_3ormore %>% 
  mutate(Species=stringr::str_replace_all(Species, c("A.artemisi"="Ambrosia artemisiifolia",
                                                     "Astolonifera"="Agrostis stolonifera",
                                                     "Atrachg"="Elymus trachycaulus",
                                             "Binermis"="Bromus inermis",
                                             "Bcarianatus"="Bromus carianatus",
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
# Add anthers per flower if not already there
# Add sex system data
sp_dat_anth <- sp_dat_ss %>% 
  mutate(Anth_per_flw = case_when(source == "JF2001" ~ 3, # all JF2001 species are grasses with 3 anthers
                                  Species %in% c("Carex communis",
                                                 "Carex hirtifolia","Carex plantaginea",
                                                 "Carex pedunculata",
                                                 "Carex stipata",
                                                 "Scirpus microcarpus",
                                                 "Dicanthelium sp 1", 
                                                 "Dicanthelium sp 2", "Phleum pratense", 
                                                 "Setaria viridis", "Schizachne purpurascens",
                                                 "Elymus repens", "Agropyron trachycaulum",
                                                 "Festuca rubra", "Festuca pratensis", 
                                                 "Festuca campestris",
                                                 "Avenula hookeri", "Hierochloe odorata",
                                                 "Koeleria cristata","Phalaris arundinacea",
                                                 "Poa juncifolia", "Poa secunda subsp. secunda",
                                                 "Stipa columbiana") ~ 3,
                                  Species %in% c("Plantago lanceolata") ~ 4,
                                  Species %in% c("Amaranthus retroflexus",
                                                 "Ambrosia artemisiifolia",
                                                 "Chenopodium album") ~ 5,
                                  Species %in% c("Rumex acetosella", "Rumex crispus") ~ 6,
                                  Species %in% c("Thalictrum dioicum") ~ Anth_per_sample))

# View(sp_dat_ss)
# View(sp_dat_anth[which(is.na(sp_dat_anth$Anth_per_flw)),]) # none NA

## Add families
# read in file with species names and families
famdat <- read.csv("raw-data/species_list_sizenum.csv", stringsAsFactors = T)
sp_dat_fam <- left_join(sp_dat_anth, famdat, by=c("Species"="species"))
# View(sp_dat_fam[which(is.na(sp_dat_fam$family)),]) # none NA


# Remove any repeat lines if they exist
sp_dat2 <- sp_dat_fam %>% distinct() 

sp_dat3 <- sp_dat2 %>% mutate(Pol_flw = Avg_pol_anth*Anth_per_flw,
                            Sd_pol_flw = Sd_pol_anth*Anth_per_flw)

# After looking at ridgeplot distributions of pollen production it looks like there are some clear outliers in some species - like Phleum pratense
# prodfull %>% filter(Species == "Phleum pratense")
# looks like Ind C11 is one, with 16161 grains/anther! The rest are all <10,000. I'll remove it
sp_dat3 <- sp_dat3 %>% 
  # remove outliers 
  filter(!(Species=="Phleum pratense"&Ind=="C11")) %>% 
  filter(!(Species=="Elymus innovatus"&Ind=="EI01"))
# prodfull_filt %>% filter(Species == "Phleum pratense") # make sure it's gone

write.csv(sp_dat3, "processed-data/size-prod-all.csv", row.names = F)

# Some species were collected by both me and Jannice. She used automated pollen counting 
# with a particle counter. It counts all particles that go through it vs small subsample I 
# counted in mine, makes less assumptions about how uniformly pollen is suspended in the sample -- I'll keep her data where there's overlap. 

# Are any species in my production data NOT also included in the JF data?
(unique(sp_dat3$Species[which(sp_dat3$source == "CS2021")]))[!(unique(sp_dat3$Species[which(sp_dat3$source == "CS2021")])) %in% unique(sp_dat3$Species[which(sp_dat3$source == "JF2001" | sp_dat3$source == "JF2004" )])]
# Amaranthus retroflexus is the only one
# Exclude all pollen production data from CS2021 except Amaranthus retroflexus
sizeprod_norep <- sp_dat3 %>% filter(source != "CS2021" |  Species == "Amaranthus retroflexus")

write.csv(sizeprod_norep, "processed-data/size-prod-norep.csv", row.names = F)


## Do the Cyperaceae produce more pollen or smaller pollen? 
# cyper <- sizeprod_norep %>% filter(family=="Cyperaceae") %>% 
#   group_by(Species) %>% 
#   summarize(m_avg_diam = mean(Avg_diam),
#             sd_avg_diam = sd(Avg_diam),
#             n = n()) %>% 
#   summarize(m_diam = mean(m_avg_diam))
# summary(cyper)
# noncyper <- sizeprod_norep %>% filter(family!="Cyperaceae") %>% 
#   group_by(Species) %>% 
#   summarize(m_avg_diam = mean(Avg_diam, na.rm=T),
#             sd_avg_diam = sd(Avg_diam),
#             n = n()) %>% 
#   summarize(m_diam = mean(m_avg_diam, na.rm=T))
# 
# cyper2 <- sizeprod_norep %>% filter(family=="Cyperaceae") %>% 
#   group_by(Species) %>% 
#   summarize(m_pol_anth = mean(Avg_pol_anth),
#             sd_pol_anth = sd(Avg_pol_anth),
#             n = n()) %>% 
#   summarize(m_pol = mean(m_pol_anth))
# noncyper2 <- sizeprod_norep %>% filter(family!="Cyperaceae") %>% 
#   group_by(Species) %>% 
#   summarize(m_pol_anth = mean(Avg_pol_anth),
#             sd_pol_anth = sd(Avg_pol_anth),
#             n = n()) %>% 
#   summarize(m_pol = mean(m_pol_anth))
