## Combine CS anther size data into one file
## Claire Smith
## Last updated: 22 June 2023

# Create size-CS2021-all.csv, a file with all the CS2021 pollen size data compiled together

##########################################################################################
# Load libraries
library(tidyverse)

##########################################################################################

# For each species
# *Load in data
# *Separate out individual and rep labels
# *Combine into dataframe with columns Ind, Rep, Species, Site

# Aart - Ambrosia artemisiifolia
# Load in data
aart_fbf_raw <- read.csv("raw-data/size-A.art.area.FBF 2021-08-14.unfiltered.csv")
aart_lp_raw <- read.csv("raw-data/size-A.art.area.LP 2021-08-13.unfiltered.csv")

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
aart_fbf_ind_rep <- gsub("^Aart_Area_(\\d+-\\d+).+", "\\1", aart_fbf_raw$Filename)
aart_fbf_ind_rep_split <- str_split_fixed(aart_fbf_ind_rep, "-", 2)
aart_fbf <- as.data.frame(aart_fbf_ind_rep_split)
names(aart_fbf) <- c("Ind","Rep")
aart_fbf$Area <- aart_fbf_raw$Area
aart_fbf$Species <- "Ambrosia artemisiifolia"
aart_fbf$Site <- "FBF"
aart_fbf$Date <- "2021-08-14"

aart_lp_ind_rep <- gsub("^Aart_Area_(\\d+-\\d+).+", "\\1", aart_lp_raw$Filename)
aart_lp_ind_rep_split <- str_split_fixed(aart_lp_ind_rep, "-", 2)
aart_lp <- as.data.frame(aart_lp_ind_rep_split)
names(aart_lp) <- c("Ind","Rep")
aart_lp$Area <- aart_lp_raw$Area
aart_lp$Species <- "Ambrosia artemisiifolia"
aart_lp$Site <- "LP"
aart_lp$Date <- "2021-08-13"

aart <- rbind(aart_fbf, aart_lp)

# aart %>% group_by(Site, Ind, Rep) %>% summarize(n=n())

# Amsp - Amaranthus retroflexus
# Load in data
amsp_fbf0717_raw <- read.csv("raw-data/size-A.sp_area.FBF20210717.unfiltered.csv")
amsp_fbf0816_raw <- read.csv("raw-data/size-A.sp_area.FBF20210816.unfiltered.csv")

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
amsp_fbf0717_ind_rep <- gsub("^Asp_Area_A?(.+-\\d+).+", "\\1", amsp_fbf0717_raw$Filename) # had to modify this because there are Xs
amsp_fbf0717_ind_rep_split <- str_split_fixed(amsp_fbf0717_ind_rep, "-", 2)
amsp_fbf0717 <- as.data.frame(amsp_fbf0717_ind_rep_split)
names(amsp_fbf0717) <- c("Ind","Rep")
amsp_fbf0717$Area <- amsp_fbf0717_raw$Area
amsp_fbf0717$Species <- "Amaranthus retroflexus"
amsp_fbf0717$Site <- "FBF"
amsp_fbf0717$Date <- "2021-07-17"

# I'll call this separate date fbf2 to distinguish in case there are repeat numbers
amsp_fbf0816_ind_rep <- gsub("^Asp_Area_A?(.+-\\d+).+", "\\1", amsp_fbf0816_raw$Filename) # had to modify this because there are Xs
amsp_fbf0816_ind_rep_split <- str_split_fixed(amsp_fbf0816_ind_rep, "-", 2)
amsp_fbf0816 <- as.data.frame(amsp_fbf0816_ind_rep_split)
names(amsp_fbf0816) <- c("Ind","Rep")
amsp_fbf0816$Area <- amsp_fbf0816_raw$Area
amsp_fbf0816$Species <- "Amaranthus retroflexus"
amsp_fbf0816$Site <- "FBF2"
amsp_fbf0816$Date <- "2021-08-16"

amsp <- rbind(amsp_fbf0717, amsp_fbf0816)

# Chsp - Chenopodium album
# Load in data
chsp_fbf1_raw <- read.csv("raw-data/size-C.sp_area.FBF20210721.unfiltered.csv")
chsp_fbf2_raw <- read.csv("raw-data/size-C.sp_area.FBF20210723.unfiltered.csv")

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
chsp_fbf1_ind_rep <- gsub("^C.sp_Area_M(\\d+-\\d+).+", "\\1", chsp_fbf1_raw$Filename)
chsp_fbf1_ind_rep_split <- str_split_fixed(chsp_fbf1_ind_rep, "-", 2)
chsp_fbf1 <- as.data.frame(chsp_fbf1_ind_rep_split)
names(chsp_fbf1) <- c("Ind","Rep")
chsp_fbf1$Area <- chsp_fbf1_raw$Area
chsp_fbf1$Species <- "Chenopodium album"
chsp_fbf1$Site <- "FBF1"
chsp_fbf1$Date <- "2021-07-21"

chsp_fbf2_ind_rep <- gsub("^Csp_Area_M(\\d+-\\d+).+", "\\1", chsp_fbf2_raw$Filename)
chsp_fbf2_ind_rep_split <- str_split_fixed(chsp_fbf2_ind_rep, "-", 2)
chsp_fbf2 <- as.data.frame(chsp_fbf2_ind_rep_split)
names(chsp_fbf2) <- c("Ind","Rep")
chsp_fbf2$Area <- chsp_fbf2_raw$Area
chsp_fbf2$Species <- "Chenopodium album"
chsp_fbf2$Site <- "FBF2"
chsp_fbf2$Date <- "2021-07-23"

chsp <- rbind(chsp_fbf1, chsp_fbf2)

# chsp %>% group_by(Site, Ind) %>% summarize(n=n())

# Plan - Plantago lanceolata
# Load in data
plan_fbf_raw <- read.csv("raw-data/size-P.lan_area.FBF20210818.unfiltered.csv")
plan_km_raw <- read.csv("raw-data/size-P.lan_areaKM20210731.unfiltered.csv")
plan_qmp_raw <- read.csv("raw-data/size-P.lan_areaQMP20210820.unfiltered.csv")

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
plan_fbf_ind_rep <- gsub("^Plan_size_.(\\d+-\\d+).+", "\\1", plan_fbf_raw$Filename)
plan_fbf_ind_rep_split <- str_split_fixed(plan_fbf_ind_rep, "-", 2)
plan_fbf <- as.data.frame(plan_fbf_ind_rep_split)
names(plan_fbf) <- c("Ind","Rep")
plan_fbf$Area <- plan_fbf_raw$Area
plan_fbf$Species <- "Plantago lanceolata"
plan_fbf$Site <- "FBF"
plan_fbf$Date <- "2021-08-18"

plan_km_ind_rep <- gsub("Plan_Area_M?(\\d+\\D)-(\\d+)_.+", "\\1-\\2", plan_km_raw$Filename)
plan_km_ind_rep_split <- str_split_fixed(plan_km_ind_rep, "-", 2)
plan_km <- as.data.frame(plan_km_ind_rep_split)
names(plan_km) <- c("Ind","Rep")
plan_km$Area <- plan_km_raw$Area
plan_km$Species <- "Plantago lanceolata"
plan_km$Site <- "KM"
plan_km$Date <- "2021-07-31"

plan_qmp_ind_rep <- gsub("Plan_area_.(\\d+-\\d+)_.+", "\\1", plan_qmp_raw$Filename)
plan_qmp_ind_rep_split <- str_split_fixed(plan_qmp_ind_rep, "-", 2)
plan_qmp <- as.data.frame(plan_qmp_ind_rep_split)
names(plan_qmp) <- c("Ind","Rep")
plan_qmp$Area <- plan_qmp_raw$Area
plan_qmp$Species <- "Plantago lanceolata"
plan_qmp$Site <- "QMP"
plan_qmp$Date <- "2021-08-20"

plan <- rbind(plan_fbf, plan_km)
plan <- rbind(plan, plan_qmp)

# plan %>% group_by(Site, Ind) %>% summarize(n=n())

# PSG - Setaria viridis
# PSG (manually measured)
PSG_raw <- read.csv("raw-data/size-PSG_area20210718.manual.csv")

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site

PSG <- data.frame(Ind=PSG_raw$X)
PSG$Area <- PSG_raw$Area
PSG$Species <- "Setaria viridis"
PSG$Site <- "FBF"
PSG$Date <- "2021-07-18"

PSG$Ind <- as.character(PSG$Ind)
PSG$Rep <- NA

# R acetosella
race_raw <- read.csv("raw-data/size-R.ace.area.unfiltered.csv")
# head(race_raw)
# View(race_raw)

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
race_LASvec <- gsub("^Race_size_Race_2uL_(L\\d+A\\D+S\\d+).+", "\\1", race_raw$Filename)
race <- data.frame(Area=race_raw$Area)
lvec <- gsub("^L(\\d+)A.*", "\\1", race_LASvec) # plant number
avec <- gsub("^L\\d+A(\\D+)S\\d+", "\\1", race_LASvec) # anther (within flower) number
svec <- gsub("^L\\d+A\\D+S(\\d+)", "\\1", race_LASvec) # subsample number

race$Ind <- lvec
race$Anther <- avec
race$Anther_rep <- svec
race$Species <- "Rumex acetosella"
race$Site <- "KM"
race$Date <- "2021-07-08"

# Tdio
tdio_raw <- read.csv("raw-data/size-TDio_filtered.csv")
# head(tdio_raw)
# View(tdio_raw)

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
tdio_LASvec <- gsub("^Tdio_size_Tdio_2-5uL_(L\\d+A\\d+S\\d+).+", "\\1", tdio_raw$Filename)
tdio <- data.frame(Area=tdio_raw$Area)
lvec <- gsub("^L(\\d+)A.*", "\\1", tdio_LASvec)
avec <- gsub("^L\\d+A(\\d+)S\\d+", "\\1", tdio_LASvec)
svec <- gsub("^L\\d+A\\d+S(\\d+)", "\\1", tdio_LASvec)

tdio$Ind <- lvec
tdio$Anther <- avec
tdio$Anther_rep <- svec
tdio$Species <- "Thalictrum dioicum"
tdio$Site <- "BI"
tdio$Date <- "2021-05-07"

# Combining species into 1 file
# names(aart)
# names(amsp)
# names(chsp)
# names(plan)
# names(PSG)

# Only take "Ind", "Area", "Species", "Site"
keep_vec <- c("Species","Ind","Site","Area")
rb1 <- rbind(aart[keep_vec],amsp[keep_vec])
rb2 <- rbind(rb1, chsp[keep_vec])
rb3 <- rbind(rb2, plan[keep_vec])
rb4 <- rbind(rb3, PSG[keep_vec])
rb5 <- rbind(rb4, tdio[keep_vec])
rb6 <- rbind(rb5, race[keep_vec])

# summary(rb5)
# subset(rb5, Area>600)

size_dat <- rb6

# Make pollen diam column
size_dat$Diam <- 2*sqrt(size_dat$Area/pi)

# Add source column
size_dat <- size_dat %>% mutate(source="CS2021")

# Add sex system info
size_dat_ss <- size_dat %>% 
  mutate(Sex_sys = case_when(source == "JF2001" ~ "hermaphroditic",
                             Species %in% c("Rumex acetosella","Thalictrum dioicum") ~ "dioecious",
                             Species %in% c("Ambrosia artemisiifolia", "Amaranthus retroflexus", "Carex communis",
                                            "Carex hirtifolia","Carex plantaginea",
                                            "Carex pedunculata",
                                            "Carex stipata",
                                            "Rumex crispus",
                                            "Scirpus microcarpus") ~ "monoecious",
                             Species %in% c("Chenopodium album", "Dichanthelium linearifolium", 
                                            "Dichanthelium implicatum", "Phleum pratense", 
                                            "Setaria viridis", "Schizacne purpurascens",
                                            "Elymus repens", "Agropyron trachycaulum",
                                            "Festuca rubra", "Festuca pratensis", 
                                            "Festuca campestris",
                                            "Avenula hookeri", "Hierochloe odorata",
                                            "Koeleria cristata","Phalaris arundinacea",
                                            "Plantago lanceolata", 
                                            "Poa juncifolia", "Poa secunda subsp. secunda",
                                            "Stipa columbiana") ~ "hermaphroditic"))
# Make sure no species were missed -- looks ok 
# summary(as.factor(size_dat_ss$Sex_sys))

# # Check data
size <- size_dat_ss
# size %>% group_by(Species, Site) %>% 
#   summarize(n=n())

# Ssome repeat anther measurements for an individual are listed as "1A" and "1B" etc -- remove the A's and B's
# to get overall individual average pollen size
# unique(size$Ind)
size <- size %>% mutate(Ind = gsub("[ABCD]", "", Ind))
# unique(size$Ind) # looks good!


# Write file that has within-individual variation in pollen size
# Each row is a pollen size measurement within an individual
write.csv(size, "processed-data/size-CS2021-within-inds.csv", row.names = F)

# Take the mean pollen size per individual, write a file where each individual is a row
avg_size_dat <- size %>% 
  filter(Species != "Setaria viridis") %>% # S viridis doesn't have individual-level info - I'll take it out
  group_by(source, Sex_sys, Species, Site, Ind) %>% 
  summarize(Avg_area = mean(Area, na.rm=T),
            Avg_diam = mean(Diam, na.rm=T), 
            Sd_diam = sd(Diam, na.rm=T),
            N_diam = n())
write.csv(avg_size_dat, "processed-data/size-CS2021.csv", row.names=F)

# Save S. viridis data in its own file to attach to later species-level tables and analyses
SV_sp_dat <- size %>% filter(Species == "Setaria viridis") %>% 
  group_by(source, Sex_sys, Species) %>% 
  summarize(Avg_area = mean(Area, na.rm=T),
            Avg_diam = mean(Diam, na.rm=T), 
            Sd_diam = sd(Diam, na.rm=T),
            N_diam = n())
write.csv(SV_sp_dat, "processed-data/SV-sp-dat.csv", row.names=F)
