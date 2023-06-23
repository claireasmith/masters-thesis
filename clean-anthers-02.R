## Combine CS pollen number data into one file
## Claire Smith
## Last updated: 20 June 2023

# Create prod-CS2021-all.csv, a file with all the CS2021 pollen production data compiled together

############################################################################################
# Load libraries
library(tidyverse) # need dplyr and tidyr specifically

############################################################################################
# Load in data
# First take the "manual" count data - I'll add in the "automated" count data later in the script
aAA <- read.csv("raw-data/prod-Ambrosia_artemisiifolia_anthers - Sheet1.csv")
aAR <- read.csv("raw-data/prod-Amaranthus_retroflexus_anthers - Sheet1.csv")
aCA <- read.csv("raw-data/prod-Chenopodium_album_anthers - Sheet1.csv")
aPL <- read.csv("raw-data/prod-Plantago_lanceolata_anthers - Sheet1.csv")
aSV <- read.csv("raw-data/prod-Purple stigma grass - Sheet1 (1).csv")

# Date was missing for these sites in Ambrosia artemisiifolia data - add in
aAA <- aAA %>% mutate(Date = case_when(Site == "Fruition Berry Farm" ~ "14 Aug 2021",
                                       Site =="Lemoine Point" ~ "20 Sept 2021"))

# Take relevant rows, make them all matching
aAA <- dplyr::select(aAA, Species, Date, Site, Label, Anther, Subsample, Tot_tube_vol_uL, Subsample_vol_uL, Pollen_sub, Observer, N_anth_tube, Anth_per_flw)
aAR <- dplyr::select(aAR, Species, Date, Site, Label, Anther, Subsample, Tot_tube_vol_uL, Subsample_vol_uL, Pollen_sub, Observer, N_anth_tube, Anth_per_flw)
aCA <- dplyr::select(aCA, Species, Date, Site, Label, Anther, Subsample, Tot_tube_vol_uL, Subsample_vol_uL, Pollen_sub, Observer, N_anth_tube, Anth_per_flw)
aPL <- dplyr::select(aPL, Species, Date, Site, Label, Anther, Subsample, Tot_tube_vol_uL, Subsample_vol_uL, Pollen_sub, Observer, N_anth_tube, Anth_per_flw)
aSV <- dplyr::select(aSV, Species, Date, Site, Label, Anther, Subsample, Tot_tube_vol_uL, Subsample_vol_uL, Pollen_sub, Observer, N_anth_tube, Anth_per_flw)
# There are a bunch of empty rows in the AR dataframe for some reason - remove these
aAR <- filter(aAR, !is.na(Pollen_sub))

# Bind together all species anther counts into one dataframe
aDat_pre <- rbind(aAA, aAR, aCA, aPL, aSV)
# View(aDat_pre)

# Calculate average per-anther pollen production per individual 
aDat <- aDat_pre %>% 
  dplyr::mutate(polanth=Pollen_sub/Subsample_vol_uL*Tot_tube_vol_uL/N_anth_tube) 

# Add in pollen counting method
aDat <- aDat %>% mutate(count_method = "manual")

#####################################################################################################################
## Automated counts

# Rumex acetosella
autoRA_raw <- read.csv("raw-data/prod-2022-09-27_R.ace.count.filtered_2SD.csv")

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
ra_count_LASvec <- gsub("^Race_count_Race_2uL_(L\\d+A\\D+S\\d+).+", "\\1", autoRA_raw$X)
aRA_auto <- data.frame(Pollen_sub=autoRA_raw$x)
lvec <- gsub("^L(\\d+)A.*", "\\1", ra_count_LASvec) # plant number
avec <- gsub("^L\\d+A(\\D+)S\\d+", "\\1", ra_count_LASvec) # anther (within flower) number
svec <- gsub("^L\\d+A\\D+S(\\d+)", "\\1", ra_count_LASvec) # subsample number

aRA_auto$Label <- lvec
aRA_auto$Anther <- avec
aRA_auto$Subsample <- svec
aRA_auto$Species <- "Rumex acetosella"
aRA_auto$Site <- "KM" 
aRA_auto$Tot_tube_vol_uL <- 250
aRA_auto$Subsample_vol_uL <- 2
aRA_auto$N_anth_tube <- 6
aRA_auto$Anth_per_flw <- 6
aRA_auto$count_method <- "auto"
aRA_auto$Date <- "8 July 2022"
aRA_auto$Observer <- "Jeremiah"
# head(aRA_auto)
# head(aDat)

aRA_auto <- aRA_auto %>% 
  dplyr::mutate(polanth=Pollen_sub/Subsample_vol_uL*Tot_tube_vol_uL/N_anth_tube)
# names(aRA_auto)
# names(aDat)

aRA_auto_inorder <- aRA_auto[names(aDat)]

aDat_plusRA <- rbind(aDat, aRA_auto_inorder)

# Thalictrum dioicum 
autoTD_raw <- read.csv("raw-data/prod-2022-09-27_T.dio.count.filtered_2SD.csv")

# Separate out individual and rep labels
# Combine into dataframe with Ind, Rep, Species, Site
TD_count_LASvec <- gsub("^Tdio_count_Tdio_2-5uL_(L\\d+A\\d+S\\d+).+", "\\1", autoTD_raw$X)
aTD_auto <- data.frame(Pollen_sub=autoTD_raw$x)
lvec <- gsub("^L(\\d+)A.*", "\\1", TD_count_LASvec) # plant number
avec <- gsub("^L\\d+A(\\d+)S\\d+", "\\1", TD_count_LASvec) # anther (within flower) number
svec <- gsub("^L\\d+A\\d+S(\\d+)", "\\1", TD_count_LASvec) # subsample number

aTD_auto$Label <- lvec
aTD_auto$Anther <- avec
aTD_auto$Subsample <- svec
aTD_auto$Species <- "Thalictrum dioicum"
aTD_auto$Site <- "BI"
aTD_auto$Tot_tube_vol_uL <- 450
aTD_auto$Subsample_vol_uL <- 2.5
aTD_auto$N_anth_tube <- 1

aTD_auto$count_method <- "auto"
aTD_auto$Date <- "7 May 2022"
aTD_auto$Observer <- "Jeremiah"

aTD_auto <- aTD_auto %>% mutate(Anth_per_flw = case_when(Label == 1 ~ 28,
                                                         Label == 2 ~ 32,
                                                         Label == 3 ~ 26,
                                                         Label == 4 ~ 19,
                                                         Label == 5 ~ 22,
                                                         Label == 6 ~ 27,
                                                         Label == 8 ~ 16,
                                                         Label == 9 ~ 27, 
                                                         Label == 10 ~ 26))

aTD_auto <- aTD_auto %>% 
  dplyr::mutate(polanth=Pollen_sub/Subsample_vol_uL*Tot_tube_vol_uL)

aTD_auto_inorder <- aTD_auto[names(aDat)]

aDat_full <- rbind(aDat_plusRA, aTD_auto_inorder)
# View(aDat_full)

# Add source column
aDat_full <- aDat_full %>% mutate(source="CS2021")

# Add sex system info
aDat_full_ss <- aDat_full %>% 
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
# summary(as.factor(aDat_full_ss$Sex_sys))
# aDat_full_ss[which(is.na(aDat_full_ss$Sex_sys)),]

# Check data
prod <- aDat_full_ss
# prod %>% group_by(Species, Site, Date, Label) %>% 
#   summarize(n=n()) %>% 
#   group_by(Species, Site, Date) %>% 
#   summarize(n=n())

# Some species have multiple sampling dates from the same sites. I know that some of the numbers repeat - I want
# to make sure individuals from the same site but different dates don't get combined. I also know in the
# size data some of these sites have been called FBF2 e.g. for the second date to distinguish them. I'll do that too. 
# later A. retro entry (16 Aug 2021) needs to become FBF2, later C. album entry (23 Jul 2021) needs to become FBF2

prod$Site[which(prod$Species=="Amaranthus retroflexus"&prod$Date=="16 Aug 2021")] <- "FBF2"
prod$Site[which(prod$Species=="Chenopodium album"&prod$Date=="22 Jul 2021")] <- "FBF1"
prod$Site[which(prod$Species=="Chenopodium album"&prod$Date=="23 Jul 2021")] <- "FBF2"
# Replace full site names (Lemoine Point, Fruition Berry Farm) with their codes for consistency
# unique(aDat$Site)
sites1 <- gsub("Lemoine Point", "LP", prod$Site)
# unique(sites1)
prod$Site <- sites1
sites2 <- gsub("Fruition Berry Farm", "FBF", prod$Site)
# unique(sites2)
prod$Site <- sites2

# Species names:
# unique(aDat$Species)
# Replace code "Purple stigma grass" with actual species name, "Setaria viridis"
spec <- gsub("Purple stigma grass", "Setaria viridis", prod$Species)
# unique(test)
prod$Species <- spec

# Check it all worked: 
# prod %>% group_by(Species, Site, Date, Label) %>% 
#   summarize(n=n()) %>% 
#   group_by(Species, Site, Date) %>% 
#   summarize(n=n())

# Add source = "CS2021" to data
prod <- prod %>% mutate(source="CS2021")

# Write file that includes within-individual variation in pollen prod
write.csv(prod, "processed-data/prod-CS2021-within-inds.csv", row.names=F)

# Calculate per-individual pollen production per anther
prod_avg <- prod %>% 
  group_by(source, Sex_sys, Species, Date, Site, Label, Observer, Anth_per_flw, count_method) %>% 
  summarize(Avg_pol_anth=mean(polanth, na.rm=T))

write.csv(prod_avg, "processed-data/prod-CS2021.csv", row.names=F)

