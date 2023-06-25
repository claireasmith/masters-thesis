## Combine JF and CS pollen capture data
## Claire Smith
## 18 June 2023

# Create stig-all.csv, a file with all stigma data compiled together

##########################################################################################
# Load libraries
library(tidyverse)
library(readxl)
##########################################################################################
# Load files
CS2021 <- read.csv("processed-data/stig-CS2021.csv")
JF2001 <- read.csv("raw-data/JF2001-all.csv")
JF2004 <- read_xls("raw-data/JF2004-all.xls")

CS2021_pre <- CS2021 %>% mutate(Stigma_length=NA, source="CS2021")
sel_vec <- names(CS2021_pre) # I want the columns in the JF datasets to eventually match the ones in the CS one - plus stigma length
# which is in the JF2004 dataset

# I'll go through both JF2001 and JF2004, rename columns that match those in sel_vec
# sel_vec
# [1] "Species"         "Date"            "Site"            "Plant"           "Flower"          "Flw_pollen"      "Stigmas_per_flw"
# [8] "Flower_height"   "Infl_max"        "Infl_min"        "D1"              "D2"              "D3"              "D4"             
# [15] "D5"              "Stigma_length"  "source"  

# JF2001: 
# head(JF2001)
# names(JF2001)
# View(JF2001)
# Species is a code - will need to replace with actual species names. Date would be jday (tho not sure why it's in the 1000s!)
# no site, Plant is ind. pos is position of flower on the plant - would correspond to flower. 
# H1 is max infl height, H2 is min
JF2001_pre <- JF2001 %>% 
  mutate(Species=spp, Date=jday, Site=NA, Plant=ind, Flower=pos, Flw_pollen=npollen, Stigmas_per_flw=nstigma, Flower_height=NA,
         Infl_max=H1, Infl_min=H2, Stigma_length=NA, source="JF2001") %>% 
  select(all_of(sel_vec))

# JF2004: 
# head(JF2004)
# names(JF2004)
# View(JF2004)
JF2004_pre <- JF2004 %>% 
  mutate(Species=spp, Date=jday, Site=NA, Plant=plant, Flower=position, Flw_pollen=npollen, Stigmas_per_flw=NA, Flower_height=NA, 
         Infl_max=height, Infl_min=NA, D1=d1, D2=d2, D3=d3, D4=d4, D5=d5, Stigma_length=stigl, source="JF2004") %>% 
  select(all_of(sel_vec)) %>% 
  distinct() #no duplicate rows
# nrow(JF2001_pre)


# Add in full names for species with shortened names
# unique(CS2021_pre$Species)
CS2021_pre$Species <- gsub("Purple stigma grass", "Setaria viridis", CS2021_pre$Species)
CS2021_pre$Species <- gsub("Amaranthus sp", "Amaranthus retroflexus", CS2021_pre$Species)
CS2021_pre$Species <- gsub("Dicanthelium sp 1", "Dichanthelium linearifolium", CS2021_pre$Species)
CS2021_pre$Species <- gsub("Dicanthelium sp 2", "Dichanthelium implicatum", CS2021_pre$Species)
# unique(JF2001_pre$Species)
JF2001_pre$Species <- gsub("arepens", "Elymus repens", JF2001_pre$Species)
JF2001_pre$Species <- gsub("binermis", "Bromus inermis", JF2001_pre$Species)
JF2001_pre$Species <- gsub("einnovatus", "Elymus innovatus", JF2001_pre$Species)
JF2001_pre$Species <- gsub("fcamp", "Festuca campestris", JF2001_pre$Species)
JF2001_pre$Species <- gsub("ppratense", "Phleum pratense", JF2001_pre$Species)
# unique(JF2004_pre$Species)
JF2004_pre <- JF2004_pre %>% # I'll use stringr's str_replace_all for this one because there are so many!
  mutate(Species=stringr::str_replace_all(Species, c("A.artemisi"="Ambrosia artemisiifolia",
                                               "C.album"="Chenopodium album",
                                               "C.communis"="Carex communis",
                                               "C.hirtifol"="Carex hirtifolia",
                                               "C.peduncul"="Carex pedunculata",
                                               "C.plantagi"="Carex plantaginea",
                                               "C.stipata"="Carex stipata",
                                               "P.lanceola"="Plantago lanceolata",
                                               "R.acetosel"="Rumex acetosella",
                                               "R.crispus"="Rumex crispus",
                                               "S.microcar"="Scirpus microcarpus",
                                               "S.purp"="Schizacne purpurascens",
                                               "T.dioicum"="Thalictrum dioicum" ))) %>% 
  distinct() # make sure there are no duplicate rows

# nrow(JF2004_pre)
# Combine all the data together into one complete dataset, stig-all.csv
stig_all <- rbind(rbind(CS2021_pre, JF2001_pre), JF2004_pre)
# Add sex system info
stig_all_ss <- stig_all %>% 
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
# Make sure no species were missed -- looks ok now
# summary(as.factor(stig_all_ss$Sex_sys))
# stig_all_ss[which(is.na(stig_all_ss$Sex_sys)),]
write.csv(stig_all_ss, "processed-data/stig-all.csv", row.names=F)


# I also want a version that has no repeat species. Remove repeats of mine from JF data - except Chenopodium album. Replace
# the CS2021 C. album with the JF2004 C. album. I collected mine too early. 
CS2021_norep <- CS2021_pre %>% filter(Species != "Chenopodium album")
unique(CS2021_norep$Species)
JF2001_norep <- JF2001_pre %>% filter(!(Species %in% unique(CS2021_norep$Species)))
unique(JF2001_norep$Species)
JF2004_norep <- JF2004_pre %>% filter(!(Species %in% unique(CS2021_norep$Species)))
unique(JF2004_norep$Species)
stig_norep <- rbind(rbind(CS2021_norep, JF2001_norep), JF2004_norep)

# Add sex system info
stig_norep_ss <- stig_norep %>% 
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
# summary(as.factor(stig_norep_ss$Sex_sys))
# stig_all_ss[which(is.na(stig_all_ss$Sex_sys)),]

write.csv(stig_norep_ss, "processed-data/stig-no-rep-spp.csv", row.names=F)

