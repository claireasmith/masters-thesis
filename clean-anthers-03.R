## Combine CS size and number data 
## Claire Smith
## Last updated: 26 July 2023

# Create size-prod-CS2021.csv, which combines size-CS2021.csv and prod-CS2021.csv.
# Each row is an individual, containing its average pollen diam (in um) and 
# average pollen count per anther. 

############################################################################################
# Load libraries
library(tidyverse) # need dplyr and tidyr specifically

# Be sure to have run clean-anthers-01.R and clean-anthers-02.R before running this to make 
# sure files are up to date
############################################################################################
# Load data
prod <- read.csv("processed-data/prod-CS2021.csv", stringsAsFactors = T)
head(prod)
summary(prod)
size <- read.csv("processed-data/size-CS2021.csv", stringsAsFactors = T)
head(size)
summary(size)

# Take a look
# View(prod)
# View(size)
# In production data, "label" is the same as "ind" but sometimes will have letters A or M corresponding to "male" or
# "anthers" that's an artifact of how I labelled tubes when I sampled. Remove these letters so that the format
# matches the "ind" column in the size data. 
# unique(prod$Label)
prod <- prod %>% mutate(Ind = gsub("\\D", "", Label))
# unique(prod$Ind) # looks good!


# Join together the data
sizeprod <- full_join(size, prod, by=c("source", "Sex_sys", "Species", "Site", "Ind"))
# Data not complete for Setaria viridis (missing individual labels from size data)
sizeprod <- sizeprod %>% filter(Species != "Setaria viridis")

# Take out incomplete rows
sizeprod <- sizeprod %>% 
  filter(!is.na(Avg_area)&!is.na(Avg_pol_anth))

# Take a look at sample size per species: 
# sizeprod %>% group_by(Species) %>% 
#   summarize(n=n())
# A tibble: 6 × 2
# Species                     n
# <fct>                   <int>
#   1 Amaranthus retroflexus     16
# 2 Ambrosia artemisiifolia    10
# 3 Chenopodium album          19
# 4 Plantago lanceolata        28
# 5 Rumex acetosella            9
# 6 Thalictrum dioicum          7

# summary(sizeprod)

write.csv(sizeprod, "processed-data/size-prod-CS2021.csv", row.names=F)
