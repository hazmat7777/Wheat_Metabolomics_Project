library(tidyverse)

metadata <- read.csv("../data/EnvData_md.csv")

glimpse(metadata)

envdata = read.csv("../data/GreenMicrobiome_FarmKits_EnvData.csv")

glimpse(envdata)

colnames(envdata)

for(f in colnames(envdata)){
    print("*******NEW_COLUMN*********")
    print(f)
    print(unique(envdata[, f]))
}

## inconsistencies

# missing values

    # oats, winter_oats- not mer
    # wheat_group2, 2nd_wheat
    # maize, maize_forage - not mergeable
    # feed_barley, winter_barley, spring_barley, spring_malting_barley

# duplicates?
    # oilseed_rape and OSR- mergeable
    # wheat, winter_wheat, seed_wheat- mergeable?
    # vining_peas, combining_peas- - not mergeable
# not sure
    # pre_pack- prepackaged things
    # fallow_MLS
    # group_1/2/3

## problems
    # crop_rotation_companion y2-4 are NA


# WHEAT CULTIVAR
# [1] "wheat_cultivar"
#  [1] ""                                        
#  [2] "Common wheat/bread wheat"                
#  [3] "Crusoe; extase"                          
#  [4] "Extase; Astronomer"                      
#  [5] "Crusoe; Extase"                          
#  [6] "Common wheat/bread wheat; Heritage wheat"
#  [7] "Extase; extase"                          
#  [8] "Champion; extase"                        
#  [9] "X days"                                  
# [10] "Crusoe; skyscraper; firefly"             
# [11] "Extase; Cranium"                         
# [12] "Crusoe; skyscraper; "                    
# [13] "Crusoe; gravity; crusoe"                 
# [14] "Var Edgar"                               
# [15] "Var Palladium "                          
# [16] "Extase; astronomer" 


# fertilisation rate useless:
# Variable"                                                     
# [3] "Nitogen 3 times in the Spring (March, April, May)"            
# [4] "As needed"                                                    
# [5] "Per crop requirements"                                        
# [6] "220Kile N/hectare; 20 kilo S/hectare"                         
# [7] "3 sp. Applications totalling 190kgN/ha"                       
# [8] "240kg/Ha split in three"                                      
# [9] "Nitrogen- March, April, May; Phosphate-March; Potassium-April"


# ?
# [1] "tillage_method"
# [1] ""                          "Reduced tillage"          
# [3] "Minimum tillage"           "Ploughed, Primary Tilling"
# [5] "Primary tilling"           "Strip tillage"            
# [7] "Ploughed"                  "Direct drilling "    

    # merge min and reduced tillage
    # merge ploughed and primary tilling and both

### cleaning data
colnames(envdata)

envdata[envdata == "OSR"] <- "oilseed_rape"

unique(envdata[, "crop_rotation_primary_y1"])


# # remove rows for which values are empty
    # think I can't do this bc I need to merge it to the other phyloseq object
# envdata2 <- envdata  %>% 
#     filter(site_name != "") # use sitename being empty as an indicator (see excel)

# dim(envdata2) # lost 10 rows




# probably need to make strings factors
