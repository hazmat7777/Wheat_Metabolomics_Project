df <- read.csv("../data/20250327_Soil_Bact_OB.csv")

# this is the practice metabolomics dataset from 10 (?) samples
# using msdial

ncol(df)
colnames(df)
?order()

library(dplyr)
glimpse(df)

length(unique(df$Formula))

dim(df2) # 225 columns, 37187 samples
length(unique(df$Calc..MW)) # all unique masses

df2<- df[with(df, order(Calc..MW)), ] # ordered by MW

df2$Calc..MW[1]
min(df$Calc..MW) # it worked

