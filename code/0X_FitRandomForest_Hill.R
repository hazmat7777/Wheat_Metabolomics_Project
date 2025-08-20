# Script to run random forest modelling on FutureReef Hill diversity with anthropogenic stress data 

# Packages ----------------------------------
library(tidyverse)
options(device = "quartz")
library(ranger)
#install.packages("mlr")
library(mlr)
install.packages("tuneRanger")
library(tuneRanger)

# Functions ----------------------------------
createRarefiedReplicateList<-function(data){
    replicates<-unique(data$rarefaction_replicate)
    replicate_l<-list()
    for (r in replicates) {
        replicate_l[[r]]<-filter(data, rarefaction_replicate==r) %>%
            dplyr::select(-"rarefaction_replicate")
    }
    return(replicate_l)
}

getRcvSplit <-function(data, testProportion){
    n<-dim(data)[1]
    indices<-sample(n, testProportion*n )
    return(list(
            "train"=data[-indices,],
            "test"=data[indices,]))
}

rsq1 <- function (x, y) cor(x, y) ^ 2
rsq2 <- function(preds, actual) {
    rss <- sum((preds - actual) ^ 2)
    tss <- sum((actual - mean(actual)) ^ 2)
    rsq <- 1 - rss/tss
    return(rsq)
}

# 0
df_rf_COI_td_hill_0<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X0", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))


# 0.2
df_rf_COI_td_hill_0.2<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X0.2", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))


# 0.4
df_rf_COI_td_hill_0.4<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X0.4", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 0.6
df_rf_COI_td_hill_0.6<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X0.6", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 0.8
df_rf_COI_td_hill_0.8<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X0.8", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 1
df_rf_COI_td_hill_1<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X1", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 1.2
df_rf_COI_td_hill_1.2<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X1.2", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 1.4
df_rf_COI_td_hill_1.4<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X1.4", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 1.6
df_rf_COI_td_hill_1.6<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X1.6", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 1.8
df_rf_COI_td_hill_1.8<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X1.8", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 2
df_rf_COI_td_hill_2<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X2", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 2.2
df_rf_COI_td_hill_2.2<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X2.2", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 2.4
df_rf_COI_td_hill_2.4<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X2.4", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 2.6
df_rf_COI_td_hill_2.6<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X2.6", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 2.8
df_rf_COI_td_hill_2.8<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X2.8", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 3
df_rf_COI_td_hill_3<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X3", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 3.2
df_rf_COI_td_hill_3.2<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X3.2", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 3.4
df_rf_COI_td_hill_3.4<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X3.4", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 3.6
df_rf_COI_td_hill_3.4<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X3.6", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 3.8
df_rf_COI_td_hill_3.4<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X3.8", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# 4
df_rf_COI_td_hill_3.4<-read.csv("Outputs/RF/COI_hilldiv_stress_PCR_touchdown_blast_query85.csv") %>%
    as_tibble() %>%
    dplyr::select(
        "lengthDeployment", "depthOfBottomInMeters", "decimalLatitude",
        "X4", "proportionshallow", "distance.km.", "grav_NC_raw",  "sediment_raw" , "pop_count_raw" , "reef_value_raw",
        "cumul_score", "region",
        "scorecn", "scoreth", "scoretr", "score", "scorecy",
        "global_effluent_2015_open_N", "global_effluent_2015_septic_N", "global_effluent_2015_treated_N",
        "mean.ARMS.variable.mean_DHW","mean.ARMS.variable.mean_SST","mean.ARMS.variable.mean_SST.anomaly",    
        "mean.ARMS.variable.mean_SST.variability","mean.ARMS.variable.mean_cloud", "mean.ARMS.variable.mean_wind", 
        "rarefaction_replicate") %>%
    filter(complete.cases(.))

# Set parameters ----------------------------------
tuneRandomForest<-FALSE

# Structure data into a list of dataframes ----------------------------------
repDf_COI_td_l_0<-createRarefiedReplicateList(df_rf_COI_td_hill_0)
repDf_COI_td_l_0.2<-createRarefiedReplicateList(df_rf_COI_td_hill_0.2)
repDf_COI_td_l_0.4<-createRarefiedReplicateList(df_rf_COI_td_hill_0.4)
repDf_COI_td_l_0.6<-createRarefiedReplicateList(df_rf_COI_td_hill_0.6)
repDf_COI_td_l_0.8<-createRarefiedReplicateList(df_rf_COI_td_hill_0.8)
repDf_COI_td_l_1<-createRarefiedReplicateList(df_rf_COI_td_hill_1)
repDf_COI_td_l_1.2<-createRarefiedReplicateList(df_rf_COI_td_hill_1.2)
repDf_COI_td_l_1.4<-createRarefiedReplicateList(df_rf_COI_td_hill_1.4)
repDf_COI_td_l_1.6<-createRarefiedReplicateList(df_rf_COI_td_hill_1.6)
repDf_COI_td_l_1.8<-createRarefiedReplicateList(df_rf_COI_td_hill_1.8)
repDf_COI_td_l_2<-createRarefiedReplicateList(df_rf_COI_td_hill_2)
repDf_COI_td_l_2.2<-createRarefiedReplicateList(df_rf_COI_td_hill_2.2)
repDf_COI_td_l_2.4<-createRarefiedReplicateList(df_rf_COI_td_hill_2.4)
repDf_COI_td_l_2.6<-createRarefiedReplicateList(df_rf_COI_td_hill_2.6)
repDf_COI_td_l_2.8<-createRarefiedReplicateList(df_rf_COI_td_hill_2.8)
repDf_COI_td_l_3<-createRarefiedReplicateList(df_rf_COI_td_hill_3)
repDf_COI_td_l_3.2<-createRarefiedReplicateList(df_rf_COI_td_hill_3.2)
repDf_COI_td_l_3.4<-createRarefiedReplicateList(df_rf_COI_td_hill_3.4)
repDf_COI_td_l_3.6<-createRarefiedReplicateList(df_rf_COI_td_hill_3.6)
repDf_COI_td_l_3.8<-createRarefiedReplicateList(df_rf_COI_td_hill_3.8)
repDf_COI_td_l_4<-createRarefiedReplicateList(df_rf_COI_td_hill_4)

# Q=0
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_0[[1]], target = "X0")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)

    # Mean of best 5 % of the results
    # mtry=5
    # min.nod.size=29
    # sample.fraction=ca.0.2
    res
}

# Q=0.2
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_0.2[[1]], target = "X0.2")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=0.4
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_0.4[[1]], target = "X0.4")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=0.6
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_0.6[[1]], target = "X0.6")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=0.8
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_0.8[[1]], target = "X0.8")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=1
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_1[[1]], target = "X1")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=1.2
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_1.2[[1]], target = "X1.2")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=1.4
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_1.4[[1]], target = "X1.4")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=1.6
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_1.6[[1]], target = "X1.6")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=1.8
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_1.8[[1]], target = "X1.8")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=2
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_2[[1]], target = "X2")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=2.2
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_2.2[[1]], target = "X2.2")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=2.4
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_2.4[[1]], target = "X2.4")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=2.6
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_2.6[[1]], target = "X2.6")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=2.8
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_2.8[[1]], target = "X2.8")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=3
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_3[[1]], target = "X3")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=3.2
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_3.2[[1]], target = "X3.2")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=3.4
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_3.4[[1]], target = "X3.4")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=3.6
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_3.6[[1]], target = "X3.6")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=3.8
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_3.8[[1]], target = "X3.8")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

# Q=4
if (tuneRandomForest) {
    task = mlr::makeClassifTask(data = repDf_COI_td_l_4[[1]], target = "X4")
    tuneRanger::estimateTimeTuneRanger(task)
    res = tuneRanger::tuneRanger(
        task, 
        measure = list(multiclass.brier), 
        num.trees = 5001,
        num.threads = 4, 
        iters = 70, 
        save.file.path = NULL)
    res
}

### getRcvSplit ----------------------------------

rcvSplitDf_COI_td_l_0<-lapply(repDf_COI_td_l_0, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_0.2<-lapply(repDf_COI_td_l_0.2, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_0.4<-lapply(repDf_COI_td_l_0.4, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_0.6<-lapply(repDf_COI_td_l_0.6, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_0.8<-lapply(repDf_COI_td_l_0.8, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_1<-lapply(repDf_COI_td_l_1, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_1.2<-lapply(repDf_COI_td_l_1.2, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_1.4<-lapply(repDf_COI_td_l_1.4, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_1.6<-lapply(repDf_COI_td_l_1.6, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_1.8<-lapply(repDf_COI_td_l_1.8, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_2<-lapply(repDf_COI_td_l_2, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_2.2<-lapply(repDf_COI_td_l_2.2, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_2.4<-lapply(repDf_COI_td_l_2.4, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_2.6<-lapply(repDf_COI_td_l_2.6, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_2.8<-lapply(repDf_COI_td_l_2.8, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_3<-lapply(repDf_COI_td_l_3, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_3.2<-lapply(repDf_COI_td_l_3.2, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_3.4<-lapply(repDf_COI_td_l_3.4, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_3.6<-lapply(repDf_COI_td_l_3.6, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_3.8<-lapply(repDf_COI_td_l_3.8, getRcvSplit, 0.2)
rcvSplitDf_COI_td_l_4<-lapply(repDf_COI_td_l_4, getRcvSplit, 0.2)

#################################################################
### Fit RCV random forests ### ----------------------------------

# Q = 0, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_0<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_0)){
    rf<-ranger(X0~.,
    data = rcvSplitDf_COI_td_l_0[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_0[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_0[[i]][["test"]]$X0)
    varImp_COI_td_l_0[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=1200, y=500) +
  labs(x="Predicted Q = 0", y="Q = 0")+
  theme_classic()+
  ggtitle("Q = 0 (touchdown PCR only)")

varImp_COI_td_df_0<-bind_rows(varImp_COI_td_l_0) %>%
    pivot_longer(cols=everything(), values_to="Importance", names_to="Variable")

varImp_COI_td_df_0 %>%
    group_by(Variable) %>%
    reframe(Importance=Importance, meanImp=mean(Importance)) %>%
    mutate(Variable = fct_reorder(Variable, desc(-meanImp))) %>%
    ggplot(data=., aes(x=Importance, y=Variable))+
        geom_boxplot()+
        theme_classic()


# Q = 0.2, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_0.2<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_0.2)){
    rf<-ranger(X0.2~.,
    data = rcvSplitDf_COI_td_l_0.2[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_0.2[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_0.2[[i]][["test"]]$X0.2)
    varImp_COI_td_l_0.2[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=1200, y=500) +
  labs(x="Predicted Q = 0.2", y="Q = 0.2")+
  theme_classic()+
  ggtitle("Q = 0.2 (touchdown PCR only)")


# Q = 0.4, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_0.4<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_0.4)){
    rf<-ranger(X0.4~.,
    data = rcvSplitDf_COI_td_l_0.4[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_0.4[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_0.4[[i]][["test"]]$X0.4)
    varImp_COI_td_l_0.4[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=1200, y=500) +
  labs(x="Predicted Q = 0.4", y="Q = 0.4")+
  theme_classic()+
  ggtitle("Q = 0.4 (touchdown PCR only)")



# Q = 0.6, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_0.6<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_0.6)){
    rf<-ranger(X0.6~.,
    data = rcvSplitDf_COI_td_l_0.6[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_0.6[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_0.6[[i]][["test"]]$X0.6)
    varImp_COI_td_l_0.4[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=1200, y=500) +
  labs(x="Predicted Q = 0.6", y="Q = 0.6")+
  theme_classic()+
  ggtitle("Q = 0.6 (touchdown PCR only)")


# Q = 0.8, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_0.8<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_0.8)){
    rf<-ranger(X0.8~.,
    data = rcvSplitDf_COI_td_l_0.8[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_0.8[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_0.8[[i]][["test"]]$X0.8)
    varImp_COI_td_l_0.8[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=1200, y=500) +
  labs(x="Predicted Q = 0.8", y="Q = 0.8")+
  theme_classic()+
  ggtitle("Q = 0.8 (touchdown PCR only)")


# Q = 1, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_1<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_1)){
    rf<-ranger(X1~.,
    data = rcvSplitDf_COI_td_l_1[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_1[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_1[[i]][["test"]]$X1)
    varImp_COI_td_l_1[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=150, y=300) +
  labs(x="Predicted Q = 1", y="Q = 1")+
  theme_classic()+
  ggtitle("Q = 1 (touchdown PCR only)")

#varImp_COI_td_df_1<-bind_rows(varImp_COI_td_l_1) %>%
#    pivot_longer(cols=everything(), values_to="Importance", names_to="Variable")
#varImp_COI_td_df_1 %>%
#    group_by(Variable) %>%
#    reframe(Importance=Importance, meanImp=mean(Importance)) %>%
#    mutate(Variable = fct_reorder(Variable, desc(-meanImp))) %>%
#    ggplot(data=., aes(x=Importance, y=Variable))+
#        geom_boxplot()+
#        theme_classic()


# Q = 1.2, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_1.2<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_1.2)){
    rf<-ranger(X1.2~.,
    data = rcvSplitDf_COI_td_l_1.2[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_1.2[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_1.2[[i]][["test"]]$X1.2)
    varImp_COI_td_l_1.2[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=150, y=300) +
  labs(x="Predicted Q = 1.2", y="Q = 1.2")+
  theme_classic()+
  ggtitle("CQ = 1.2 (touchdown PCR only)")


# Q = 1.4, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_1.4<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_1.4)){
    rf<-ranger(X1.4~.,
    data = rcvSplitDf_COI_td_l_1.4[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_1.4[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_1.4[[i]][["test"]]$X1.4)
    varImp_COI_td_l_1.4[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=150, y=300) +
  labs(x="Predicted Q = 1.4", y="Q = 1.4")+
  theme_classic()+
  ggtitle("CQ = 1.4 (touchdown PCR only)")


# Q = 1.6, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_1.6<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_1.6)){
    rf<-ranger(X1.6~.,
    data = rcvSplitDf_COI_td_l_1.6[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_1.6[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_1.6[[i]][["test"]]$X1.6)
    varImp_COI_td_l_1.6[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=150, y=300) +
  labs(x="Predicted Q = 1.6", y="Q = 1.6")+
  theme_classic()+
  ggtitle("CQ = 1.6 (touchdown PCR only)")


# Q = 1.8, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_1.8<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_1.8)){
    rf<-ranger(X1.8~.,
    data = rcvSplitDf_COI_td_l_1.8[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_1.8[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_1.8[[i]][["test"]]$X1.8)
    varImp_COI_td_l_1.8[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=150, y=300) +
  labs(x="Predicted Q = 1.8", y="Q = 1.8")+
  theme_classic()+
  ggtitle("CQ = 1.8 (touchdown PCR only)")


# Q = 2, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_2<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_2)){
    rf<-ranger(X2~.,
    data = rcvSplitDf_COI_td_l_2[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_2[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_2[[i]][["test"]]$X2)
    varImp_COI_td_l_2[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2)), x=150, y=300) +
  labs(x="Predicted Q = 2", y="Q = 2")+
  theme_classic()+
  ggtitle("CQ = 2 (touchdown PCR only)")



# Q = 2.2, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_2.2<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_2.2)){
    rf<-ranger(X2.2~.,
    data = rcvSplitDf_COI_td_l_2.2[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=2.29,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_2.2[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_2.2[[i]][["test"]]$X2.2)
    varImp_COI_td_l_2.2[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2.2)), x=150, y=300) +
  labs(x="Predicted Q = 2.2", y="Q = 2.2")+
  theme_classic()+
  ggtitle("CQ = 2.2 (touchdown PCR only)")


  # Q = 2.4, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_2.4<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_2.4)){
    rf<-ranger(X2.4~.,
    data = rcvSplitDf_COI_td_l_2.4[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=2.49,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_2.4[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_2.4[[i]][["test"]]$X2.4)
    varImp_COI_td_l_2.4[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2.4)), x=150, y=300) +
  labs(x="Predicted Q = 2.4", y="Q = 2.4")+
  theme_classic()+
  ggtitle("CQ = 2.4 (touchdown PCR only)")


  # Q = 2.6, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_2.6<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_2.6)){
    rf<-ranger(X2.6~.,
    data = rcvSplitDf_COI_td_l_2.6[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=2.69,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_2.6[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_2.6[[i]][["test"]]$X2.6)
    varImp_COI_td_l_2.6[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2.6)), x=150, y=300) +
  labs(x="Predicted Q = 2.6", y="Q = 2.6")+
  theme_classic()+
  ggtitle("CQ = 2.6 (touchdown PCR only)")


  # Q = 2.8, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_2.8<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_2.8)){
    rf<-ranger(X2.8~.,
    data = rcvSplitDf_COI_td_l_2.8[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=2.89,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_2.8[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_2.8[[i]][["test"]]$X2.8)
    varImp_COI_td_l_2.8[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 2.8)), x=50, y=50) +
  labs(x="Predicted Q = 2.8", y="Q = 2.8")+
  theme_classic()+
  ggtitle("CQ = 2.8 (touchdown PCR only)")


  # Q = 3, touchdown blast query 85 ----------------------------------
predictions<-c()
observations<-c()

varImp_COI_td_l_3<-list()
for (i in 1:length(rcvSplitDf_COI_td_l_3)){
    rf<-ranger(X3~.,
    data = rcvSplitDf_COI_td_l_3[[i]][["train"]], 
    importance = 'permutation',
    local.importance = TRUE,
    scale.permutation.importance = TRUE,
    mtry = 5,
    min.node.size=39,
    num.trees=5001)

    p_rf<-predict(rf, rcvSplitDf_COI_td_l_3[[i]][["test"]])
    predictions<-c(predictions, p_rf$predictions)
    observations<-c(observations, rcvSplitDf_COI_td_l_3[[i]][["test"]]$X3)
    varImp_COI_td_l_3[[i]]<-rf$variable.importance
}

po_COI_td_df<-data.frame("preds"=predictions, "obs"=observations)
ggplot(po_COI_td_df, aes(x=preds, y=obs))+
  geom_point()+
  geom_smooth(method="lm")+
  annotate("text", label=paste0("R2=",round(rsq2(po_COI_td_df$preds, po_COI_td_df$obs), 3)), x=150, y=300) +
  labs(x="Predicted Q = 3", y="Q = 3")+
  theme_classic()+
  ggtitle("CQ = 3 (touchdown PCR only)")
