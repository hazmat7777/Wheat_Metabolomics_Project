if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("MAIT")


browseVignettes("MAIT")