#### Wheat Metabolomics Project

### Data
`prepared_ms`- filtered, normalized mass spec DF using `MSPrep` package in R.

I used the following parameters to filter and normalize the data:

prepared_ms <- msPrepare(ms_formulae_prep, 
                            cvMax = 0.50,
                            minPropPresent = 2/3,
                            compVars = c("Calc..MW", "RT..min."),
                            sampleVars = c("subject_id","replicate"),
                            colExtraText = "Soil.raw.._Peak.Rating..20250318_",
                            separator = "-",
                            missingValue = 0, # placeholder. I have loads though (how many? I should check)
                            imputeMethod = "bpca", # what is this?
                            nPcs = 3,
                            normalizeMethod = "median",
                            transform = "none"
                            )


### Results
pathways.out- 'pathway intensities'. Higher pathway abundance corresponds to greater gene family abundance (more copies of the genes involved in that pathway) in that particular sample, which reflects a greater relative metabolic activity 