sort(colnames(boost_data))[400:500]
colnames(boost_data)[grep("Pyoluteorin", colnames(boost_data))]

# metabs of interest
#1-(4-hydroxyphenyl)ethanone, name "Massbank.LU021452.4.Hydroxyacetophenone.4..Hydroxyacetophenone.1..4.hydroxyphenyl.ethanone"
#"Phloroglucinol"
#"phenazine-1-carboxylic acid"- not found
#  pyoluteorin - not found

metabolites <- c("Massbank.LU021452.4.Hydroxyacetophenone.4..Hydroxyacetophenone.1..4.hydroxyphenyl.ethanone", "Phloroglucinol")

# new fn to plot boxplots for metabolites of interest against GT
plot_metabolites_vs_gt <- function(metabolites, data, target_col = "gt"){
    require(ggplot2)
    
    for(metab in metabolites) {
        if(metab %in% colnames(data)) {
        plot <- ggplot(data, aes_string(x = target_col, y = metab)) +
            geom_boxplot() +
            labs(x = "GT Status", y = paste("Distribution of", metab)) +
            theme_bw()
        
        ggsave(filename = paste0("../results/plots/metab_", metab, "_vs_gt.png"), plot = plot, width = 8, height = 6)
        } else {
        message(paste("Metabolite", metab, "not found in data. Skipping."))
        }
    }
}

# Call the function with the metabolites of interest
plot_metabolites_vs_gt(metabolites, fk_metabolom_gt_t)

# do an anova for each metabolite against GT
anova_results <- lapply(metabolites, function(metab) {
    if(metab %in% colnames(fk_metabolom_gt_t)) {
        print(paste("Running ANOVA for:", metab))
        model <- aov(as.formula(paste(metab, "~ gt")), data = fk_metabolom_gt_t)
        return(summary(model))
    } else {
        message(paste("Metabolite", metab, "not found in data. Skipping."))
        return(NULL)
    }
})

anova_results