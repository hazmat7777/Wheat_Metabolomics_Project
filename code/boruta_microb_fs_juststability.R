# Clean Microbiome Stability Results with Improved Taxonomic Names
library(phyloseq)

# Load your existing stability results and phyloseq object
microbiome_stability_results <- read.csv("../results/microbiome_stability_results.csv", stringsAsFactors = FALSE)
ps_genus <- readRDS("../data/metabarcoding/ps_16S_highdiv_genus_relative.rds")

cat("Loaded stability results with", nrow(microbiome_stability_results), "features\n")

# Function to improve taxonomic names using phyloseq tax_table
improve_microbiome_names <- function(feature_names, ps_object) {
  improved_names <- character(length(feature_names))
  tax_tab <- tax_table(ps_object)
  
  for(i in seq_along(feature_names)) {
    feature <- feature_names[i]
    
    # Check if feature exists in tax table (by taxa_names or by matching)
    if(feature %in% taxa_names(ps_object)) {
      tax_info <- tax_tab[feature, ]
      
      # Extract taxonomic levels
      genus <- as.character(tax_info[,"Genus"])
      family <- as.character(tax_info[,"Family"])
      order <- as.character(tax_info[,"Order"])
      class <- as.character(tax_info[,"Class"])
      
      # Check if genus starts with "Incertae" or is NA/empty
      if(!is.na(genus) && grepl("^Incertae", genus)) {
        # Use family name with "Unclassified" prefix
        if(!is.na(family) && family != "" && !grepl("^Incertae", family)) {
          improved_names[i] <- paste("Unclassified", family)
        } else if(!is.na(order) && order != "" && !grepl("^Incertae", order)) {
          improved_names[i] <- paste("Unclassified", order)
        } else if(!is.na(class) && class != "" && !grepl("^Incertae", class)) {
          improved_names[i] <- paste("Unclassified", class)
        } else {
          improved_names[i] <- "Unclassified bacterium"
        }
      } else if(is.na(genus) || genus == "" || genus == "NA") {
        # No genus info, use family
        if(!is.na(family) && family != "" && !grepl("^Incertae", family)) {
          improved_names[i] <- paste("Unclassified", family)
        } else if(!is.na(order) && order != "" && !grepl("^Incertae", order)) {
          improved_names[i] <- paste("Unclassified", order)
        } else {
          improved_names[i] <- "Unclassified bacterium"
        }
      } else {
        # Use genus name, clean it up
        clean_genus <- genus
        # Remove common prefixes/suffixes that make names messy
        clean_genus <- gsub("\\[|\\]", "", clean_genus)  # Remove brackets
        clean_genus <- gsub("_.*", "", clean_genus)       # Remove everything after underscore
        clean_genus <- trimws(clean_genus)
        improved_names[i] <- clean_genus
      }
    } else {
      # Feature not found in tax table, try to clean the original name
      clean_name <- feature
      clean_name <- gsub("\\.", " ", clean_name)        # Replace dots with spaces
      clean_name <- gsub("^X", "", clean_name)          # Remove leading X
      clean_name <- gsub("  +", " ", clean_name)        # Remove multiple spaces
      clean_name <- trimws(clean_name)                  # Remove leading/trailing spaces
      improved_names[i] <- clean_name
    }
  }
  
  return(improved_names)
}

# Get feature names from the Original_Feature_Name column
original_features <- microbiome_stability_results$Original_Feature_Name

# Improve the names
cat("Improving taxonomic names...\n")
improved_names <- improve_microbiome_names(original_features, ps_genus)

# Create updated stability export with improved names
microbiome_stability_improved <- data.frame(
  Improved_Taxonomic_Name = improved_names,
  Original_Feature_Name = microbiome_stability_results$Original_Feature_Name,
  Times_Selected = microbiome_stability_results$Times_Selected,
  Selection_Rate = microbiome_stability_results$Selection_Rate,
  stringsAsFactors = FALSE
)

# Show examples of improvements
cat("\n=== EXAMPLES OF NAME IMPROVEMENTS ===\n")
# Find examples where names changed
changed_indices <- which(microbiome_stability_improved$Improved_Taxonomic_Name != microbiome_stability_improved$Original_Feature_Name)
if(length(changed_indices) > 0) {
  n_examples <- min(10, length(changed_indices))
  for(i in 1:n_examples) {
    idx <- changed_indices[i]
    cat("Original:", microbiome_stability_improved$Original_Feature_Name[idx], 
        "\n-> Improved:", microbiome_stability_improved$Improved_Taxonomic_Name[idx], "\n\n")
  }
} else {
  cat("No significant name changes found\n")
}

# Look for "Incertae" cases specifically
incertae_examples <- grep("Unclassified", microbiome_stability_improved$Improved_Taxonomic_Name)
if(length(incertae_examples) > 0) {
  cat("Examples of 'Incertae Sedis' → 'Unclassified Family' conversions:\n")
  for(i in head(incertae_examples, 5)) {
    cat("→", microbiome_stability_improved$Improved_Taxonomic_Name[i], "\n")
  }
}

# Export improved results
output_file <- "../results/microbiome_stability_results_improved_names.csv"
write.csv(microbiome_stability_improved, output_file, row.names = FALSE)
cat("\nImproved stability table exported to:", output_file, "\n")

# Show summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Total features:", nrow(microbiome_stability_improved), "\n")
cat("Features with improved names:", sum(microbiome_stability_improved$Improved_Taxonomic_Name != microbiome_stability_improved$Original_Feature_Name), "\n")
cat("Features with 'Unclassified' prefix:", sum(grepl("^Unclassified", microbiome_stability_improved$Improved_Taxonomic_Name)), "\n")
cat("Features selected ≥50% of time:", sum(microbiome_stability_improved$Selection_Rate >= 0.5), "\n")
cat("Features selected ≥75% of time:", sum(microbiome_stability_improved$Selection_Rate >= 0.75), "\n")

# Show top 10 stable features with improved names
cat("\n=== TOP 10 MOST STABLE FEATURES (IMPROVED NAMES) ===\n")
top10_improved <- head(microbiome_stability_improved[, c("Improved_Taxonomic_Name", "Times_Selected", "Selection_Rate")], 10)
print(top10_improved)

# Create a summary of taxonomic name types
cat("\n=== TAXONOMIC NAME TYPE SUMMARY ===\n")
name_types <- data.frame(
  Category = c("Proper genus names", "Unclassified (Family level)", "Unclassified (Order level)", 
               "Unclassified (Class level)", "Unclassified bacterium", "Other"),
  Count = c(
    sum(!grepl("^Unclassified", microbiome_stability_improved$Improved_Taxonomic_Name)),
    sum(grepl("^Unclassified [A-Z][a-z]+aceae$", microbiome_stability_improved$Improved_Taxonomic_Name)),
    sum(grepl("^Unclassified [A-Z][a-z]+ales$", microbiome_stability_improved$Improved_Taxonomic_Name)),
    sum(grepl("^Unclassified [A-Z][a-z]+ia$", microbiome_stability_improved$Improved_Taxonomic_Name)),
    sum(microbiome_stability_improved$Improved_Taxonomic_Name == "Unclassified bacterium"),
    sum(grepl("^Unclassified", microbiome_stability_improved$Improved_Taxonomic_Name) & 
        !grepl("aceae$|ales$|ia$|bacterium$", microbiome_stability_improved$Improved_Taxonomic_Name))
  )
)
print(name_types)

cat("\nScript completed successfully!\n")
cat("Improved names file:", output_file, "\n")