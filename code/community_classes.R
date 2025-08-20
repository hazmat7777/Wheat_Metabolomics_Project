# Save this as a new clean file
library(tidyverse)
library(phyloseq)
#install.packages("xgboost")
library(xgboost)
library(caret)

ps_highdiv_relative <- readRDS("../data/metabarcoding/ps_16S_highdiv_relative.rds")
ps_genus <- tax_glom(ps_highdiv_relative, "Genus", NArm = TRUE)
taxa_names(ps_genus) <- make.unique(as.character(tax_table(ps_genus)[, "Genus"]))

cat("ESVs after tax_glom:", ntaxa(ps_genus), "\n")

#####

# Load required libraries
library(phyloseq)
library(cluster)
library(vegan)
library(ggplot2)
library(dplyr)
if (!require(fpc, quietly = TRUE)) {
  install.packages("fpc")
  library(fpc)
}
# Function to calculate Jensen-Shannon divergence
jensen_shannon_divergence <- function(p, q) {
  # Ensure no zeros (add small pseudocount)
  p <- p + 1e-10
  q <- q + 1e-10
  
  # Normalize to sum to 1
  p <- p / sum(p)
  q <- q / sum(q)
  
  # Calculate M (average distribution)
  m <- (p + q) / 2
  
  # Calculate JS divergence
  js <- 0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m))
  
  return(sqrt(js))  # Return JS distance (square root of divergence)
}

# Function to compute all-against-all JS divergence matrix
compute_js_matrix <- function(otu_table) {
  n_samples <- nrow(otu_table)
  js_matrix <- matrix(0, nrow = n_samples, ncol = n_samples)
  rownames(js_matrix) <- rownames(otu_table)
  colnames(js_matrix) <- rownames(otu_table)
  
  cat("Computing JS divergence matrix for", n_samples, "samples...\n")
  
  for (i in 1:n_samples) {
    if (i %% 50 == 0) cat("Progress:", i, "/", n_samples, "\n")
    
    for (j in i:n_samples) {
      if (i == j) {
        js_matrix[i, j] <- 0
      } else {
        js_dist <- jensen_shannon_divergence(otu_table[i, ], otu_table[j, ])
        js_matrix[i, j] <- js_dist
        js_matrix[j, i] <- js_dist  # Symmetric matrix
      }
    }
  }
  
  return(js_matrix)
}

# Function to compute Calinski-Harabasz index
calinski_harabasz <- function(distance_matrix, clusters) {
  # Convert to distance object if not already
  if (!inherits(distance_matrix, "dist")) {
    distance_matrix <- as.dist(distance_matrix)
  }
  
  # Use cluster.stats from fpc package or calculate manually
  if (require(fpc, quietly = TRUE)) {
    stats <- cluster.stats(distance_matrix, clusters)
    return(stats$ch)
  } else {
    # Manual calculation
    n <- length(clusters)
    k <- length(unique(clusters))
    
    if (k == 1 || k == n) return(0)
    
    # This is a simplified version - install fpc package for full calculation
    coords <- cmdscale(distance_matrix, k = min(n-1, 10))
    ch <- calinhara(coords, clusters)
    return(ch)
  }
}

# Main analysis function
identify_community_classes <- function(ps_genus, k_range = 2:10) {
  
  # Extract OTU table (samples x taxa)
  otu_table <- as.data.frame(otu_table(ps_genus))
  if (taxa_are_rows(ps_genus)) {
    otu_table <- t(otu_table)
  }
  
  cat("OTU table dimensions:", dim(otu_table), "\n")
  cat("Number of samples:", nrow(otu_table), "\n")
  cat("Number of genera:", ncol(otu_table), "\n")
  
  # Convert to relative abundances
  otu_rel <- otu_table / rowSums(otu_table)
  
  # Compute Jensen-Shannon divergence matrix
  js_matrix <- compute_js_matrix(otu_rel)
  
  # Convert to distance object
  js_dist <- as.dist(js_matrix)
  
  # Test different numbers of clusters
  ch_scores <- numeric(length(k_range))
  pam_results <- list()
  
  cat("\nTesting different numbers of clusters...\n")
  
  for (i in seq_along(k_range)) {
    k <- k_range[i]
    cat("Testing k =", k, "\n")
    
    # Perform PAM clustering
    pam_result <- pam(js_dist, k = k, diss = TRUE)
    pam_results[[i]] <- pam_result
    
    # Calculate Calinski-Harabasz index
    ch_scores[i] <- calinski_harabasz(js_matrix, pam_result$clustering)
    
    cat("  CH index:", ch_scores[i], "\n")
  }
  
  # Find optimal k
  optimal_k_idx <- which.max(ch_scores)
  optimal_k <- k_range[optimal_k_idx]
  optimal_pam <- pam_results[[optimal_k_idx]]
  
  cat("\nOptimal number of clusters:", optimal_k, "\n")
  cat("Best CH index:", max(ch_scores), "\n")
  
  # Create results summary
  results <- list(
    js_matrix = js_matrix,
    js_distance = js_dist,
    k_range = k_range,
    ch_scores = ch_scores,
    optimal_k = optimal_k,
    optimal_clustering = optimal_pam,
    cluster_assignments = optimal_pam$clustering
  )
  
  # Plot CH scores
  ch_plot <- data.frame(k = k_range, CH = ch_scores)
  
  p1 <- ggplot(ch_plot, aes(x = k, y = CH)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
    labs(title = "Calinski-Harabasz Index vs Number of Clusters",
         x = "Number of clusters (k)",
         y = "Calinski-Harabasz Index") +
    theme_minimal()
  
  print(p1)
  
  # Print cluster summary
  cat("\nCluster sizes:\n")
  print(table(optimal_pam$clustering))
  
  return(results)
}

# Run the analysis
cat("Starting community class identification...\n")
community_classes <- identify_community_classes(ps_genus, k_range = 2:8)

# Add cluster assignments to phyloseq metadata
sample_data(ps_genus)$community_class <- community_classes$cluster_assignments

# Test relationship with take-all status
if ("gt" %in% colnames(sample_data(ps_genus))) {
  cat("\nTesting association with take-all status...\n")
  
  cluster_gt_table <- table(
    sample_data(ps_genus)$community_class,
    sample_data(ps_genus)$gt
  )
  
  print(cluster_gt_table)
  
  # Chi-square test
  if (all(cluster_gt_table >= 5)) {
    chi_test <- chisq.test(cluster_gt_table)
    cat("Chi-square test p-value:", chi_test$p.value, "\n")
  } else {
    fisher_test <- fisher.test(cluster_gt_table, simulate.p.value = TRUE)
    cat("Fisher's exact test p-value:", fisher_test$p.value, "\n")
  }
}

# Save results
saveRDS(community_classes, "community_classes_results.rds")

cat("\nAnalysis complete! Results saved to community_classes_results.rds\n")