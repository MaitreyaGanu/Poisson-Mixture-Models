#_______________________________________________________________________________
# POISSON MIXTURE MODELS AND OTHER CLUSTERING METHODS ON PBMC3K - FIXED K
#_______________________________________________________________________________

# Ensure all necessary packages are installed and loaded
# install.packages(c("Seurat", "SeuratData", "mclust", "ggplot2", "dplyr", "RColorBrewer", "dbscan", "patchwork", "reshape2"))

# Load libraries
library(Seurat)
library(SeuratData) # For pbmc3k dataset
library(mclust)     # For adjustedRandIndex and Mclust
library(ggplot2)    # For plotting
library(dplyr)      # For data manipulation
library(RColorBrewer)# For plot colors
library(dbscan)     # For DBSCAN clustering
library(patchwork)  # For combining plots
library(reshape2)   # For melt function

# Set seed for reproducibility
set.seed(12345)

#_______________________________________________________________________________
# 1. LOAD AND PREPROCESS DATA
#_______________________________________________________________________________
cat("Loading and preprocessing PBMC3k data...\n")

# Load PBMC3k data using SeuratData
# If not already installed, run: InstallData("pbmc3k")
pbmc <- LoadData("pbmc3k", type = "pbmc3k.final") # Using pbmc3k.final for annotations

# Get raw counts
counts <- GetAssayData(pbmc, assay = "RNA", slot = "counts")
X_raw <- t(as.matrix(counts)) # Cells as rows, genes as columns

cat(paste("Raw data dimensions: ", nrow(X_raw), " cells, ", ncol(X_raw), " genes.\n"))

#_______________________________________________________________________________
# PRELIMINARY EXPLORATION OF RAW DATA
#_______________________________________________________________________________
cat("\nPerforming preliminary exploration of raw data (X_raw)...\n")

# 1. Histogram of non-zero counts in X_raw
non_zero_counts_raw <- X_raw[X_raw > 0]
plot_hist_raw_counts <- ggplot(data.frame(counts = non_zero_counts_raw), aes(x = counts)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000), labels = scales::comma) +
  labs(title = "Histogram of Non-Zero Counts in Raw Data (X_raw)",
       x = "Counts (log10 scale)", y = "Frequency") +
  theme_minimal()
print(plot_hist_raw_counts)
cat("Summary of non-zero counts in X_raw:\n")
print(summary(non_zero_counts_raw))

# 2. Bar plot of total counts for top N most expressed genes in X_raw
total_counts_per_gene_raw <- colSums(X_raw)
top_n_genes <- 30
top_genes_df_raw <- data.frame(
  gene = names(sort(total_counts_per_gene_raw, decreasing = TRUE)[1:top_n_genes]),
  total_counts = as.numeric(sort(total_counts_per_gene_raw, decreasing = TRUE)[1:top_n_genes])
)
top_genes_df_raw$gene <- factor(top_genes_df_raw$gene, levels = top_genes_df_raw$gene)

plot_top_genes_raw <- ggplot(top_genes_df_raw, aes(x = gene, y = total_counts)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = paste("Top", top_n_genes, "Most Expressed Genes in Raw Data (X_raw)"),
       x = "Gene", y = "Total Counts")
print(plot_top_genes_raw)

# 3. Histogram of library sizes (total counts per cell) from X_raw
library_sizes_raw <- rowSums(X_raw)
plot_hist_library_sizes_raw <- ggplot(data.frame(library_size = library_sizes_raw), aes(x = library_size)) +
  geom_histogram(binwidth = 100, fill = "lightgreen", color = "black") +
  scale_x_continuous(labels = scales::comma) +
  labs(title = "Histogram of Library Sizes (Total Counts per Cell) in Raw Data (X_raw)",
       x = "Library Size", y = "Number of Cells") +
  theme_minimal()
print(plot_hist_library_sizes_raw)
cat("Summary of library sizes in X_raw:\n")
print(summary(library_sizes_raw))


# Use top 1000 variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
var_genes <- VariableFeatures(pbmc)
X <- X_raw[, var_genes]
N <- nrow(X) # Number of cells
D <- ncol(X) # Number of genes

cat(paste("Data for PMM (X): ", N, " cells, ", D, " variable genes.\n"))

# True cluster labels from Seurat annotations
if ("seurat_annotations" %in% colnames(pbmc@meta.data)) {
  true_labels <- factor(pbmc$seurat_annotations)
} else {
  warning("Seurat annotations not found in pbmc$seurat_annotations. Using pbmc$seurat_clusters if available.")
  if("seurat_clusters" %in% colnames(pbmc@meta.data)) {
    true_labels <- factor(pbmc$seurat_clusters)
  } else {
    stop("Cannot find true cluster labels (seurat_annotations or seurat_clusters) in Seurat object.")
  }
}
K_true <- length(levels(true_labels)) # Should be 9 for pbmc3k.final
cat(paste("Using fixed number of clusters K =", K_true, "(from Seurat annotations)\n"))
print(table(true_labels))


#_______________________________________________________________________________
# 2. HELPER FUNCTION: LOG POISSON PMF
#_______________________________________________________________________________
log_poisson_pmf_sum <- function(x_row, lambda_k) {
  sum(x_row * log(lambda_k + 1e-10) - lambda_k - lgamma(x_row + 1))
}

#_______________________________________________________________________________
# 3. POISSON MIXTURE MODEL - EXPECTATION-MAXIMIZATION (EM)
#_______________________________________________________________________________
cat("\nStarting PMM with EM algorithm for K =", K_true, "...\n")

fit_pmm_em <- function(X, K, max_iter = 100, tol = 1e-5) {
  N <- nrow(X)
  D <- ncol(X)
  pi_k <- rep(1 / K, K)
  initial_assignments <- sample(1:K, N, replace = TRUE)
  lambda <- matrix(1, nrow = K, ncol = D)
  for(k_init in 1:K){
    if(sum(initial_assignments == k_init) > 0){
      lambda[k_init,] <- colMeans(X[initial_assignments == k_init, , drop = FALSE]) + 1e-5
    } else {
      lambda[k_init,] <- colMeans(X) + runif(D, 0, 0.1) + 1e-5
    }
    lambda[k_init, lambda[k_init,] <= 0] <- 1e-5
  }
  logliks <- numeric(max_iter)
  current_loglik <- -Inf 
  for (iter in 1:max_iter) {
    log_r_nk_unnormalized <- matrix(0, nrow = N, ncol = K)
    for (k_idx in 1:K) {
      log_prob_X_given_lambda_k <- apply(X, 1, log_poisson_pmf_sum, lambda_k = lambda[k_idx, ])
      log_r_nk_unnormalized[, k_idx] <- log(pi_k[k_idx] + 1e-10) + log_prob_X_given_lambda_k
    }
    max_log_r_nk <- apply(log_r_nk_unnormalized, 1, max)
    log_sum_exp_r_nk <- max_log_r_nk + log(rowSums(exp(log_r_nk_unnormalized - max_log_r_nk)))
    r_nk <- exp(log_r_nk_unnormalized - log_sum_exp_r_nk)
    r_nk <- r_nk / rowSums(r_nk)
    N_k <- colSums(r_nk)
    pi_k <- N_k / N
    pi_k[pi_k < 1e-10] <- 1e-10
    pi_k <- pi_k / sum(pi_k)
    for (k_idx in 1:K) {
      if (N_k[k_idx] > 1e-5) {
        lambda[k_idx, ] <- colSums(r_nk[, k_idx] * X) / N_k[k_idx]
      } else {
        lambda[k_idx, ] <- colMeans(X) + runif(D, 0, 0.1) + 1e-5
      }
      lambda[k_idx, lambda[k_idx,] <= 0] <- 1e-5
    }
    new_loglik <- sum(log_sum_exp_r_nk) 
    logliks[iter] <- new_loglik
    if (iter > 1 && !is.na(new_loglik) && !is.na(logliks[iter-1])) { 
      if (abs(new_loglik - logliks[iter-1]) < tol * abs(logliks[iter-1]) || abs(new_loglik - logliks[iter-1]) < tol ) {
        current_loglik <- new_loglik 
        break
      }
    }
    current_loglik <- new_loglik 
  }
  lambda[lambda < 1e-10] <- 1e-10
  list(lambda = lambda, pi_k = pi_k, responsibilities = r_nk,
       logLik_trajectory = logliks[1:iter], logLik = current_loglik, K = K)
}

em_model_fixed_k <- fit_pmm_em(X, K_true, max_iter = 150, tol = 1e-5)
clusters_em <- apply(em_model_fixed_k$responsibilities, 1, which.max)
loglik_em_fixed_k <- em_model_fixed_k$logLik

cat(paste("EM (K=", K_true, "): LogLik =", round(loglik_em_fixed_k, 2), "\n"))
cat(paste("ARI for EM (K=", K_true, ") vs Seurat Annotations: ",
          round(adjustedRandIndex(true_labels, as.factor(clusters_em)), 3), "\n"))

#_______________________________________________________________________________
# 4. POISSON MIXTURE MODEL - GIBBS SAMPLING (BAYESIAN)
#_______________________________________________________________________________
cat("\nStarting PMM with Gibbs Sampling (Bayesian) for K =", K_true, "...\n")

fit_pmm_gibbs_bayesian <- function(X, K, n_iter = 500, burn_in = 100,
                                   alpha_0 = 1.0, beta_0 = 1.0,
                                   delta_0 = 1.0) {
  N <- nrow(X)
  D <- ncol(X)
  Z <- sample(1:K, N, replace = TRUE)
  pi_k <- as.numeric(table(factor(Z, levels = 1:K))) / N
  if(any(pi_k == 0)) pi_k[pi_k==0] <- 1e-5; pi_k <- pi_k/sum(pi_k)
  lambda <- matrix(0, nrow = K, ncol = D)
  for (k_idx in 1:K) {
    cells_in_k <- (Z == k_idx)
    if (sum(cells_in_k) > 0) {
      lambda[k_idx, ] <- colMeans(X[cells_in_k, , drop = FALSE]) + runif(D, 0, 0.1)
    } else {
      lambda[k_idx, ] <- rgamma(D, shape = alpha_0, rate = beta_0)
    }
    lambda[k_idx, lambda[k_idx,] <= 1e-10] <- 1e-10
  }
  for (iter in 1:n_iter) {
    log_prob_Z_unnorm <- matrix(0, nrow = N, ncol = K)
    for (k_idx in 1:K) {
      log_prob_X_given_lambda_k <- apply(X, 1, log_poisson_pmf_sum, lambda_k = lambda[k_idx, ])
      log_prob_Z_unnorm[, k_idx] <- log(pi_k[k_idx] + 1e-10) + log_prob_X_given_lambda_k
    }
    prob_Z <- exp(log_prob_Z_unnorm - apply(log_prob_Z_unnorm, 1, max, na.rm=TRUE))
    prob_Z <- prob_Z / rowSums(prob_Z, na.rm=TRUE)
    prob_Z[is.na(prob_Z)] <- 1/K
    for (i in 1:N) {
      if(all(is.na(prob_Z[i,])) || sum(prob_Z[i,], na.rm=TRUE) == 0) {
        Z[i] <- sample(1:K, 1)
      } else {
        Z[i] <- sample(1:K, 1, prob = prob_Z[i, ])
      }
    }
    N_k <- as.numeric(table(factor(Z, levels = 1:K)))
    gamma_samples <- rgamma(K, shape = N_k + delta_0, rate = 1)
    pi_k <- gamma_samples / sum(gamma_samples)
    pi_k[pi_k < 1e-10] <- 1e-10
    pi_k <- pi_k / sum(pi_k)
    for (k_idx in 1:K) {
      cells_in_k_mask <- (Z == k_idx)
      N_k_current <- sum(cells_in_k_mask)
      if (N_k_current > 0) {
        sum_X_in_k <- colSums(X[cells_in_k_mask, , drop = FALSE])
        alpha_posterior <- alpha_0 + sum_X_in_k
        beta_posterior <- beta_0 + N_k_current
        lambda[k_idx, ] <- rgamma(D, shape = alpha_posterior, rate = beta_posterior)
      } else {
        lambda[k_idx, ] <- rgamma(D, shape = alpha_0, rate = beta_0)
      }
      lambda[k_idx, lambda[k_idx,] <= 1e-10] <- 1e-10
    }
  }
  log_r_nk_unnormalized_final <- matrix(0, nrow = N, ncol = K)
  for (k_idx in 1:K) {
    log_prob_X_given_lambda_k_final <- apply(X, 1, log_poisson_pmf_sum, lambda_k = lambda[k_idx, ])
    log_r_nk_unnormalized_final[, k_idx] <- log(pi_k[k_idx] + 1e-10) + log_prob_X_given_lambda_k_final
  }
  max_log_r_nk_final <- apply(log_r_nk_unnormalized_final, 1, max, na.rm=TRUE)
  log_sum_exp_r_nk_final <- max_log_r_nk_final + log(rowSums(exp(log_r_nk_unnormalized_final - max_log_r_nk_final), na.rm=TRUE))
  final_loglik <- sum(log_sum_exp_r_nk_final, na.rm=TRUE)
  lambda[lambda < 1e-10] <- 1e-10
  list(clusters = Z, lambda = lambda, pi_k = pi_k,
       logLik = final_loglik, K = K,
       alpha_0 = alpha_0, beta_0 = beta_0)
}

default_alpha_0 <- 1.0
default_beta_0 <- 1.0
gibbs_model_fixed_k <- fit_pmm_gibbs_bayesian(X, K_true,
                                              n_iter = 600, burn_in = 200,
                                              alpha_0 = default_alpha_0, beta_0 = default_beta_0)
clusters_gibbs <- gibbs_model_fixed_k$clusters
loglik_gibbs_fixed_k <- gibbs_model_fixed_k$logLik

cat(paste("Gibbs (K=", K_true, "): LogLik =", round(loglik_gibbs_fixed_k, 2), "\n"))
cat(paste("ARI for Gibbs (K=", K_true, ") vs Seurat Annotations: ",
          round(adjustedRandIndex(true_labels, as.factor(clusters_gibbs)), 3), "\n"))
#_______________________________________________________________________________
# 6. DIMENSION REDUCTION & OTHER CLUSTERING METHODS
#_______________________________________________________________________________
cat("\nPerforming dimension reduction...\n")

# Standard Seurat workflow for PCA and UMAP (run on var_genes)
pbmc <- ScaleData(pbmc, features = var_genes, verbose = FALSE) # Scale only variable features for PCA/UMAP
pbmc <- RunPCA(pbmc, features = var_genes, npcs = 30, verbose = FALSE) # Using 30 PCs e.g.
pbmc <- RunUMAP(pbmc, dims = 1:10, verbose = FALSE) # Using 10 PCs for UMAP

# Ensure cell name consistency for adding cluster labels
cells_in_X <- rownames(X) # These are the cells used for PMM modeling
cells_in_pbmc_metadata <- rownames(pbmc@meta.data)
if(!all(cells_in_X %in% cells_in_pbmc_metadata)){
  warning("Some cells in matrix X are not found in pbmc metadata for PMM. This should not happen if X comes from pbmc.")
}

# Add PMM cluster assignments to Seurat object metadata
# Match cells from X (used for PMM) to the full pbmc object
pbmc$PMM_EM_Clusters <- NA # Initialize column
idx_match_em <- match(cells_in_X, cells_in_pbmc_metadata)
pbmc$PMM_EM_Clusters[idx_match_em] <- factor(clusters_em)

pbmc$PMM_Gibbs_Clusters <- NA # Initialize column
idx_match_gibbs <- match(cells_in_X, cells_in_pbmc_metadata)
pbmc$PMM_Gibbs_Clusters[idx_match_gibbs] <- factor(clusters_gibbs)

cat("\nApplying other clustering methods (KMeans, GMM, DBSCAN)...\n")

# --- K-Means Clustering (on UMAP embeddings) ---
cat("Running K-Means...\n")
umap_emb <- Embeddings(pbmc, "umap")
set.seed(123)
kmeans_result <- kmeans(umap_emb, centers = K_true, nstart = 25)
pbmc$KMeans_Clusters <- factor(kmeans_result$cluster)
# Ensure true_labels for ARI match the cells in umap_emb (which should be all cells in pbmc)
ari_kmeans <- adjustedRandIndex(true_labels[rownames(umap_emb)], pbmc$KMeans_Clusters[rownames(umap_emb)])
cat(paste("ARI for K-Means (K=", K_true, ") vs Seurat Annotations: ", round(ari_kmeans, 3), "\n"))

# --- Gaussian Mixture Model (GMM) (on PCA embeddings) ---
cat("Running Gaussian Mixture Model (GMM)...\n")
# Use the same number of PCs as used for UMAP, e.g., 10
num_pcs_for_gmm <- 10 
pca_emb_gmm <- Embeddings(pbmc, "pca")[, 1:num_pcs_for_gmm]
gmm_model <- Mclust(pca_emb_gmm, G = K_true, verbose = FALSE) # Specify K_true clusters
pbmc$GMM_Clusters <- factor(gmm_model$classification)
# Ensure true_labels for ARI match the cells in pca_emb_gmm (all cells in pbmc)
ari_gmm <- adjustedRandIndex(true_labels[rownames(pca_emb_gmm)], pbmc$GMM_Clusters[rownames(pca_emb_gmm)])
cat(paste("ARI for GMM (K=", K_true, ") vs Seurat Annotations: ", round(ari_gmm, 3), "\n"))
cat(paste("GMM found model:", gmm_model$modelName, "with", gmm_model$G, "clusters.\n"))

# --- DBSCAN Clustering (on UMAP embeddings) ---
cat("Running DBSCAN...\n")
# DBSCAN does not take K as an input. Parameters eps and minPts determine cluster number.
dbscan_params <- list(eps = 0.5, minPts = 10) # Using previously chosen parameters
dbscan_result <- dbscan::dbscan(umap_emb, eps = dbscan_params$eps, minPts = dbscan_params$minPts)
pbmc$DBSCAN_Clusters <- factor(dbscan_result$cluster) # Cluster 0 is noise
num_dbscan_clusters <- length(unique(dbscan_result$cluster[dbscan_result$cluster != 0]))
cat(paste("DBSCAN (eps=", dbscan_params$eps, ", minPts=", dbscan_params$minPts, ") found ", num_dbscan_clusters, " clusters (excluding noise).\n"))
# For ARI, subset true_labels to only those cells that were not classified as noise by DBSCAN
non_noise_dbscan_indices <- which(pbmc$DBSCAN_Clusters[rownames(umap_emb)] != 0)
if(length(non_noise_dbscan_indices) > 0 && num_dbscan_clusters > 1) {
  ari_dbscan <- adjustedRandIndex(true_labels[rownames(umap_emb)][non_noise_dbscan_indices],
                                  pbmc$DBSCAN_Clusters[rownames(umap_emb)][non_noise_dbscan_indices])
  cat(paste("ARI for DBSCAN vs Seurat Annotations: ", round(ari_dbscan, 3), " (on ", num_dbscan_clusters," clusters, excluding noise points).\n"))
} else {
  cat("DBSCAN classified all points as noise or resulted in too few clusters for ARI calculation.\n")
}

#_______________________________________________________________________________
# 7. UMAP VISUALIZATIONS
#_______________________________________________________________________________
cat("\nGenerating UMAP plots...\n")

plot_umap_seurat <- DimPlot(pbmc, reduction = "umap", group.by = "seurat_annotations", repel = TRUE, label.size = 3) +
  ggtitle("Seurat Annotations (True Clusters)")
plot_umap_em <- DimPlot(pbmc, reduction = "umap", group.by = "PMM_EM_Clusters", repel = TRUE, label.size = 3) +
  ggtitle(paste("PMM EM (K =", K_true, ")"))
plot_umap_gibbs <- DimPlot(pbmc, reduction = "umap", group.by = "PMM_Gibbs_Clusters", repel = TRUE, label.size = 3) +
  ggtitle(paste("PMM Gibbs (K =", K_true, ")"))
plot_umap_kmeans <- DimPlot(pbmc, reduction = "umap", group.by = "KMeans_Clusters", repel = TRUE, label.size = 3) +
  ggtitle(paste("KMeans (K =", K_true, ")"))
plot_umap_gmm <- DimPlot(pbmc, reduction = "umap", group.by = "GMM_Clusters", repel = TRUE, label.size = 3) +
  ggtitle(paste("GMM (K =", K_true, ")"))
plot_umap_dbscan <- DimPlot(pbmc, reduction = "umap", group.by = "DBSCAN_Clusters", repel = TRUE, label.size = 3) +
  ggtitle(paste("DBSCAN (", num_dbscan_clusters, " clusters found)", sep=""))

# Arrange plots in a 2x3 grid
combined_umaps <- (plot_umap_seurat | plot_umap_em | plot_umap_gibbs) /
  (plot_umap_kmeans | plot_umap_gmm | plot_umap_dbscan)

print(combined_umaps + plot_layout(guides = "collect") & theme(legend.position = 'bottom'))

