# ============================================================================
#         ATAC-ONLY CLUSTERING BENCHMARKING SCRIPT (FINALIZED V5)
# ============================================================================
# This script performs a comprehensive benchmarking analysis of different
# clustering methods on ATAC-seq data. It compares standalone algorithms
# with a Poisson Mixture Model initialized by these methods and includes
# hyperparameter tuning for a random initialization based on the Gamma and
# Dirichlet distributions.
#
# V5 Updates:
# - Added UMAP visualization for ground truth labels.
# - Included standalone methods (k-means, Seurat, DBSCAN) in the
#   Log-Likelihood vs. ARI scatter plot.
# - Added ARI value labels to the main performance bar plot.
# - Added a new plot showing the log-likelihood convergence for EM methods.
# - Added a new confusion matrix heatmap for the best-performing method.
# ============================================================================

# ------------------------------ SETUP ------------------------------------
# Set a seed for reproducibility of random processes.
set.seed(42)

# ------------------------------ LIBRARIES ------------------------------------
cat("INFO: Loading required packages...\n")
# Ensure BiocManager is installed for managing Bioconductor packages.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Define the list of required packages for this analysis.
required_packages <- c(
  "Seurat", "SeuratData", "Signac", "EnsDb.Hsapiens.v86", "ggplot2",
  "cowplot", "dplyr", "progress", "mclust", "aricode", "dbscan", "RColorBrewer",
  "ggrepel", "gtools", "tidyr"
)

# Loop through the packages, install if missing, and load them.
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("INFO: Installing package:", pkg, "\n")
    # Use BiocManager for Bioconductor packages.
    if (pkg %in% c("Signac", "EnsDb.Hsapiens.v86")) {
      BiocManager::install(pkg, update = FALSE)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ----------------------------- DATA LOADING & PREPROCESSING -------------------------------
# This section handles loading the PBMC multiome dataset, focusing on the ATAC-seq
# component. It includes standard preprocessing steps like TF-IDF transformation,
# singular value decomposition (SVD), and gene activity scoring.

# Create a progress bar for the data loading and setup phase.
pb_setup <- progress_bar$new(
  format = "  [:bar] :current/:total :message ETA: :eta",
  total = 9, width = 60
)

# Load pbmcMultiome data if not already in the environment
pb_setup$tick(tokens = list(message = "Loading ATAC data..."))
if (!exists("pbmc.atac")) {
  # The SeuratData object is large, so this might take time.
  # If not installed, run: SeuratData::InstallData("pbmcMultiome")
  pbmc.atac <- SeuratData::LoadData("pbmcMultiome", "pbmc.atac")
}

# Filter out low-quality cells.
pbmc.atac <- subset(pbmc.atac, seurat_annotations != "filtered")
pb_setup$tick(tokens = list(message = "Getting annotations..."))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
Annotation(pbmc.atac) <- annotations

pb_setup$tick(tokens = list(message = "Running TF-IDF..."))
pbmc.atac <- RunTFIDF(pbmc.atac)

pb_setup$tick(tokens = list(message = "Finding top features..."))
pbmc.atac <- FindTopFeatures(pbmc.atac, min.cutoff = "q0")

pb_setup$tick(tokens = list(message = "Running SVD..."))
pbmc.atac <- RunSVD(pbmc.atac)

pb_setup$tick(tokens = list(message = "Computing gene activity..."))
gene.activities <- GeneActivity(pbmc.atac)
pbmc.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

pb_setup$tick(tokens = list(message = "Normalizing Gene Activity + PCA..."))
pbmc.atac <- NormalizeData(pbmc.atac, assay = "ACTIVITY", verbose = FALSE)
pbmc.atac <- FindVariableFeatures(pbmc.atac, assay = "ACTIVITY", verbose = FALSE)
pbmc.atac <- ScaleData(pbmc.atac, assay = "ACTIVITY", do.scale = FALSE, verbose = FALSE)
pbmc.atac <- RunPCA(pbmc.atac, assay = "ACTIVITY", verbose = FALSE)
pbmc.atac <- RunUMAP(pbmc.atac, reduction = "pca", dims = 1:30, verbose = FALSE)

pb_setup$tick(tokens = list(message = "Finding neighbors + Seurat clusters..."))
pbmc.atac <- FindNeighbors(pbmc.atac, reduction = "pca", dims = 1:30, graph.name = "pca_nn", verbose = FALSE)
pbmc.atac <- FindClusters(pbmc.atac, graph.name = "pca_nn", resolution = 0.5, verbose = FALSE)

pb_setup$tick(tokens = list(message = "Finalizing data preparation..."))
cat("\nINFO: Data loading and preprocessing complete.\n")


# --------------------------- CLUSTERING INPUTS ------------------------------
cat("INFO: Preparing input matrix for custom clustering...\n")
# Ground truth labels for evaluation.
true_labels <- pbmc.atac$seurat_annotations
# Seurat's default clustering results.
seurat_clusters <- Idents(pbmc.atac)
# Number of clusters (K) based on true cell types.
K <- length(unique(true_labels))
# Select the top 2000 variable genes for clustering.
top2000_genes <- head(VariableFeatures(pbmc.atac[["ACTIVITY"]]), 2000)
# Create the input matrix 'X' from gene activity scores.
X <- as.matrix(GetAssayData(pbmc.atac, assay = "ACTIVITY", slot = "data")[top2000_genes, ])
# Preprocess the matrix: set negative values to 0, round, and remove zero-sum genes.
X[X < 0] <- 0
X <- round(X)
X <- X[rowSums(X) > 0, ]


# ------------------------ EM + INITIALIZATION FUNCTIONS -------------------------------

#' Poisson Mixture Model using Expectation-Maximization
#'
#' @param X Input data matrix (genes x cells).
#' @param K Number of clusters.
#' @param init_lambda Initial lambda parameters (Poisson rates).
#' @param init_pi Initial pi parameters (mixing proportions).
#' @param max_iter Maximum number of iterations.
#' @param tol Convergence tolerance.
#' @param method_name A name for the method for progress bar display.
#' @return A list containing cluster labels, log-likelihood history, and final parameters.
poisson_em <- function(X, K, init_lambda, init_pi, max_iter = 100, tol = 1e-5, method_name = "") {
  N <- ncol(X); D <- nrow(X)
  pi <- init_pi; lambda <- init_lambda
  loglik <- numeric(max_iter)
  
  pb <- progress_bar$new(
    format = paste0("  EM (", method_name, ") [:bar] :percent ETA: :eta"),
    total = max_iter, width = 60, clear = FALSE
  )
  
  for (iter in 1:max_iter) {
    pb$tick()
    # E-step: Calculate responsibilities
    log_r <- sapply(1:K, function(k) log(pi[k] + 1e-10) + colSums(dpois(X, lambda[, k], log = TRUE)))
    log_r_stable <- sweep(log_r, 1, apply(log_r, 1, max), "-")
    r <- exp(log_r_stable)
    r <- r / rowSums(r)
    
    # M-step: Update parameters
    Nk <- colSums(r)
    pi <- Nk / N
    lambda <- (X %*% r) / (matrix(rep(Nk, each = D), D) + 1e-10)
    lambda[lambda < 1e-4] <- 1e-4 # Add floor to prevent numerical issues
    
    # Calculate log-likelihood
    loglik[iter] <- sum(sapply(1:N, function(n) max(log_r[n, ]) + log(sum(exp(log_r[n, ] - max(log_r[n, ]))))))
    
    # Check for convergence
    if (iter > 1 && abs(loglik[iter] - loglik[iter - 1]) < tol * abs(loglik[iter - 1])) {
      loglik <- loglik[1:iter]
      pb$update(1) # Complete the progress bar
      break
    }
  }
  labels <- apply(r, 1, which.max)
  list(labels = labels, loglik = loglik, pi = pi, lambda = lambda)
}


#' Generate Initial Parameters for the EM Algorithm from existing clusters
#'
#' @param X Input data matrix (genes x cells).
#' @param K Number of clusters.
#' @param method The initialization method ('kmeans', 'seurat', 'dbscan', 'true_labels').
#' @param init_clusters A vector of initial cluster assignments.
#' @return A list containing initial pi, lambda, and the (potentially adjusted) K.
get_em_init <- function(X, K, method, init_clusters = NULL) {
  cat(paste0("INFO: Initializing EM with '", method, "'...\n"))
  N <- ncol(X); D <- nrow(X)
  
  if (method == "kmeans") {
    pb <- progress_bar$new(format = "  [:bar] Running k-means for init ETA: :eta", total = 1, width = 60)
    init_clust <- kmeans(t(X), centers = K, nstart = 20)$cluster
    pb$tick()
  } else if (method %in% c("seurat", "dbscan", "true_labels")) {
    if (is.null(init_clusters)) stop("Initial clusters must be provided for this method.")
    init_clust <- as.numeric(factor(init_clusters))
    if (length(unique(init_clust)) != K) {
      warning(paste("Number of initial clusters for", method, "does not match K. Adjusting K."), immediate. = TRUE)
      K <- length(unique(init_clust))
    }
  } else {
    stop("Unsupported initialization method.")
  }
  
  pi <- as.vector(table(factor(init_clust, levels = 1:K))) / N
  pi[is.na(pi)] <- 0
  lambda <- sapply(1:K, function(k) {
    cells_in_clust <- init_clust == k
    if (sum(cells_in_clust) > 1) rowMeans(X[, cells_in_clust, drop = FALSE])
    else if (sum(cells_in_clust) == 1) X[, cells_in_clust]
    else runif(D, 0.1, 1)
  })
  lambda[lambda < 1e-4] <- 1e-4
  list(pi = pi, lambda = lambda, K = K)
}

#' Calculate Log-Likelihood for a given clustering under a Poisson Mixture Model
#'
#' @param X Input data matrix (genes x cells).
#' @param labels A vector of cluster assignments for each cell.
#' @return The total log-likelihood of the data given the clustering.
calculate_poisson_loglik <- function(X, labels) {
  N <- ncol(X); D <- nrow(X)
  
  # Ensure labels are numeric factors starting from 1
  labels_factor <- factor(labels)
  K <- length(levels(labels_factor))
  labels_numeric <- as.numeric(labels_factor)
  
  # Calculate pi and lambda from the given labels
  pi <- as.vector(table(labels_numeric)) / N
  lambda <- sapply(1:K, function(k) {
    cells_in_clust <- labels_numeric == k
    if (sum(cells_in_clust) > 1) rowMeans(X[, cells_in_clust, drop = FALSE])
    else if (sum(cells_in_clust) == 1) X[, cells_in_clust]
    else rep(1e-4, D) # Should not happen if K is derived from labels
  })
  lambda[lambda < 1e-4] <- 1e-4
  
  # Calculate the log-likelihood
  log_probs <- sapply(1:K, function(k) log(pi[k] + 1e-10) + colSums(dpois(X, lambda[, k], log = TRUE)))
  
  # Log-sum-exp for numerical stability
  loglik <- sum(sapply(1:N, function(n) {
    max_log_prob <- max(log_probs[n, ])
    max_log_prob + log(sum(exp(log_probs[n, ] - max_log_prob)))
  }))
  
  return(loglik)
}


# -------------------------- COMPREHENSIVE ANALYSIS FUNCTION ----------------------------

#' Run a comprehensive clustering analysis and generate visualizations.
#'
#' @param X The input data matrix (genes x cells).
#' @param true_labels The ground truth cell type labels.
#' @param seurat_labels The labels from the standard Seurat clustering workflow.
#' @param K The target number of clusters.
#' @param analysis_obj The Seurat object containing UMAP and PCA reductions.
#' @param dataset_name A name for the dataset to use in filenames.
run_comprehensive_analysis <- function(X, true_labels, seurat_labels, K, analysis_obj, dataset_name) {
  cat(paste0("\n\n========== COMPREHENSIVE ANALYSIS: ", dataset_name, " ==========\n"))
  true_labels_numeric <- as.numeric(factor(true_labels))
  results <- list(); aris <- list(); logliks <- list()
  D <- nrow(X) # Number of genes
  
  # --- 1. Standalone Methods ---
  cat("\nINFO: Running standalone clustering methods...\n")
  pb_standalone <- progress_bar$new(format = "  [:bar] Standalone Methods ETA: :eta", total = 2, width = 60)
  results$kmeans <- kmeans(t(X), K, nstart = 20)$cluster
  aris$kmeans <- aricode::ARI(results$kmeans, true_labels_numeric)
  logliks$kmeans <- calculate_poisson_loglik(X, results$kmeans)
  pb_standalone$tick()
  
  results$seurat <- as.numeric(factor(seurat_labels))
  aris$seurat <- aricode::ARI(results$seurat, true_labels_numeric)
  logliks$seurat <- calculate_poisson_loglik(X, results$seurat)
  
  pca_coords <- analysis_obj@reductions$pca@cell.embeddings[, 1:30]
  eps_grid <- c(1.5, 2, 2.5, 3, 3.5); minPts_grid <- c(10, 15, 20, 25, 30)
  dbscan_params <- expand.grid(eps = eps_grid, minPts = minPts_grid)
  dbscan_aris <- sapply(1:nrow(dbscan_params), function(i) {
    db_res <- dbscan(pca_coords, eps = dbscan_params$eps[i], minPts = dbscan_params$minPts[i])
    aricode::ARI(db_res$cluster, true_labels_numeric)
  })
  best_param_idx <- which.max(dbscan_aris)
  best_params <- dbscan_params[best_param_idx, ]
  cat(paste0("INFO: Best DBSCAN params: eps=", best_params$eps, ", minPts=", best_params$minPts, " (ARI=", round(max(dbscan_aris), 3), ")\n"))
  results$dbscan <- dbscan(pca_coords, eps = best_params$eps, minPts = best_params$minPts)$cluster
  aris$dbscan <- max(dbscan_aris)
  logliks$dbscan <- calculate_poisson_loglik(X, results$dbscan)
  pb_standalone$tick()
  
  # Add ground truth for plotting
  results$truth <- true_labels
  aris$truth <- 1.0
  logliks$truth <- calculate_poisson_loglik(X, true_labels_numeric)
  
  # --- 2. EM with Deterministic Initializations ---
  cat("INFO: Running Poisson EM with deterministic initializations...\n")
  
  init_kmeans <- get_em_init(X, K, method = "kmeans")
  results$em_kmeans <- poisson_em(X, init_kmeans$K, init_kmeans$lambda, init_kmeans$pi, method_name = "k-means init")
  aris$em_kmeans <- aricode::ARI(results$em_kmeans$labels, true_labels_numeric)
  logliks$em_kmeans <- tail(results$em_kmeans$loglik, 1)
  
  init_seurat <- get_em_init(X, K, method = "seurat", init_clusters = results$seurat)
  results$em_seurat <- poisson_em(X, init_seurat$K, init_seurat$lambda, init_seurat$pi, method_name = "Seurat init")
  aris$em_seurat <- aricode::ARI(results$em_seurat$labels, true_labels_numeric)
  logliks$em_seurat <- tail(results$em_seurat$loglik, 1)
  
  dbscan_clusters_for_init <- results$dbscan
  if(any(dbscan_clusters_for_init == 0)) dbscan_clusters_for_init[dbscan_clusters_for_init == 0] <- max(dbscan_clusters_for_init) + 1
  K_dbscan <- length(unique(dbscan_clusters_for_init))
  init_dbscan <- get_em_init(X, K_dbscan, method = "dbscan", init_clusters = dbscan_clusters_for_init)
  results$em_dbscan <- poisson_em(X, init_dbscan$K, init_dbscan$lambda, init_dbscan$pi, method_name = "DBSCAN init")
  aris$em_dbscan <- aricode::ARI(results$em_dbscan$labels, true_labels_numeric)
  logliks$em_dbscan <- tail(results$em_dbscan$loglik, 1)
  
  init_truth <- get_em_init(X, K, method = "true_labels", init_clusters = true_labels)
  results$em_truth <- poisson_em(X, init_truth$K, init_truth$lambda, init_truth$pi, method_name = "Truth init")
  aris$em_truth <- aricode::ARI(results$em_truth$labels, true_labels_numeric)
  logliks$em_truth <- tail(results$em_truth$loglik, 1)
  
  # --- 3. EM with Random Initialization (Gamma/Dirichlet hyperparameter tuning) ---
  cat("INFO: Tuning EM with random Gamma/Dirichlet-based initializations...\n")
  alpha_grid <- c(0.1, 0.5, 1, 2, 5)
  beta_grid <- c(0.1, 0.5, 1, 2, 5)
  random_params <- expand.grid(alpha = alpha_grid, beta = beta_grid)
  random_tuning_results <- data.frame()
  
  best_random_ari <- -1
  best_random_result <- NULL
  best_random_params <- NULL
  
  pb_random <- progress_bar$new(format = "  Random Init Tuning [:bar] :percent ETA: :eta", total = nrow(random_params), width = 60)
  
  for(i in 1:nrow(random_params)) {
    params <- random_params[i, ]
    # Initialize pi from a Dirichlet(1,1,...,1) distribution and lambda from Gamma
    rand_pi <- gtools::rdirichlet(1, rep(1, K))[1, ]
    rand_lambda <- sapply(1:K, function(k) rgamma(D, shape = params$alpha, rate = params$beta))
    rand_lambda[rand_lambda < 1e-4] <- 1e-4
    
    current_result <- poisson_em(X, K, rand_lambda, rand_pi, method_name = paste0("a=", params$alpha, ", b=", params$beta))
    current_ari <- aricode::ARI(current_result$labels, true_labels_numeric)
    
    random_tuning_results <- rbind(random_tuning_results, data.frame(alpha = params$alpha, beta = params$beta, ARI = current_ari))
    
    if (current_ari > best_random_ari) {
      best_random_ari <- current_ari
      best_random_result <- current_result
      best_random_params <- params
    }
    pb_random$tick()
  }
  
  cat(paste0("INFO: Best random init found with alpha=", best_random_params$alpha, ", beta=", best_random_params$beta, " (ARI=", round(best_random_ari, 3), ")\n"))
  results$em_random_best <- best_random_result
  aris$em_random_best <- best_random_ari
  logliks$em_random_best <- tail(best_random_result$loglik, 1)
  
  # --- 4. Generate and Save Plots ---
  cat("INFO: Generating and saving visualizations...\n")
  output_dir <- file.path(Sys.getenv("HOME"), "Desktop", "ATAC_Clustering_Results")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # PLOT 1: ARI Barplot
  ari_df <- data.frame(Method = names(aris), ARI = unlist(aris))
  ari_df$Type <- case_when(
    grepl("^em_random", ari_df$Method) ~ "EM (Best Random Init)",
    grepl("^em_", ari_df$Method) ~ "EM (Informed Init)",
    ari_df$Method == "truth" ~ "Ground Truth",
    TRUE ~ "Standalone"
  )
  p_ari <- ggplot(ari_df, aes(x = reorder(Method, ARI), y = ARI, fill = Type)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.8) +
    geom_text(aes(label = round(ARI, 3)), hjust = -0.2, size = 3.5, color = "black") +
    coord_flip(ylim = c(0, 1.1)) + theme_minimal(base_size = 14) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = paste("Clustering Performance (ARI) -", dataset_name), x = "Clustering Method", y = "Adjusted Rand Index (ARI)") +
    theme(legend.position = "bottom", plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(output_dir, paste0(dataset_name, "_1_ari_barplot.png")), plot = p_ari, width = 12, height = 8, dpi = 300)
  
  # PLOT 2: Final Log-Likelihood vs. ARI
  ll_vs_ari_df <- data.frame(Method = names(logliks), Final_LogLikelihood = unlist(logliks), ARI = unlist(aris[names(logliks)]))
  p_ll_ari <- ggplot(ll_vs_ari_df, aes(x = Final_LogLikelihood, y = ARI, color = Method)) +
    geom_point(size = 5, alpha = 0.8) +
    ggrepel::geom_text_repel(aes(label = Method), size = 3.5, box.padding = 0.5, max.overlaps = Inf) +
    theme_minimal(base_size = 12) +
    labs(title = "Final Log-Likelihood vs. Clustering Performance", x = "Final Log-Likelihood (under Poisson Model)", y = "Adjusted Rand Index (ARI)") +
    theme(legend.position = "none", plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(output_dir, paste0(dataset_name, "_2_loglik_vs_ari.png")), plot = p_ll_ari, width = 10, height = 8, dpi = 300)
  
  # PLOT 3: Log-Likelihood Convergence
  em_logliks_full <- list(
    em_kmeans = results$em_kmeans$loglik,
    em_seurat = results$em_seurat$loglik,
    em_dbscan = results$em_dbscan$loglik,
    em_truth = results$em_truth$loglik,
    em_random_best = results$em_random_best$loglik
  )
  loglik_df <- purrr::map_dfr(em_logliks_full, ~data.frame(loglik = .x, iter = 1:length(.x)), .id = "Method")
  p_convergence <- ggplot(loglik_df, aes(x = iter, y = loglik, color = Method, group = Method)) +
    geom_line(linewidth = 1.2) + geom_point(size = 2) +
    theme_minimal(base_size = 12) + scale_color_brewer(palette = "Set1") +
    labs(title = "EM Log-Likelihood Convergence", x = "Iteration", y = "Log-Likelihood") +
    theme(legend.position = "bottom", plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(output_dir, paste0(dataset_name, "_3_convergence_plot.png")), plot = p_convergence, width = 10, height = 8, dpi = 300)
  
  # PLOT 4: UMAP Visualizations with Legends
  umap_df <- as.data.frame(analysis_obj@reductions$umap@cell.embeddings)
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  
  map_cluster_to_truth <- function(pred_labels, true_labels) {
    df <- data.frame(pred = pred_labels, true = true_labels)
    mapping <- df %>% group_by(pred) %>% summarise(top_true = names(which.max(table(true)))) %>% mutate(new_label = paste0(pred, ": ", top_true))
    new_labels <- as.character(pred_labels)
    for(i in 1:nrow(mapping)) { new_labels[pred_labels == mapping$pred[i]] <- mapping$new_label[i] }
    return(factor(new_labels))
  }
  
  all_methods <- names(aris)
  for(method in all_methods) {
    raw_labels <- if(grepl("^em_", method)) results[[method]]$labels else results[[method]]
    if(method == "dbscan" && any(raw_labels == 0)) { raw_labels[raw_labels == 0] <- "Noise" }
    
    if(method == "truth") {
      umap_df[[method]] <- factor(raw_labels)
    } else {
      umap_df[[method]] <- map_cluster_to_truth(raw_labels, true_labels)
    }
  }
  
  plot_list <- list()
  for (method in all_methods) {
    n_colors <- length(unique(umap_df[[method]]))
    color_palette <- colorRampPalette(RColorBrewer::brewer.pal(min(n_colors, 9), "Set1"))(n_colors)
    p <- ggplot(umap_df, aes_string(x = "UMAP_1", y = "UMAP_2", color = method)) +
      geom_point(size = 0.3, alpha = 0.8) + scale_color_manual(values = color_palette, name = "Cluster Label") +
      theme_void() + labs(title = method) +
      theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"), legend.position = "right",
            legend.key.size = unit(0.4, "cm"), legend.text = element_text(size=8),
            legend.title = element_text(size=10, face="bold"), plot.background = element_rect(fill = "white", color = NA))
    plot_list[[method]] <- p
  }
  
  p_umap_grid <- cowplot::plot_grid(plotlist = plot_list, ncol = 3)
  ggsave(file.path(output_dir, paste0(dataset_name, "_4_umap_grid.png")), plot = p_umap_grid, width = 24, height = (ceiling(length(plot_list)/3) * 6), dpi = 300, limitsize = FALSE)
  cat("INFO: UMAP grid saved. Note: This image may be very large to accommodate all legends.\n")
  
  # PLOT 5: Confusion Matrix for Best Method
  best_method_name <- ari_df %>% filter(Method != "truth" & Method != "em_truth") %>% arrange(desc(ARI)) %>% head(1) %>% pull(Method)
  best_labels <- if(grepl("^em_", best_method_name)) results[[best_method_name]]$labels else results[[best_method_name]]
  
  conf_mat_df <- as.data.frame(table(Predicted = best_labels, Truth = true_labels))
  
  p_conf_mat <- ggplot(conf_mat_df, aes(x = Truth, y = Predicted, fill = Freq)) +
    geom_tile(color = "black") +
    geom_text(aes(label = Freq), color = "white", size = 3) +
    scale_fill_gradient(low = "#0d3b66", high = "#faf0ca") +
    theme_minimal(base_size = 12) +
    labs(title = paste("Confusion Matrix: Best Method (", best_method_name, ") vs. Ground Truth"),
         x = "Ground Truth Labels", y = "Predicted Cluster Labels") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.background = element_rect(fill = "white", color = NA))
  ggsave(file.path(output_dir, paste0(dataset_name, "_5_confusion_matrix.png")), plot = p_conf_mat, width = 12, height = 10, dpi = 300)
  
  # --- 5. Save RDS results ---
  cat("INFO: Saving final results object...\n")
  saveRDS(
    list(results = results, ARIs = aris, logliks = logliks, dbscan_best_params = best_params, 
         random_em_best_params = best_random_params, random_tuning_results = random_tuning_results),
    file.path(output_dir, paste0(dataset_name, "_comprehensive_results.rds"))
  )
  cat(paste0("\n✅ Completed comprehensive analysis for ", dataset_name, "\n"))
}

# ----------------------------- RUN THE ANALYSIS ----------------------------
run_comprehensive_analysis(
  X = X,
  true_labels = true_labels,
  seurat_labels = seurat_clusters,
  K = K,
  analysis_obj = pbmc.atac,
  dataset_name = "ATAC_GeneActivity_Top2000"
)
