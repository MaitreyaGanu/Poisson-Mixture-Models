library(flexmix)
library(mclust)
library(ggplot2)

set.seed(42)

# Parameters
N <- 1000
K_true <- 5
true_features <- 5
alpha <- 5
beta <- 1
EM_INIT <- 15
K_range <- 1:10

# Log-Sum-Exp Trick
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Synthetic Dataset Generator
generate_dataset <- function(noise_features) {
  D <- true_features + noise_features
  z <- sample(1:K_true, N, replace = TRUE)
  lambda_mat <- matrix(rgamma(K_true * D, shape = alpha, rate = beta), nrow = K_true)
  
  if (noise_features > 0) {
    noise_idx <- (true_features + 1):(true_features + noise_features)
    shared_lambda <- rgamma(noise_features, shape = alpha, rate = beta)
    lambda_mat[, noise_idx] <- matrix(rep(shared_lambda, each = K_true), nrow = K_true)
  }
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    X[n, ] <- rpois(D, lambda = lambda_mat[z[n], ])
  }
  
  list(X = X, z = z)
}

# Custom EM Algorithm with BIC-based K selection
run_my_em_autoK <- function(X, K_range, max_iter = 100, tol = 1e-4, restarts = EM_INIT) {
  N <- nrow(X)
  D <- ncol(X)
  best_bic <- Inf
  best_result <- NULL
  
  for (K in K_range) {
    best_ll <- -Inf
    best_model <- NULL
    
    for (r in 1:restarts) {
      pi_k <- rep(1/K, K)
      lambda_kd <- matrix(rgamma(K * D, shape = alpha, rate = beta), nrow = K)
      ll_old <- -Inf
      
      for (iter in 1:max_iter) {
        log_resp <- matrix(0, N, K)
        for (k in 1:K) {
          log_resp[, k] <- log(pi_k[k]) + rowSums(dpois(X, lambda = matrix(lambda_kd[k, ], N, D, byrow = TRUE), log = TRUE))
        }
        log_sum <- apply(log_resp, 1, logsumexp)
        log_resp <- sweep(log_resp, 1, log_sum, "-")
        resp <- exp(log_resp)
        
        N_k <- colSums(resp)
        if (any(N_k == 0)) break  # avoid division by zero
        
        pi_k <- N_k / N
        lambda_kd <- (t(resp) %*% X) / N_k
        
        ll_new <- sum(log_sum)
        if (abs(ll_new - ll_old) < tol) break
        ll_old <- ll_new
      }
      
      if (!is.nan(ll_new) && ll_new > best_ll) {
        best_ll <- ll_new
        cluster_assignments <- apply(resp, 1, which.max)
        actual_K <- length(unique(cluster_assignments))
        best_model <- list(K = K, clusters = cluster_assignments, ll = ll_new, pi = pi_k, lambda = lambda_kd, K_actual = actual_K)
      }
    }
    
    if (!is.null(best_model)) {
      bic <- -2 * best_model$ll + (K * D + K - 1) * log(N)
      if (bic < best_bic) {
        best_bic <- bic
        best_result <- best_model
      }
    }
  }
  
  return(best_result)
}

# FlexMix with BIC and filtered non-empty clusters
run_flexmix_autoK <- function(X, K_range) {
  data_df <- data.frame(X)
  models <- list()
  bics <- numeric(length(K_range))
  
  for (i in seq_along(K_range)) {
    k <- K_range[i]
    model <- tryCatch({
      flexmix(X ~ 1, data = data_df, k = k, model = FLXMCmvpois(), control = list(verbose = 0))
    }, error = function(e) NULL)
    
    models[[i]] <- model
    bics[i] <- if (!is.null(model)) BIC(model) else Inf
  }
  
  best_idx <- which.min(bics)
  best_model <- models[[best_idx]]
  
  if (!is.null(best_model)) {
    cl <- clusters(best_model)
    actual_K <- length(unique(cl[cl > 0]))
    return(list(clusters = cl, K = K_range[best_idx], K_actual = actual_K))
  } else {
    return(NULL)
  }
}

# Experiment loop
noise_levels <- seq(0, 95, by = 5)
results <- data.frame(
  noise_features = noise_levels,
  ARI_myEM = NA,
  ARI_flexmix = NA,
  time_myEM = NA,
  time_flexmix = NA,
  K_myEM = NA,
  K_myEM_actual = NA,
  K_flexmix = NA,
  K_flexmix_actual = NA
)

pb <- txtProgressBar(min = 0, max = length(noise_levels), style = 3)

for (i in seq_along(noise_levels)) {
  nf <- noise_levels[i]
  cat("\nNoise Features:", nf, "\n")
  data <- generate_dataset(nf)
  X <- data$X
  z_true <- data$z
  
  # My EM
  t1 <- Sys.time()
  em_result <- run_my_em_autoK(X, K_range)
  t2 <- Sys.time()
  if (!is.null(em_result)) {
    results$ARI_myEM[i] <- adjustedRandIndex(em_result$clusters, z_true)
    results$time_myEM[i] <- as.numeric(difftime(t2, t1, units = "secs"))
    results$K_myEM[i] <- em_result$K
    results$K_myEM_actual[i] <- em_result$K_actual
  }
  
  # FlexMix
  t1 <- Sys.time()
  fx_result <- run_flexmix_autoK(X, K_range)
  t2 <- Sys.time()
  if (!is.null(fx_result)) {
    results$ARI_flexmix[i] <- adjustedRandIndex(fx_result$clusters, z_true)
    results$time_flexmix[i] <- as.numeric(difftime(t2, t1, units = "secs"))
    results$K_flexmix[i] <- fx_result$K
    results$K_flexmix_actual[i] <- fx_result$K_actual
  }
  
  # Update progress bar
  setTxtProgressBar(pb, i)
}

close(pb)  # Close progress bar after loop


# Results table
print(results)

# Plot ARI
ggplot(results, aes(x = noise_features)) +
  geom_line(aes(y = ARI_myEM, color = "My EM")) +
  geom_line(aes(y = ARI_flexmix, color = "FlexMix")) +
  labs(title = "ARI vs. Noise", x = "Noisy features", y = "ARI") +
  scale_color_manual(values = c("My EM" = "blue", "FlexMix" = "red")) +
  theme_minimal()

# Plot Runtime
ggplot(results, aes(x = noise_features)) +
  geom_line(aes(y = time_myEM, color = "My EM")) +
  geom_line(aes(y = time_flexmix, color = "FlexMix")) +
  labs(title = "Runtime vs. Noise", x = "Noisy features", y = "Time (s)") +
  scale_color_manual(values = c("My EM" = "blue", "FlexMix" = "red")) +
  theme_minimal()

# Plot Detected K (non-empty clusters only)
ggplot(results, aes(x = noise_features)) +
  geom_line(aes(y = as.numeric(K_myEM_actual), color = "My EM")) +
  geom_point(aes(y = as.numeric(K_myEM_actual), color = "My EM")) +
  geom_line(aes(y = as.numeric(K_flexmix_actual), color = "FlexMix (non-empty clusters)")) +
  geom_point(aes(y = as.numeric(K_flexmix_actual), color = "FlexMix (non-empty clusters)")) +
  labs(
    title = "Selected Number of Clusters vs. Noise",
    x = "Noisy features",
    y = "Detected Clusters"
  ) +
  scale_color_manual(values = c("My EM" = "blue", "FlexMix (non-empty clusters)" = "red")) +
  theme_minimal()
##########################################################################################################################################################

library(ggplot2)
library(pheatmap)
library(mclust)

set.seed(42)

# Parameters
N <- 1000
K_true <- 5
true_features <- 5
noise_features <- 50
alpha <- 5
beta <- 1
runs <- 50

# Data Generator Function
generate_dataset <- function(noise_features) {
  D <- true_features + noise_features
  z <- sample(1:K_true, N, replace = TRUE)
  lambda_mat <- matrix(rgamma(K_true * D, shape = alpha, rate = beta), nrow = K_true)
  
  if (noise_features > 0) {
    noise_idx <- (true_features + 1):(true_features + noise_features)
    shared_lambda <- rgamma(noise_features, shape = alpha, rate = beta)
    lambda_mat[, noise_idx] <- matrix(rep(shared_lambda, each = K_true), nrow = K_true)
  }
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    X[n, ] <- rpois(D, lambda = lambda_mat[z[n], ])
  }
  
  colnames(X) <- c(paste0("Informative_", 1:true_features),
                   paste0("Noise_", 1:noise_features))
  
  df <- as.data.frame(X)
  df$TrueCluster <- z
  return(df)
}

# Generate and Save One Dataset for Export
df_one <- generate_dataset(noise_features)
write.csv(df_one, "~/Desktop/poisson_dataset_50_noise.csv", row.names = FALSE)

# Save True Cluster Labels Separately
write.csv(data.frame(TrueCluster = df_one$TrueCluster), "~/Desktop/true_labels_50_noise.csv", row.names = FALSE)

# Visual Heatmaps
X_only <- df_one[, 1:(true_features + noise_features)]

# Mixed (random row order)
pheatmap(X_only[sample(nrow(X_only)), ], main = "Mixed Heatmap (Shuffled Rows)")

# Clustered (sorted by TrueCluster)
sorted_idx <- order(df_one$TrueCluster)
pheatmap(X_only[sorted_idx, ], main = "Clustered Heatmap (Sorted by True Cluster)")

# 5 Runs for Variance Analysis
results <- data.frame(Run = 1:runs, ARI = NA)

for (i in 1:runs) {
  set.seed(100 + i)
  df <- generate_dataset(noise_features)
  z_true <- df$TrueCluster
  
  # Use Gaussian Mixture on informative features only
  model <- Mclust(df[, 1:true_features], G = K_true)
  z_pred <- model$classification
  
  results$ARI[i] <- adjustedRandIndex(z_true, z_pred)
}

# Boxplot of ARI over 5 runs
ggplot(results, aes(x = "ARI over 50 runs", y = ARI)) +
  geom_boxplot(fill = "skyblue") +
  geom_jitter(width = 0.1) +
  labs(title = "Clustering Stability over 50 Simulations", x = "", y = "Adjusted Rand Index") +
  theme_minimal()



library(pheatmap)

# Load or generate your dataset
set.seed(42)
df <- generate_dataset(50)  # 5 informative + 50 noise features
X <- df[, 1:(true_features + noise_features)]

# Mixed Heatmap (shuffled rows, no clustering)
pheatmap(X[sample(nrow(X)), ],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Mixed Heatmap (Shuffled Rows)")

# Clustered Heatmap (sorted by true cluster, no clustering)
sorted_idx <- order(df$TrueCluster)
pheatmap(X[sorted_idx, ],
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Clustered Heatmap (True Cluster Order)")

################################################################################################################

# CHATGPT REFINED VERSION


library(flexmix)
library(mclust)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

# Parameters
N <- 1000
K_true <- 5
true_features <- 5
alpha <- 5
beta <- 1
EM_INIT <- 15
K_range <- 1:10
REPEATS <- 5  # Number of datasets per noise level

# Log-Sum-Exp Trick
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Synthetic Dataset Generator
generate_dataset <- function(noise_features) {
  D <- true_features + noise_features
  z <- sample(1:K_true, N, replace = TRUE)
  lambda_mat <- matrix(rgamma(K_true * D, shape = alpha, rate = beta), nrow = K_true)
  
  if (noise_features > 0) {
    noise_idx <- (true_features + 1):(true_features + noise_features)
    shared_lambda <- rgamma(noise_features, shape = alpha, rate = beta)
    lambda_mat[, noise_idx] <- matrix(rep(shared_lambda, each = K_true), nrow = K_true)
  }
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    X[n, ] <- rpois(D, lambda = lambda_mat[z[n], ])
  }
  
  list(X = X, z = z)
}

# Custom EM
run_my_em_autoK <- function(X, K_range, max_iter = 100, tol = 1e-4, restarts = EM_INIT) {
  N <- nrow(X)
  D <- ncol(X)
  best_bic <- Inf
  best_result <- NULL
  
  for (K in K_range) {
    best_ll <- -Inf
    best_model <- NULL
    
    for (r in 1:restarts) {
      pi_k <- rep(1/K, K)
      lambda_kd <- matrix(rgamma(K * D, shape = alpha, rate = beta), nrow = K)
      ll_old <- -Inf
      
      for (iter in 1:max_iter) {
        log_resp <- matrix(0, N, K)
        for (k in 1:K) {
          log_resp[, k] <- log(pi_k[k]) + rowSums(dpois(X, lambda = matrix(lambda_kd[k, ], N, D, byrow = TRUE), log = TRUE))
        }
        log_sum <- apply(log_resp, 1, logsumexp)
        log_resp <- sweep(log_resp, 1, log_sum, "-")
        resp <- exp(log_resp)
        
        N_k <- colSums(resp)
        if (any(N_k == 0)) break
        
        pi_k <- N_k / N
        lambda_kd <- (t(resp) %*% X) / N_k
        ll_new <- sum(log_sum)
        if (abs(ll_new - ll_old) < tol) break
        ll_old <- ll_new
      }
      
      if (!is.nan(ll_new) && ll_new > best_ll) {
        best_ll <- ll_new
        cluster_assignments <- apply(resp, 1, which.max)
        actual_K <- length(unique(cluster_assignments))
        best_model <- list(K = K, clusters = cluster_assignments, ll = ll_new, pi = pi_k, lambda = lambda_kd, K_actual = actual_K)
      }
    }
    
    if (!is.null(best_model)) {
      bic <- -2 * best_model$ll + (K * D + K - 1) * log(N)
      if (bic < best_bic) {
        best_bic <- bic
        best_result <- best_model
      }
    }
  }
  
  return(best_result)
}

# FlexMix
run_flexmix_autoK <- function(X, K_range) {
  data_df <- data.frame(X)
  models <- list()
  bics <- numeric(length(K_range))
  
  for (i in seq_along(K_range)) {
    k <- K_range[i]
    model <- tryCatch({
      flexmix(X ~ 1, data = data_df, k = k, model = FLXMCmvpois(), control = list(verbose = 0))
    }, error = function(e) NULL)
    
    models[[i]] <- model
    bics[i] <- if (!is.null(model)) BIC(model) else Inf
  }
  
  best_idx <- which.min(bics)
  best_model <- models[[best_idx]]
  
  if (!is.null(best_model)) {
    cl <- clusters(best_model)
    actual_K <- length(unique(cl[cl > 0]))
    return(list(clusters = cl, K = K_range[best_idx], K_actual = actual_K))
  } else {
    return(NULL)
  }
}

# Run Experiments with 5 replicates per noise level
noise_levels <- seq(0, 95, by = 5)
all_results <- data.frame()

for (nf in noise_levels) {
  for (rep in 1:REPEATS) {
    cat("Noise:", nf, "| Replicate:", rep, "\n")
    data <- generate_dataset(nf)
    X <- data$X
    z_true <- data$z
    
    # Save true labels for professor
    if (nf == 0 && rep == 1) {
      write.csv(data.frame(z_true), "true_labels_first_run.csv", row.names = FALSE)
    }
    
    # My EM
    t1 <- Sys.time()
    em_result <- run_my_em_autoK(X, K_range)
    t2 <- Sys.time()
    ARI_myEM <- if (!is.null(em_result)) adjustedRandIndex(em_result$clusters, z_true) else NA
    time_myEM <- as.numeric(difftime(t2, t1, units = "secs"))
    K_myEM <- if (!is.null(em_result)) em_result$K else NA
    K_myEM_actual <- if (!is.null(em_result)) em_result$K_actual else NA
    
    # FlexMix
    t1 <- Sys.time()
    fx_result <- run_flexmix_autoK(X, K_range)
    t2 <- Sys.time()
    ARI_fx <- if (!is.null(fx_result)) adjustedRandIndex(fx_result$clusters, z_true) else NA
    time_fx <- as.numeric(difftime(t2, t1, units = "secs"))
    K_fx <- if (!is.null(fx_result)) fx_result$K else NA
    K_fx_actual <- if (!is.null(fx_result)) fx_result$K_actual else NA
    
    all_results <- rbind(all_results, data.frame(
      noise_features = nf,
      replicate = rep,
      ARI_myEM = ARI_myEM,
      ARI_flexmix = ARI_fx,
      time_myEM = time_myEM,
      time_flexmix = time_fx,
      K_myEM = K_myEM,
      K_myEM_actual = K_myEM_actual,
      K_flexmix = K_fx,
      K_flexmix_actual = K_fx_actual
    ))
  }
}

# Save results
write.csv(all_results, "all_results_boxplot_ready.csv", row.names = FALSE)

# Plot Boxplots
ggplot(all_results, aes(x = factor(noise_features), y = ARI_myEM)) +
  geom_boxplot(fill = "lightblue") +
  labs(title = "ARI of My EM vs Noise (Boxplot)", x = "Noise Features", y = "ARI") +
  theme_minimal()

ggplot(all_results, aes(x = factor(noise_features), y = ARI_flexmix)) +
  geom_boxplot(fill = "salmon") +
  labs(title = "ARI of FlexMix vs Noise (Boxplot)", x = "Noise Features", y = "ARI") +
  theme_minimal()




