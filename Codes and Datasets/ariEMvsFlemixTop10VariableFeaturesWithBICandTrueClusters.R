#-------------------- Libraries --------------------#
library(flexmix)
library(mclust)
library(ggplot2)
library(tidyr)

set.seed(42)

#-------------------- Parameters --------------------#
N <- 1000
K_true <- 5
true_features <- 5
noise_features <- 95
D <- true_features + noise_features
alpha <- 5
beta <- 1
K_range <- 1:10

logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#-------------------- Dataset Generator --------------------#
generate_dataset <- function() {
  z <- sample(1:K_true, N, replace = TRUE)
  lambda_mat <- matrix(rgamma(K_true * D, shape = alpha, rate = beta), nrow = K_true)
  
  # Shared noise features
  noise_idx <- (true_features + 1):D
  shared_lambda <- rgamma(noise_features, shape = alpha, rate = beta)
  lambda_mat[, noise_idx] <- matrix(rep(shared_lambda, each = K_true), nrow = K_true)
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    X[n, ] <- rpois(D, lambda = lambda_mat[z[n], ])
  }
  
  list(X = X, z = z)
}

#-------------------- Your Poisson EM --------------------#
run_my_em <- function(X, K, max_iter = 100, tol = 1e-4) {
  N <- nrow(X)
  D <- ncol(X)
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
  
  clusters <- apply(resp, 1, which.max)
  list(ll = ll_new, clusters = clusters, pi = pi_k, lambda = lambda_kd)
}

#-------------------- EM with BIC --------------------#
run_my_em_bic <- function(X, K_range) {
  N <- nrow(X)
  D <- ncol(X)
  best_bic <- Inf
  best_result <- NULL
  best_K <- NULL
  
  for (K in K_range) {
    result <- run_my_em(X, K)
    bic <- -2 * result$ll + (K * D + K - 1) * log(N)
    if (bic < best_bic) {
      best_bic <- bic
      best_result <- result
      best_K <- K
    }
  }
  
  best_result$K <- best_K
  return(best_result)
}

#-------------------- FlexMix fixed K --------------------#
run_flexmix_K <- function(X, K) {
  data_df <- data.frame(X)
  model <- tryCatch({
    flexmix(X ~ 1, data = data_df, k = K, model = FLXMCmvpois(), control = list(verbose = 0))
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    cl <- clusters(model)
    return(list(clusters = cl, K = K))
  } else {
    return(NULL)
  }
}

#-------------------- FlexMix with BIC --------------------#
run_flexmix_bic <- function(X, K_range) {
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
    return(list(clusters = cl, K = K_range[best_idx]))
  } else {
    return(NULL)
  }
}

#-------------------- Run Full Experiment --------------------#
data <- generate_dataset()
X <- data$X
z_true <- data$z

# Feature selection: Top 10 most variable features
variances <- apply(X, 2, var)
top_10_idx <- order(variances, decreasing = TRUE)[1:10]
X_top10 <- X[, top_10_idx]

# Count how many informative features are present
informative_in_top10 <- sum(1:true_features %in% top_10_idx)
cat("Top 10 highly variable features include", informative_in_top10, "true informative features\n")

# Run BIC models once
em_bic_model <- run_my_em_bic(X_top10, K_range)
ari_em_bic <- adjustedRandIndex(em_bic_model$clusters, z_true)
em_bic_K <- em_bic_model$K
cat("EM BIC selected K:", em_bic_K, "\n")

fx_bic_model <- run_flexmix_bic(X_top10, K_range)
ari_fx_bic <- if (!is.null(fx_bic_model)) adjustedRandIndex(fx_bic_model$clusters, z_true) else NA
fx_bic_K <- if (!is.null(fx_bic_model)) fx_bic_model$K else NA
cat("FlexMix BIC selected K:", fx_bic_K, "\n")

# Evaluate models for all K in K_range
results_df <- data.frame()

for (K in K_range) {
  em_fixed <- run_my_em(X_top10, 5)
  ari_em_fixed <- adjustedRandIndex(em_fixed$clusters, z_true)
  
  fx_fixed <- run_flexmix_K(X_top10, 5)
  ari_fx_fixed <- if (!is.null(fx_fixed)) adjustedRandIndex(fx_fixed$clusters, z_true) else NA
  
  results_df <- rbind(results_df, data.frame(
    K = K,
    EM_K5 = ari_em_fixed,
    EM_BIC = ari_em_bic,
    FlexMix_K5 = ari_fx_fixed,
    FlexMix_BIC = ari_fx_bic
  ))
}

#-------------------- Plot --------------------#
results_long <- pivot_longer(results_df, cols = -K, names_to = "Model", values_to = "ARI")

ggplot(results_long, aes(x = K, y = ARI, color = Model, group = Model)) +
  geom_line(size = 1.1) +
  geom_point(size = 2.2) +
  labs(
    title = "ARI Comparison vs K using Top 10 Variable Features",
    subtitle = paste(
      informative_in_top10, "true informative features in Top 10 |",
      "EM_BIC chose K =", em_bic_K, "|",
      "FlexMix_BIC chose K =", fx_bic_K
    ),
    x = "K (number of clusters)",
    y = "Adjusted Rand Index (ARI)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14)
##################################################################################################################################
















#==================== Libraries =====================#
library(flexmix)
library(aricode)  # for ARI
library(ggplot2)
library(tidyr)
library(progress)
library(umap)

set.seed(42)

#==================== Parameters ====================#
N <- 1000
K_true <- 5
true_features <- 5
noise_features <- 95
D <- true_features + noise_features
alpha <- 5
beta <- 1
K_range <- 1:10
RESTARTS <- 50

#==================== LogSumExp =====================#
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#==================== Dataset Generator =============#
generate_dataset <- function() {
  z <- sample(1:K_true, N, replace = TRUE)
  lambda_mat <- matrix(rgamma(K_true * D, shape = alpha, rate = beta), nrow = K_true)
  
  noise_idx <- (true_features + 1):D
  shared_lambda <- rgamma(noise_features, shape = alpha, rate = beta)
  lambda_mat[, noise_idx] <- matrix(rep(shared_lambda, each = K_true), nrow = K_true)
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    X[n, ] <- rpois(D, lambda = lambda_mat[z[n], ])
  }
  list(X = X, z = z)
}

#==================== Poisson EM =====================#
run_my_em <- function(X, K, max_iter = 100, tol = 1e-4, restarts = RESTARTS) {
  N <- nrow(X)
  D <- ncol(X)
  best_ll <- -Inf
  best_result <- NULL
  
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
      best_result <- list(ll = ll_new, clusters = apply(resp, 1, which.max), pi = pi_k, lambda = lambda_kd)
    }
  }
  return(best_result)
}

#==================== EM BIC =====================#
run_my_em_bic <- function(X, K_range) {
  N <- nrow(X)
  D <- ncol(X)
  best_bic <- Inf
  best_result <- NULL
  best_K <- NULL
  
  for (K in K_range) {
    result <- run_my_em(X, K)
    bic <- -2 * result$ll + (K * D + K - 1) * log(N)
    if (bic < best_bic) {
      best_bic <- bic
      best_result <- result
      best_K <- K
    }
  }
  
  best_result$K <- best_K
  return(best_result)
}

#==================== FlexMix =====================#
run_flexmix_K <- function(X, K) {
  data_df <- data.frame(X)
  model <- tryCatch({
    flexmix(X ~ 1, data = data_df, k = K, model = FLXMCmvpois(), control = list(verbose = 0))
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    cl <- clusters(model)
    return(list(clusters = cl, K = K))
  } else {
    return(NULL)
  }
}

run_flexmix_bic <- function(X, K_range) {
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
    return(list(clusters = cl, K = K_range[best_idx]))
  } else {
    return(NULL)
  }
}

#==================== Run Experiment =====================#
data <- generate_dataset()
X <- data$X
z_true <- data$z

variances <- apply(X, 2, var)
top_10_idx <- order(variances, decreasing = TRUE)[1:10]
X_top10 <- X[, top_10_idx]
informative_in_top10 <- sum(1:true_features %in% top_10_idx)
cat("Top 10 highly variable features include", informative_in_top10, "true informative features\n")

results_df <- data.frame()
pb <- progress_bar$new(total = length(K_range), format = "Progress [:bar] :percent ETA: :eta")

for (K in K_range) {
  pb$tick()
  
  em_fixed <- run_my_em(X_top10, 5)
  ari_em_fixed <- ARI(em_fixed$clusters, z_true)
  
  em_bic_model <- run_my_em_bic(X_top10, K_range)
  ari_em_bic <- ARI(em_bic_model$clusters, z_true)
  
  fx_fixed <- run_flexmix_K(X_top10, 5)
  ari_fx_fixed <- if (!is.null(fx_fixed)) ARI(fx_fixed$clusters, z_true) else NA
  
  fx_bic <- run_flexmix_bic(X_top10, K_range)
  ari_fx_bic <- if (!is.null(fx_bic)) ARI(fx_bic$clusters, z_true) else NA
  
  results_df <- rbind(results_df, data.frame(
    K = K,
    EM_K5 = ari_em_fixed,
    EM_BIC = ari_em_bic,
    FlexMix_K5 = ari_fx_fixed,
    FlexMix_BIC = ari_fx_bic
  ))
}

#==================== Plot ARI vs K =====================#
results_long <- pivot_longer(results_df, cols = -K, names_to = "Model", values_to = "ARI")

ggplot(results_long, aes(x = K, y = ARI, color = Model, group = Model)) +
  geom_line(size = 1.1) +
  geom_point(size = 2.2) +
  labs(
    title = "ARI Comparison vs K using Top 10 Variable Features",
    subtitle = paste(informative_in_top10, "true informative features in Top 10 | EM_BIC chose K =", em_bic_model$K,
                     "| FlexMix_BIC chose K =", fx_bic$K),
    x = "K (number of clusters)",
    y = "Adjusted Rand Index (ARI)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14)

#==================== UMAP Visualization =====================#
umap_out <- umap(X_top10)$layout
umap_df <- data.frame(
  UMAP1 = umap_out[,1],
  UMAP2 = umap_out[,2],
  True = as.factor(z_true),
  EM_K5 = as.factor(em_fixed$clusters),
  EM_BIC = as.factor(em_bic_model$clusters),
  FlexMix_K5 = as.factor(fx_fixed$clusters),
  FlexMix_BIC = as.factor(fx_bic$clusters)
)

plot_cluster <- function(label, title) {
  ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = .data[[label]])) +
    geom_point(alpha = 0.6, size = 1.2) +
    labs(title = title, color = "Cluster") +
    theme_minimal(base_size = 14)
}

plot_cluster("True", "True Clusters")
plot_cluster("EM_K5", "EM (K=5)")
plot_cluster("EM_BIC", paste0("EM (BIC K=", em_bic_model$K, ")"))
plot_cluster("FlexMix_K5", "FlexMix (K=5)")
plot_cluster("FlexMix_BIC", paste0("FlexMix (BIC K=", fx_bic$K, ")"))


#___________________________________________________________________________________________________

# re running the above code with increased sample size and changing parameters

#==================== Libraries =====================#
library(flexmix)
library(aricode)  # for ARI
library(ggplot2)
library(tidyr)
library(progress)

set.seed(43)

#==================== Parameters =====================#
N <- 1000
K_true <- 5
true_features <- 5
noise_features <- 95
D <- true_features + noise_features
alpha <- 5
beta <- 1
K_range <- 1:10
RESTARTS <- 50

#==================== LogSumExp =====================#
logsumexp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#==================== Dataset Generator =============#
generate_dataset <- function() {
  z <- sample(1:K_true, N, replace = TRUE)
  lambda_mat <- matrix(rgamma(K_true * D, shape = alpha, rate = beta), nrow = K_true)
  
  noise_idx <- (true_features + 1):D
  shared_lambda <- rgamma(noise_features, shape = alpha, rate = beta)
  lambda_mat[, noise_idx] <- matrix(rep(shared_lambda, each = K_true), nrow = K_true)
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    X[n, ] <- rpois(D, lambda = lambda_mat[z[n], ])
  }
  
  list(X = X, z = z)
}

#==================== Poisson EM =====================#
run_my_em <- function(X, K, max_iter = 100, tol = 1e-6, restarts = RESTARTS) {
  N <- nrow(X)
  D <- ncol(X)
  best_ll <- -Inf
  best_result <- NULL
  
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
      best_result <- list(ll = ll_new, clusters = apply(resp, 1, which.max), pi = pi_k, lambda = lambda_kd)
    }
  }
  
  return(best_result)
}

#==================== EM BIC =====================#
run_my_em_bic <- function(X, K_range) {
  N <- nrow(X)
  D <- ncol(X)
  best_bic <- Inf
  best_result <- NULL
  best_K <- NULL
  
  for (K in K_range) {
    result <- run_my_em(X, K)
    bic <- -2 * result$ll + (K * D + K - 1) * log(N)
    if (bic < best_bic) {
      best_bic <- bic
      best_result <- result
      best_K <- K
    }
  }
  
  best_result$K <- best_K
  return(best_result)
}

#==================== FlexMix =====================#
run_flexmix_K <- function(X, K) {
  data_df <- data.frame(X)
  model <- tryCatch({
    flexmix(X ~ 1, data = data_df, k = K, model = FLXMCmvpois(), control = list(verbose = 0))
  }, error = function(e) NULL)
  
  if (!is.null(model)) {
    cl <- clusters(model)
    return(list(clusters = cl, K = K))
  } else {
    return(NULL)
  }
}

run_flexmix_bic <- function(X, K_range) {
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
    return(list(clusters = cl, K = K_range[best_idx]))
  } else {
    return(NULL)
  }
}

#==================== Run Experiment =====================#
data <- generate_dataset()
X <- data$X
z_true <- data$z

variances <- apply(X, 2, var)
top_20_idx <- order(variances, decreasing = TRUE)[1:10]
X_top20 <- X[, top_20_idx]
informative_in_top20 <- sum(1:true_features %in% top_20_idx)
cat("Top 20 highly variable features include", informative_in_top20, "true informative features\n")

# BIC-based models (run ONCE)
em_bic_model <- run_my_em_bic(X_top20, K_range)
ari_em_bic <- ARI(em_bic_model$clusters, z_true)

fx_bic <- run_flexmix_bic(X_top20, K_range)
ari_fx_bic <- if (!is.null(fx_bic)) ARI(fx_bic$clusters, z_true) else NA

results_df <- data.frame()
pb <- progress_bar$new(total = length(K_range), format = "Progress [:bar] :percent ETA: :eta")

# Loop through K only for fixed-K runs
for (K in K_range) {
  pb$tick()
  
  em_fixed <- run_my_em(X_top20, 5)
  ari_em_fixed <- ARI(em_fixed$clusters, z_true)
  
  fx_fixed <- run_flexmix_K(X_top20, 5)
  ari_fx_fixed <- if (!is.null(fx_fixed)) ARI(fx_fixed$clusters, z_true) else NA
  
  results_df <- rbind(results_df, data.frame(
    K = K,
    EM_K5 = ari_em_fixed,
    EM_BIC = ari_em_bic,          # same every row
    FlexMix_K5 = ari_fx_fixed,
    FlexMix_BIC = ari_fx_bic      # same every row
  ))
}

#==================== Plot ARI vs K =====================#
results_long <- pivot_longer(results_df, cols = -K, names_to = "Model", values_to = "ARI")

ggplot(results_long, aes(x = K, y = ARI, color = Model, group = Model)) +
  geom_line(size = 1.1) +
  geom_point(size = 2.2) +
  labs(
    title = "ARI Comparison vs K using Top 10 Variable Features",
    subtitle = paste(informative_in_top20, "true informative features in Top 10 | EM_BIC chose K =", em_bic_model$K,
                     "| FlexMix_BIC chose K =", fx_bic$K),
    x = "K (number of clusters)",
    y = "Adjusted Rand Index (ARI)",
    color = "Model"
  ) +
  theme_minimal(base_size = 14)


