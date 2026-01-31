# Required Libraries
library(ggplot2)
library(pheatmap)
library(mclust)
library(flexmix)
library(gridExtra)

set.seed(42)

# === Gamma-Poisson Illustration ===
alpha <- 5
beta <- 1
lambda_sample <- rgamma(1, shape = alpha, rate = beta)

# Plot Gamma(5,1) PDF with Sampled Lambda
x_vals <- seq(0.01, 25, length.out = 500)
gamma_df <- data.frame(
  x = x_vals,
  density = dgamma(x_vals, shape = alpha, rate = beta)
)

p1 <- ggplot(gamma_df, aes(x, density)) +
  geom_line(color = "blue", size = 1.2) +
  geom_vline(xintercept = lambda_sample, linetype = "dashed", color = "red", linewidth = 1) +
  labs(title = "Gamma(5,1) PDF with Sampled λ",
       x = "λ", y = "Density")

# Poisson(λ) PMF
k_vals <- 0:25
pois_df <- data.frame(
  x = k_vals,
  prob = dpois(k_vals, lambda = lambda_sample)
)

p2 <- ggplot(pois_df, aes(x, prob)) +
  geom_bar(stat = "identity", fill = "orange", color = "black") +
  labs(title = paste0("Poisson(λ ≈ ", round(lambda_sample, 2), ") PMF"),
       x = "x", y = "P(X = x)")

grid.arrange(p1, p2, ncol = 2)

# === Data Generator ===
generate_dataset <- function(noise_features = 0, true_features = 5, N = 1000, K = 5, alpha = 5, beta = 1) {
  D <- true_features + noise_features
  z <- sample(1:K, N, replace = TRUE)
  lambda_mat <- matrix(rgamma(K * D, shape = alpha, rate = beta), nrow = K)
  
  if (noise_features > 0) {
    noise_idx <- (true_features + 1):(true_features + noise_features)
    shared_lambda <- rgamma(noise_features, shape = alpha, rate = beta)
    lambda_mat[, noise_idx] <- matrix(rep(shared_lambda, each = K), nrow = K)
  }
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    X[n, ] <- rpois(D, lambda = lambda_mat[z[n], ])
  }
  
  df <- as.data.frame(X)
  colnames(df) <- c(paste0("Informative_", 1:true_features), paste0("Noise_", 1:noise_features))
  df$TrueCluster <- z
  return(df)
}

# === Heatmap with 0 Noise Features ===
df_clean <- generate_dataset(noise_features = 0)
X_clean <- df_clean[, 1:5]  # only informative
sorted_idx <- order(df_clean$TrueCluster)

# Heatmap (rows sorted by true cluster)
pheatmap(X_clean[sorted_idx, ],
         cluster_rows = FALSE, cluster_cols = FALSE,
         main = "Heatmap: Only Informative Features (Sorted by True Cluster)")

# === Run Clustering and Compute ARI ===
z_true <- df_clean$TrueCluster

# EM via mclust
model <- Mclust(X_clean, G = 5)
z_pred <- model$classification

# ARI: Adjusted Rand Index (measures similarity of predicted clustering to true labels)
# Ranges from -1 (completely wrong) to 1 (perfect match). 0 = random clustering.
ari_val <- adjustedRandIndex(z_true, z_pred)

cat("Adjusted Rand Index (ARI) on 0-noise dataset:", round(ari_val, 3), "\n")
