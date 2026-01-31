# Load libraries
library(MCMCpack)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(tidyr)
library(dplyr)

set.seed(123)

# Step 1: Generate Data
generate_data <- function(N, K, alpha0, beta0) {
  pi_true <- rdirichlet(1, rep(1, K))[1, ]
  lambda_true <- rgamma(K, shape = alpha0, rate = beta0)
  z <- sample(1:K, N, replace = TRUE, prob = pi_true)
  x <- rpois(N, lambda = lambda_true[z])
  return(list(x = x, z = z, pi = pi_true, lambda = lambda_true))
}

# Parameters
K <- 4
N <- 500
alpha0 <- 2
beta0 <- 0.0001
data <- generate_data(N, K, alpha0, beta0)

# Step 2: EM Algorithm with safe convergence
EM_poisson <- function(x, K = 3, max_iter = 100, tol = 1e-6) {
  N <- length(x)
  pi_k <- rep(1/K, K)
  lambda_k <- quantile(x, probs = seq(0.2, 0.8, length.out = K))
  ll_old <- -Inf
  
  for (iter in 1:max_iter) {
    gamma_nk <- sapply(1:K, function(k) pi_k[k] * dpois(x, lambda_k[k]))
    gamma_nk <- gamma_nk / rowSums(gamma_nk)
    
    Nk <- colSums(gamma_nk)
    pi_k <- Nk / N
    lambda_k <- colSums(gamma_nk * x) / Nk
    
    ll_new <- sum(log(rowSums(sapply(1:K, function(k) pi_k[k] * dpois(x, lambda_k[k])))))
    
    cat(sprintf("EM Iteration %d: Log-Likelihood = %.6f\n", iter, ll_new))
    
    if (abs(ll_new - ll_old) < tol) break
    ll_old <- ll_new
  }
  return(list(pi = pi_k, lambda = lambda_k))
}

em_fit <- EM_poisson(data$x, K)

# Step 3: Gibbs Sampling
gibbs_poisson <- function(x, K = 3, alpha0 = 2, beta0 = 0.01, iter = 1000) {
  N <- length(x)
  z <- sample(1:K, N, replace = TRUE)
  lambda_samples <- matrix(0, nrow = iter, ncol = K)
  pi_samples <- matrix(0, nrow = iter, ncol = K)
  
  for (it in 1:iter) {
    for (k in 1:K) {
      nk <- sum(z == k)
      sum_xk <- sum(x[z == k])
      lambda_samples[it, k] <- rgamma(1, shape = alpha0 + sum_xk, rate = beta0 + nk)
    }
    nk <- table(factor(z, levels = 1:K))
    pi_samples[it, ] <- rdirichlet(1, rep(1, K) + nk)
    
    for (n in 1:N) {
      probs <- pi_samples[it, ] * dpois(x[n], lambda_samples[it, ])
      probs <- probs / sum(probs)
      z[n] <- sample(1:K, 1, prob = probs)
    }
  }
  
  return(list(
    lambda = colMeans(lambda_samples[500:iter, ]),
    pi = colMeans(pi_samples[500:iter, ])
  ))
}

gibbs_fit <- gibbs_poisson(data$x, K, alpha0, beta0)

# Step 4: Distribution Comparison
x_vals <- 0:max(data$x)
df_dens <- data.frame(
  x = x_vals,
  True = rowSums(sapply(1:K, function(k) data$pi[k] * dpois(x_vals, data$lambda[k]))),
  EM = rowSums(sapply(1:K, function(k) em_fit$pi[k] * dpois(x_vals, em_fit$lambda[k]))),
  Gibbs = rowSums(sapply(1:K, function(k) gibbs_fit$pi[k] * dpois(x_vals, gibbs_fit$lambda[k])))
)

df_dens_long <- df_dens %>%
  pivot_longer(cols = -x, names_to = "Method", values_to = "Density")

p1 <- ggplot() +
  geom_histogram(data = data.frame(x = data$x), aes(x = x, y = ..density..),
                 binwidth = 5, fill = "lightblue", alpha = 0.4, boundary = 0) +
  geom_line(data = df_dens_long, aes(x = x, y = Density, color = Method, linetype = Method), size = 1) +
  scale_color_manual(values = c("True" = "blue", "EM" = "red", "Gibbs" = "darkgreen")) +
  scale_linetype_manual(values = c("True" = "solid", "EM" = "dashed", "Gibbs" = "dotted")) +
  labs(title = "Comparison of True, EM, and Gibbs Distributions with Generated Data",
       x = "x", y = "Probability Density",
       color = "Method", linetype = "Method") +
  theme_minimal() +
  theme(legend.position = "right")

# Step 5: λ Comparison
df_lambda <- data.frame(
  Cluster = factor(rep(1:K, each = 3)),
  Method = rep(c("True", "EM", "Gibbs"), K),
  Lambda = c(data$lambda, em_fit$lambda, gibbs_fit$lambda)
)

p2 <- ggplot(df_lambda, aes(x = Cluster, y = Lambda, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("True" = "blue", "EM" = "orange", "Gibbs" = "forestgreen")) +
  labs(title = "Poisson Rates (λ) Comparison", x = "Cluster", y = "λ Value", fill = "Method") +
  theme_minimal() +
  theme(legend.position = "right")

# Step 6: π Comparison
df_pi <- data.frame(
  Cluster = factor(rep(1:K, each = 3)),
  Method = rep(c("True", "EM", "Gibbs"), K),
  Pi = c(data$pi, em_fit$pi, gibbs_fit$pi)
)

p3 <- ggplot(df_pi, aes(x = Cluster, y = Pi, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  scale_fill_manual(values = c("True" = "blue", "EM" = "orange", "Gibbs" = "forestgreen")) +
  labs(title = "Mixing Coefficients (π) Comparison", x = "Cluster", y = "π Value", fill = "Method") +
  theme_minimal() +
  theme(legend.position = "right")

# Step 7: Plot All
grid.arrange(p1, p2, p3, nrow = 3)
