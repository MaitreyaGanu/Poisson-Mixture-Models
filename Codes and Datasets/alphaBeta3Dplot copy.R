#==================== Libraries =====================#
library(ggplot2)
library(aricode)    # for ARI
library(tidyr)
library(gridExtra)
library(pbapply)    # for progress bars
library(plotly)
library(dplyr)

set.seed(43)

#==================== Parameters =====================#
N          <- 1000                # samples per replicate
K_true     <- 5
true_feats <- 5
noise_seq  <- seq(0, 95, 5)       # noise in steps of 5
alphas     <- c(1, 5, 10)         # tunable
betas      <- c(0.001, 0.01, 0.1) # tunable
n_reps     <- 50                  # MC reps per combo

init_types <- c("true", "random", "top20EM", "top20init")

#==================== Helper functions =====================#
mu_sigma_to_r_theta <- function(mu, sigma2) {
  r     <- mu^2 / (sigma2 - mu)
  theta <- mu / sigma2
  r[r <= 0]         <- 1e-6
  theta[theta <= 0] <- 1e-6
  theta[theta >= 1] <- 1 - 1e-6
  list(r = r, theta = theta)
}

compute_full_loglik <- function(X, omega, mu_kd, sigma2_kd) {
  pts <- mu_sigma_to_r_theta(mu_kd, sigma2_kd)
  ll_mat <- sapply(seq_along(omega), function(k) {
    rowSums(dnbinom(X, size = pts$r[k, ], prob = pts$theta[k, ], log = TRUE)) +
      log(omega[k])
  })
  sum(apply(ll_mat, 1, function(v) { m <- max(v); m + log(sum(exp(v - m))) }))
}

run_nb_em <- function(X, K = K_true, alpha, beta, init = "random", z_init = NULL) {
  N <- nrow(X); D <- ncol(X)
  
  if (init == "random") {
    omega     <- rep(1/K, K)
    mu_kd     <- matrix(rgamma(K * D, alpha, beta), nrow = K)
    sigma2_kd <- mu_kd + 1
    tau       <- NULL
  } else if (init == "true") {
    tau                 <- matrix(0, N, K)
    tau[cbind(1:N, z_init)] <- 1
    omega               <- colMeans(tau)
    mu_kd               <- (t(tau) %*% X) / colSums(tau)
    sigma2_kd           <- (t(tau) %*% (X^2)) / colSums(tau) - mu_kd^2
    mask0               <- (sigma2_kd <= mu_kd) | is.na(sigma2_kd)
    sigma2_kd[mask0]    <- mu_kd[mask0] + 1e-3
  } else if (init == "top20EM") {
    vars    <- apply(X, 2, var)
    top_idx <- order(vars, decreasing = TRUE)[1:min(20, ncol(X))]
    fit20   <- run_nb_em(X[, top_idx], K, alpha, beta, "random")
    return(run_nb_em(X, K, alpha, beta, "true", fit20$z_hat))
  } else if (init == "top20init") {
    vars    <- apply(X, 2, var)
    top_idx <- order(vars, decreasing = TRUE)[1:min(20, ncol(X))]
    fit20   <- run_nb_em(X[, top_idx], K, alpha, beta, "random")
    return(run_nb_em(X, K, alpha, beta, "true", fit20$z_hat))
  }
  
  for (it in 1:100) {
    pts <- mu_sigma_to_r_theta(mu_kd, sigma2_kd)
    logtau <- sapply(1:K, function(k)
      rowSums(dnbinom(X, size = pts$r[k, ], prob = pts$theta[k, ], log = TRUE)) +
        log(omega[k])
    )
    mrow <- apply(logtau, 1, max)
    tau  <- exp(logtau - mrow)
    tau  <- tau / rowSums(tau)
    
    N_k       <- colSums(tau)
    omega     <- N_k / N
    mu_kd     <- (t(tau) %*% X) / N_k
    sigma2_kd <- (t(tau) %*% (X^2)) / N_k - mu_kd^2
    mask      <- (sigma2_kd <= mu_kd) | is.na(sigma2_kd)
    sigma2_kd[mask] <- mu_kd[mask] + 1e-3
  }
  
  z_hat <- apply(tau, 1, which.max)
  list(omega = omega, mu = mu_kd, sigma2 = sigma2_kd, z_hat = z_hat)
}

#==================== Simulation grid =====================#
grid <- expand.grid(
  alpha = alphas,
  beta  = betas,
  noise = noise_seq,
  rep   = 1:n_reps,
  init  = init_types,
  stringsAsFactors = FALSE
)
grid$ARI    <- NA_real_
grid$LogLik <- NA_real_

# Progress bar
pb <- txtProgressBar(min = 0, max = nrow(grid), style = 3)

for (i in seq_len(nrow(grid))) {
  par_alpha <- grid$alpha[i]
  par_beta  <- grid$beta[i]
  n_noise   <- grid$noise[i]
  init_type <- grid$init[i]
  
  D      <- true_feats + n_noise
  z_true <- sample(1:K_true, N, replace = TRUE)
  lambda <- matrix(rgamma(K_true * true_feats, par_alpha, par_beta),
                   nrow = K_true, ncol = true_feats)
  if (n_noise > 0) {
    shared_noise <- rgamma(n_noise, par_alpha, par_beta)
    lambda <- cbind(lambda, matrix(rep(shared_noise, each = K_true), nrow = K_true))
  }
  
  X <- matrix(0, nrow = N, ncol = D)
  for (n in 1:N) {
    mu  <- lambda[z_true[n], ]
    var <- mu + 1
    pts <- mu_sigma_to_r_theta(mu, var)
    X[n, ] <- rnbinom(D, size = pts$r, prob = pts$theta)
  }
  
  fit <- run_nb_em(X, K_true, par_alpha, par_beta, init_type,
                   if (init_type == "true") z_true else NULL)
  
  if (length(fit$z_hat) != length(z_true)) {
    warning(sprintf("Skipping iteration %d due to length mismatch.", i))
    next
  }
  
  grid$ARI[i]    <- ARI(fit$z_hat, z_true)
  grid$LogLik[i] <- compute_full_loglik(X, fit$omega, fit$mu, fit$sigma2)
  setTxtProgressBar(pb, i)
}
close(pb)

#==================== Visualizations =====================#
# 1) ARI vs Noise, faceted by init
p_ari <- ggplot(grid, aes(factor(noise), ARI, fill = init)) +
  geom_boxplot() +
  facet_wrap(~init, ncol = 2) +
  labs(
    title = "ARI vs Noise for Different EM Initializations",
    x = "Number of Noise Features",
    y = "Adjusted Rand Index"
  ) +
  theme_minimal(base_size = 14)

# 2) LogLik vs Noise, faceted by init
p_ll <- ggplot(grid, aes(factor(noise), LogLik, fill = init)) +
  geom_boxplot() +
  facet_wrap(~init, ncol = 2) +
  labs(
    title = "Log‑Likelihood vs Noise for Different EM Initializations",
    x = "Number of Noise Features",
    y = "Full‑Data Log‑Likelihood"
  ) +
  theme_minimal(base_size = 14)

# 3) 3D Surface Plot of Mean ARI at Noise = 0, init = true
mean_ari <- grid %>%
  filter(init == "true", noise == 0) %>%
  group_by(alpha, beta) %>%
  summarize(MeanARI = mean(ARI), .groups = "drop")

p_heat <- plot_ly(
  mean_ari, x = ~alpha, y = ~beta, z = ~MeanARI,
  type = "surface"
) %>%
  layout(
    title = "Mean ARI over (α, β) — True Init, 0 Noise",
    scene = list(
      xaxis = list(title = "alpha"),
      yaxis = list(title = "beta"),
      zaxis = list(title = "Mean ARI")
    )
  )

# Print plots
print(p_ari)
print(p_ll)
p_heat  # interactive plot
