# 0. Libraries and Setup
# Ensure these packages are installed: install.packages(c("flexmix", "mclust", "ggplot2", "dplyr", "tidyr"))
library(flexmix)
library(mclust) # For adjustedRandIndex
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer

set.seed(123) # For reproducibility

# --------------------------------------------------------------------------
# 1. Custom Poisson EM Algorithm Implementation
# --------------------------------------------------------------------------

#' Calculate log sum exp for numerical stability
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

#' Single run of Poisson EM algorithm
#'
#' @param data Matrix of observations (samples x features)
#' @param k Number of clusters
#' @param max_iter Maximum number of iterations
#' @param tol Tolerance for convergence
#' @param epsilon_lambda Small constant to prevent lambdas from being zero
#' @return A list containing lambdas, pi, responsibilities, log_likelihood, and cluster assignments
poisson_em_single_run <- function(data, k, max_iter = 100, tol = 1e-6, epsilon_lambda = 1e-6) {
  n_samples <- nrow(data)
  n_features <- ncol(data)
  
  # Initialization: Randomly assign samples to clusters, then compute initial params
  initial_assignments <- sample(1:k, n_samples, replace = TRUE)
  pi_k <- as.numeric(table(factor(initial_assignments, levels = 1:k))) / n_samples
  pi_k[pi_k == 0] <- epsilon_lambda # Avoid pi_k being zero
  pi_k <- pi_k / sum(pi_k)
  
  lambdas_jk <- matrix(epsilon_lambda, nrow = k, ncol = n_features)
  for (j in 1:k) {
    if (sum(initial_assignments == j) > 0) {
      cluster_data <- data[initial_assignments == j, , drop = FALSE]
      lambdas_jk[j, ] <- pmax(epsilon_lambda, colMeans(cluster_data, na.rm = TRUE))
    } else { # Handle empty initial cluster
      lambdas_jk[j, ] <- pmax(epsilon_lambda, colMeans(data, na.rm = TRUE) + rnorm(n_features, 0, 0.1*mean(data, na.rm=TRUE))) # Small random lambdas relative to data mean
    }
  }
  lambdas_jk[is.na(lambdas_jk)] <- epsilon_lambda # Ensure no NAs if colMeans resulted in NA
  lambdas_jk[lambdas_jk <= 0] <- epsilon_lambda # Ensure positivity
  
  log_likelihood_old <- -Inf
  responsibilities_ij <- matrix(0, nrow = n_samples, ncol = k)
  
  for (iter in 1:max_iter) {
    # E-step: Calculate responsibilities
    log_prob_data_given_cluster <- matrix(0, nrow = n_samples, ncol = k)
    for (j in 1:k) {
      # Sum of log probabilities of each feature for a given cluster
      # dpois returns log probabilities if log = TRUE
      # Need to sum these log probabilities for all features for each sample
      log_prob_data_given_cluster[, j] <- rowSums(
        sapply(1:n_features, function(f) dpois(data[, f], lambdas_jk[j, f], log = TRUE)),
        na.rm = TRUE # sum over features, handling potential -Inf from dpois
      )
    }
    
    weighted_log_prob <- sweep(log_prob_data_given_cluster, 2, log(pi_k), "+")
    
    for (i in 1:n_samples) {
      log_sum_exp_row_i <- log_sum_exp(weighted_log_prob[i, ])
      if(is.finite(log_sum_exp_row_i)){ # Check if log_sum_exp is finite
        responsibilities_ij[i, ] <- exp(weighted_log_prob[i, ] - log_sum_exp_row_i)
      } else { # if log_sum_exp_row_i is -Inf (all elements of weighted_log_prob[i,] were -Inf)
        responsibilities_ij[i, ] <- 1/k # assign equal probability or handle as error
      }
    }
    # Handle cases where a row in responsibilities_ij might sum to 0 due to underflow/all log_prob_data_given_cluster being -Inf
    # or if weighted_log_prob results in all -Inf for a sample
    row_sums_resp <- rowSums(responsibilities_ij, na.rm = TRUE)
    problematic_rows <- (row_sums_resp < epsilon_lambda | is.na(row_sums_resp) | is.nan(row_sums_resp))
    if(any(problematic_rows)){
      responsibilities_ij[problematic_rows, ] <- 1/k # Distribute probability equally
    }
    # Ensure responsibilities sum to 1 for each sample
    responsibilities_ij <- responsibilities_ij / rowSums(responsibilities_ij, na.rm = TRUE) 
    # If any rowSum was 0 and now has NaNs, set to 1/k
    responsibilities_ij[is.nan(responsibilities_ij)] <- 1/k
    
    
    # M-step: Update parameters
    N_k <- colSums(responsibilities_ij, na.rm = TRUE)
    pi_k_new <- N_k / n_samples
    # Ensure pi_k_new sums to 1 and no component is zero, especially if a cluster becomes empty
    pi_k_new[N_k < epsilon_lambda] <- epsilon_lambda 
    pi_k_new <- pi_k_new / sum(pi_k_new)
    
    
    lambdas_jk_new <- matrix(epsilon_lambda, nrow = k, ncol = n_features)
    for (j in 1:k) {
      if (N_k[j] > epsilon_lambda) { # only update if cluster is not effectively empty
        weighted_sum_data <- colSums(responsibilities_ij[, j] * data, na.rm = TRUE)
        lambdas_jk_new[j, ] <- pmax(epsilon_lambda, weighted_sum_data / N_k[j])
      } else { # Reinitialize lambdas for nearly empty cluster
        # Take a small random sample from data to initialize lambda for this empty cluster
        sample_indices_for_reinit <- sample(1:n_samples, min(n_samples, 5), replace = FALSE)
        reinit_data_mean <- colMeans(data[sample_indices_for_reinit, ,drop=FALSE], na.rm = TRUE)
        reinit_data_mean[is.na(reinit_data_mean)] <- mean(data, na.rm=TRUE) # fallback if sample means are NA
        lambdas_jk_new[j, ] <- pmax(epsilon_lambda, reinit_data_mean + rnorm(n_features, 0, 0.1*mean(data, na.rm=TRUE)))
      }
    }
    lambdas_jk_new[is.na(lambdas_jk_new)] <- epsilon_lambda # if colMeans was NA
    lambdas_jk_new[lambdas_jk_new <= 0] <- epsilon_lambda # ensure positivity
    
    
    pi_k <- pi_k_new
    lambdas_jk <- lambdas_jk_new
    
    # Calculate log-likelihood
    current_log_likelihood <- 0
    log_prob_data_given_cluster_updated <- matrix(0, nrow = n_samples, ncol = k)
    for (j in 1:k) {
      log_prob_data_given_cluster_updated[, j] <- rowSums(
        sapply(1:n_features, function(f) dpois(data[, f], lambdas_jk[j, f], log = TRUE)),
        na.rm = TRUE
      )
    }
    for (i in 1:n_samples) {
      # log P(x_i) = log sum_j P(x_i | z_i=j, lambda_j)P(z_i=j | pi_j)
      # log P(x_i) = log sum_j exp( log P(x_i | z_i=j, lambda_j) + log pi_j )
      # This is log_sum_exp of (log_prob_data_given_cluster_updated[i,j] + log(pi_k[j])) over j
      log_likelihood_val_i <- log_sum_exp(log(pi_k) + log_prob_data_given_cluster_updated[i, ])
      if(is.finite(log_likelihood_val_i)){
        current_log_likelihood <- current_log_likelihood + log_likelihood_val_i
      } else {
        # If -Inf, it might pull down the total log-likelihood significantly
        # but it's better than NA. If it's NA, there's a deeper issue.
        current_log_likelihood <- current_log_likelihood - 1e10 # Add a large penalty
      }
    }
    
    if (is.na(current_log_likelihood) || !is.finite(current_log_likelihood)) { # Check for NA or +/- Inf
      current_log_likelihood <- -Inf # Penalize bad runs
      # warning(paste("Iteration", iter, ": Log-likelihood became NA or Inf. Breaking."))
      break 
    }
    if (abs(current_log_likelihood - log_likelihood_old) / (abs(log_likelihood_old) + 1e-3) < tol && iter > 1) { # Relative tolerance, ensure at least one iter
      break
    }
    log_likelihood_old <- current_log_likelihood
  }
  
  cluster_assignments <- apply(responsibilities_ij, 1, which.max)
  
  if(any(tabulate(cluster_assignments, nbins=k) == 0) && k > 1){ # Check if any cluster is empty
    # message(paste("Warning: k=",k, ", one or more clusters became empty. Penalizing likelihood.", sep=""))
    current_log_likelihood <- -Inf # Penalize if clusters were lost (unless k=1)
  }
  
  
  return(list(
    lambdas = lambdas_jk,
    pi = pi_k,
    responsibilities = responsibilities_ij,
    log_likelihood = current_log_likelihood,
    cluster_assignments = cluster_assignments,
    iterations = iter
  ))
}

#' Poisson EM with multiple initializations
#'
#' @param data Matrix of observations
#' @param k Number of clusters
#' @param n_inits Number of random initializations
#' @param ... Other parameters for poisson_em_single_run
#' @return The best result from multiple runs based on log-likelihood
poisson_em_multistart <- function(data, k, n_inits = 5, ...) {
  best_log_likelihood <- -Inf
  best_result <- NULL
  
  if (k <= 0) stop("Number of clusters k must be at least 1.")
  if (k == 1) { 
    n_samples <- nrow(data)
    n_features <- ncol(data)
    lambdas_jk <- matrix(pmax(1e-6, colMeans(data, na.rm = TRUE)), nrow = 1, ncol = n_features)
    lambdas_jk[is.na(lambdas_jk)] <- 1e-6 # If all data for a feature is NA
    pi_k <- 1
    responsibilities_ij <- matrix(1, nrow = n_samples, ncol = 1)
    
    log_likelihood <- sum(sapply(1:n_samples, function(sample_idx) {
      sum(dpois(data[sample_idx,], lambdas_jk[1,], log = TRUE), na.rm=TRUE) 
    }), na.rm=TRUE)
    if(is.na(log_likelihood) || !is.finite(log_likelihood)) log_likelihood <- -Inf 
    
    return(list(
      lambdas = lambdas_jk,
      pi = pi_k,
      responsibilities = responsibilities_ij,
      log_likelihood = log_likelihood,
      cluster_assignments = rep(1, n_samples),
      iterations = 1 
    ))
  }
  
  
  for (i_init in 1:n_inits) {
    current_result <- tryCatch({
      poisson_em_single_run(data, k, ...)
    }, error = function(e) {
      # message(sprintf("Error in EM run %d for k=%d: %s. Log-likelihood set to -Inf.", i_init, k, e$message))
      list(log_likelihood = -Inf, cluster_assignments = rep(1, nrow(data))) # Provide default assignments
    })
    
    if (!is.null(current_result$log_likelihood) && 
        !is.na(current_result$log_likelihood) &&
        is.finite(current_result$log_likelihood) && # Ensure it's not Inf
        current_result$log_likelihood > best_log_likelihood) {
      best_log_likelihood <- current_result$log_likelihood
      best_result <- current_result
    }
  }
  if (is.null(best_result)) { 
    # warning(paste("All EM initializations failed to produce a valid model for k =", k, ". Returning a dummy structure."))
    return(list( 
      lambdas = matrix(1e-6, nrow = k, ncol = ncol(data)),
      pi = rep(1/k, k),
      responsibilities = matrix(1/k, nrow = nrow(data), ncol = k),
      log_likelihood = -Inf,
      cluster_assignments = rep(1, nrow(data)), 
      iterations = 0
    ))
  }
  return(best_result)
}

# --------------------------------------------------------------------------
# 2. Data Generation Function
# --------------------------------------------------------------------------
#' Generate Poisson mixture data
#'
#' @param n_samples_per_cluster Number of samples per cluster
#' @param n_features Total number of features
#' @param num_noisy_features Number of features that are noisy (same lambda across clusters)
#' @param lambdas_dist_pattern Vector of lambdas for distinguishing features (one lambda per cluster)
#' @param lambda_noise_value Lambda value for noisy features
#' @param true_k True number of clusters (length of lambdas_dist_pattern)
#' @return A list containing the data matrix and true labels
generate_poisson_mixture_data <- function(n_samples_per_cluster,
                                          n_features,
                                          num_noisy_features,
                                          lambdas_dist_pattern, 
                                          lambda_noise_value,
                                          true_k) {
  
  if (length(lambdas_dist_pattern) != true_k) {
    stop("Length of lambdas_dist_pattern must be equal to true_k.")
  }
  if (num_noisy_features < 0 || num_noisy_features > n_features) {
    stop("num_noisy_features must be between 0 and n_features.")
  }
  
  
  n_total_samples <- n_samples_per_cluster * true_k
  true_labels <- rep(1:true_k, each = n_samples_per_cluster)
  
  data_matrix <- matrix(0, nrow = n_total_samples, ncol = n_features)
  
  num_dist_features <- n_features - num_noisy_features
  
  all_feature_indices <- 1:n_features
  if (num_noisy_features > 0) {
    # Ensure sample doesn't fail if n_features = 0 (edge case, not here)
    noisy_feature_indices <- sample(all_feature_indices, num_noisy_features)
  } else {
    noisy_feature_indices <- c()
  }
  dist_feature_indices <- setdiff(all_feature_indices, noisy_feature_indices)
  
  
  for (i in 1:n_total_samples) {
    current_cluster <- true_labels[i]
    
    if (num_dist_features > 0 && length(dist_feature_indices) > 0) { 
      lambda_for_dist_feature <- lambdas_dist_pattern[current_cluster]
      data_matrix[i, dist_feature_indices] <- rpois(length(dist_feature_indices), lambda_for_dist_feature)
    }
    
    if (num_noisy_features > 0 && length(noisy_feature_indices) > 0) { 
      data_matrix[i, noisy_feature_indices] <- rpois(length(noisy_feature_indices), lambda_noise_value)
    }
  }
  # Ensure no NAs in data_matrix, replace with 0 if any rpois produced NA (should not happen for lambda > 0)
  data_matrix[is.na(data_matrix)] <- 0
  return(list(data = data_matrix, true_labels = true_labels))
}


# --------------------------------------------------------------------------
# 3. Simulation Parameters
# --------------------------------------------------------------------------
N_FEATURES <- 100 
N_SAMPLES_PER_CLUSTER <- 100 
TRUE_K <- 3
LAMBDAS_DIST_PATTERN <- c(2, 7, 12) 
LAMBDA_NOISE_VALUE <- mean(LAMBDAS_DIST_PATTERN) 
N_EM_INITS <- 5 
FLEXMIX_NREP <- 5 
MAX_EM_ITER <- 100
EM_TOL <- 1e-5 
K_RANGE_BIC <- 2:5 

results_list <- list()


# --------------------------------------------------------------------------
# 4. Simulation Loop
# --------------------------------------------------------------------------
cat("Starting simulation...\n")
for (num_noisy in 0:N_FEATURES) { 
  cat(sprintf("Processing: Dataset with %d noisy features out of %d\n", num_noisy, N_FEATURES))
  
  sim_data <- generate_poisson_mixture_data(
    n_samples_per_cluster = N_SAMPLES_PER_CLUSTER,
    n_features = N_FEATURES,
    num_noisy_features = num_noisy,
    lambdas_dist_pattern = LAMBDAS_DIST_PATTERN,
    lambda_noise_value = LAMBDA_NOISE_VALUE,
    true_k = TRUE_K
  )
  data_matrix <- sim_data$data
  true_labels <- sim_data$true_labels
  
  # --- Part 1: Fixed k = TRUE_K ---
  time_myem_fixed_start <- Sys.time()
  myem_fixed_res <- poisson_em_multistart(data_matrix, k = TRUE_K, n_inits = N_EM_INITS, 
                                          max_iter = MAX_EM_ITER, tol = EM_TOL, epsilon_lambda = 1e-6)
  time_myem_fixed_end <- Sys.time()
  time_myem_fixed <- as.numeric(difftime(time_myem_fixed_end, time_myem_fixed_start, units = "secs"))
  ari_myem_fixed <- if (!is.null(myem_fixed_res) && !is.null(myem_fixed_res$cluster_assignments) && is.finite(myem_fixed_res$log_likelihood) && myem_fixed_res$log_likelihood > -Inf) {
    tryCatch(adjustedRandIndex(true_labels, myem_fixed_res$cluster_assignments), error = function(e) NA_real_)
  } else { NA_real_ }
  
  time_flexmix_fixed_start <- Sys.time()
  flexmix_fixed_model <- tryCatch({
    flexmix(data_matrix ~ 1, k = TRUE_K, model = FLXMCmvpois(), 
            control = list(nrep = FLEXMIX_NREP, iter.max = MAX_EM_ITER, tol = EM_TOL, verbose=0))
  }, error = function(e) { NULL })
  time_flexmix_fixed_end <- Sys.time()
  time_flexmix_fixed <- as.numeric(difftime(time_flexmix_fixed_end, time_flexmix_fixed_start, units = "secs"))
  ari_flexmix_fixed <- if (!is.null(flexmix_fixed_model)) {
    tryCatch(adjustedRandIndex(true_labels, clusters(flexmix_fixed_model)), error = function(e) NA_real_)
  } else { NA_real_ }
  
  # --- Part 2: Variable k, selected by BIC ---
  time_myem_bic_start <- Sys.time()
  best_bic_myem <- Inf
  best_k_myem <- NA_integer_
  final_myem_model_bic <- NULL
  
  for (k_test in K_RANGE_BIC) {
    if (k_test > nrow(data_matrix)) next 
    current_myem_res_bic <- poisson_em_multistart(data_matrix, k = k_test, n_inits = N_EM_INITS, 
                                                  max_iter = MAX_EM_ITER, tol = EM_TOL, epsilon_lambda = 1e-6)
    if (!is.null(current_myem_res_bic) && !is.na(current_myem_res_bic$log_likelihood) && is.finite(current_myem_res_bic$log_likelihood) && current_myem_res_bic$log_likelihood > -Inf ) {
      n_params <- (k_test - 1) + (k_test * N_FEATURES) 
      bic_val <- -2 * current_myem_res_bic$log_likelihood + n_params * log(nrow(data_matrix))
      if (!is.na(bic_val) && is.finite(bic_val) && bic_val < best_bic_myem) {
        best_bic_myem <- bic_val
        best_k_myem <- k_test
        final_myem_model_bic <- current_myem_res_bic
      }
    }
  }
  time_myem_bic_end <- Sys.time()
  time_myem_bic <- as.numeric(difftime(time_myem_bic_end, time_myem_bic_start, units = "secs"))
  k_selected_myem <- best_k_myem
  ari_myem_bic <- if (!is.null(final_myem_model_bic) && !is.null(final_myem_model_bic$cluster_assignments)) {
    tryCatch(adjustedRandIndex(true_labels, final_myem_model_bic$cluster_assignments), error = function(e) NA_real_)
  } else { NA_real_ }
  
  time_flexmix_bic_start <- Sys.time()
  flexmix_bic_model_step <- tryCatch({
    stepFlexmix(data_matrix ~ 1, k = K_RANGE_BIC, model = FLXMCmvpois(), 
                nrep = FLEXMIX_NREP, verbose = FALSE,
                control = list(iter.max = MAX_EM_ITER, tol = EM_TOL, verbose=0)) 
  }, error = function(e) { NULL })
  time_flexmix_bic_end <- Sys.time()
  time_flexmix_bic <- as.numeric(difftime(time_flexmix_bic_end, time_flexmix_bic_start, units = "secs"))
  
  k_selected_flexmix <- NA_integer_
  ari_flexmix_bic <- NA_real_
  best_flexmix_model_bic <- NULL
  
  if (!is.null(flexmix_bic_model_step)) {
    best_flexmix_model_bic <- tryCatch(getModel(flexmix_bic_model_step, "BIC"), error = function(e) NULL)
    
    if (!is.null(best_flexmix_model_bic)) {
      k_selected_flexmix <- best_flexmix_model_bic@k 
      cl_flex_bic <- tryCatch(clusters(best_flexmix_model_bic), error = function(e) NULL)
      if (!is.null(cl_flex_bic)) {
        ari_flexmix_bic <- tryCatch(adjustedRandIndex(true_labels, cl_flex_bic), error = function(e) NA_real_)
      }
    }
  }
  
  results_list[[length(results_list) + 1]] <- data.frame(
    NumNoisyFeatures = num_noisy,
    ARI_MyEM_FixedK = ari_myem_fixed, Time_MyEM_FixedK = time_myem_fixed,
    ARI_Flexmix_FixedK = ari_flexmix_fixed, Time_Flexmix_FixedK = time_flexmix_fixed,
    ARI_MyEM_BIC = ari_myem_bic, Time_MyEM_BIC = time_myem_bic, K_MyEM_BIC = k_selected_myem,
    ARI_Flexmix_BIC = ari_flexmix_bic, Time_Flexmix_BIC = time_flexmix_bic, K_Flexmix_BIC = k_selected_flexmix
  )
}
cat("Simulation finished.\n")

results_df <- bind_rows(results_list)

# Ensure numeric types for ARI columns just in case
results_df <- results_df %>%
  mutate(across(starts_with("ARI_"), as.numeric))

cat("\nSummary of results_df before plotting:\n")
print(summary(results_df)) # ADDED SUMMARY
cat("\n")


# --------------------------------------------------------------------------
# 5. Plotting Results
# --------------------------------------------------------------------------
cat("Generating plots...\n")

# Plot ARI vs. NumNoisyFeatures (Fixed K)
plot_ari_fixed_k <- ggplot(results_df, aes(x = NumNoisyFeatures)) +
  geom_line(aes(y = ARI_MyEM_FixedK, color = "My EM (Fixed K)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_MyEM_FixedK, color = "My EM (Fixed K)"), na.rm = TRUE) +
  geom_line(aes(y = ARI_Flexmix_FixedK, color = "Flexmix (Fixed K)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_Flexmix_FixedK, color = "Flexmix (Fixed K)"), na.rm = TRUE) +
  labs(title = paste("ARI vs. Noisy Features (True K =", TRUE_K, ")"),
       x = "Number of Noisy Features (i)", y = "Adjusted Rand Index (ARI)",
       color = "Algorithm") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (Fixed K)" = "blue", "Flexmix (Fixed K)" = "red"))
# ylim(0,1) REMOVED

print(plot_ari_fixed_k)

# Plot ARI vs. NumNoisyFeatures (BIC Selected K)
plot_ari_bic <- ggplot(results_df, aes(x = NumNoisyFeatures)) +
  geom_line(aes(y = ARI_MyEM_BIC, color = "My EM (BIC)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_MyEM_BIC, color = "My EM (BIC)"), na.rm = TRUE) +
  geom_line(aes(y = ARI_Flexmix_BIC, color = "Flexmix (BIC)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_Flexmix_BIC, color = "Flexmix (BIC)"), na.rm = TRUE) +
  labs(title = "ARI vs. Noisy Features (BIC Selected K)",
       x = "Number of Noisy Features (i)", y = "Adjusted Rand Index (ARI)",
       color = "Algorithm") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange"))
# ylim(0,1) REMOVED

print(plot_ari_bic)

# Plot Selected K by BIC vs. NumNoisyFeatures
# Create a placeholder for NA K values for plotting if you want to see them distinctly
results_df_plot_k <- results_df %>%
  mutate(
    K_MyEM_BIC_plot = ifelse(is.na(K_MyEM_BIC), min(K_RANGE_BIC, na.rm=TRUE)-0.5, K_MyEM_BIC), # Plot NAs below lowest K
    K_Flexmix_BIC_plot = ifelse(is.na(K_Flexmix_BIC), min(K_RANGE_BIC, na.rm=TRUE)-0.5, K_Flexmix_BIC)
  )

results_df_k_long <- results_df_plot_k %>%
  select(NumNoisyFeatures, K_MyEM_BIC_plot, K_Flexmix_BIC_plot) %>%
  rename(K_MyEM_BIC = K_MyEM_BIC_plot, K_Flexmix_BIC = K_Flexmix_BIC_plot) %>%
  pivot_longer(cols = c(K_MyEM_BIC, K_Flexmix_BIC), names_to = "Algorithm", values_to = "SelectedK") %>%
  mutate(Algorithm = case_when(
    Algorithm == "K_MyEM_BIC" ~ "My EM (BIC)",
    Algorithm == "K_Flexmix_BIC" ~ "Flexmix (BIC)",
    TRUE ~ Algorithm
  ))

plot_selected_k <- ggplot(results_df_k_long, aes(x = NumNoisyFeatures, y = SelectedK, color = Algorithm)) +
  geom_line(na.rm = TRUE) + 
  geom_point(na.rm = TRUE) +
  geom_hline(yintercept = TRUE_K, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(results_df$NumNoisyFeatures, na.rm=TRUE)*0.9, y = TRUE_K + 0.2, label = paste("True K =", TRUE_K), color="gray30") +
  labs(title = "Selected Number of Clusters (K) by BIC",
       x = "Number of Noisy Features (i)", y = "Selected K",
       color = "Algorithm") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange")) +
  scale_y_continuous(breaks = (min(K_RANGE_BIC,na.rm=TRUE)-1):max(K_RANGE_BIC,na.rm=TRUE), 
                     labels = function(x) ifelse(x < min(K_RANGE_BIC, na.rm=TRUE), "NA/Fail", x))


print(plot_selected_k)


# Plot Time vs. NumNoisyFeatures
results_df_time_long <- results_df %>%
  select(NumNoisyFeatures, Time_MyEM_FixedK, Time_Flexmix_FixedK, Time_MyEM_BIC, Time_Flexmix_BIC) %>%
  pivot_longer(cols = -NumNoisyFeatures, names_to = "Algorithm_Metric", values_to = "Time") %>%
  mutate(
    Algorithm = case_when(
      grepl("MyEM_FixedK", Algorithm_Metric) ~ "My EM (Fixed K)",
      grepl("Flexmix_FixedK", Algorithm_Metric) ~ "Flexmix (Fixed K)",
      grepl("MyEM_BIC", Algorithm_Metric) ~ "My EM (BIC)",
      grepl("Flexmix_BIC", Algorithm_Metric) ~ "Flexmix (BIC)",
      TRUE ~ Algorithm_Metric
    )
  )

plot_time <- ggplot(results_df_time_long, aes(x = NumNoisyFeatures, y = Time, color = Algorithm)) +
  geom_line(na.rm = TRUE) +
  geom_point(na.rm = TRUE) +
  labs(title = "Execution Time vs. Number of Noisy Features",
       x = "Number of Noisy Features (i)", y = "Time (seconds)",
       color = "Algorithm & Mode") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (Fixed K)" = "blue", "Flexmix (Fixed K)" = "red",
                                "My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange")) #+
# scale_y_log10() # Consider using log scale if time varies greatly; can cause issues if times are zero or NA

print(plot_time)


# Save results if needed
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_filename <- paste0("poisson_em_flexmix_comparison_results_", timestamp, ".csv")
# write.csv(results_df, results_filename, row.names = FALSE)
# cat(paste("Results saved to:", results_filename, "\n"))

cat("Script finished. Check plots.\n")
#------------------------------------------------------------------------------------------------------------------------------------------------
# 0. Libraries and Setup
# Ensure these packages are installed: install.packages(c("flexmix", "mclust", "ggplot2", "dplyr", "tidyr"))
library(flexmix)
library(mclust) # For adjustedRandIndex
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer

set.seed(123) # For reproducibility




















# --------------------------------------------------------------------------
# 1. Custom Poisson EM Algorithm Implementation with Gamma Prior for Lambdas
# --------------------------------------------------------------------------

#' Calculate log sum exp for numerical stability
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

#' Single run of Poisson EM algorithm with Gamma prior on lambdas
poisson_em_single_run <- function(data, k, gamma_prior_shape, gamma_prior_rate, 
                                  max_iter = 100, tol = 1e-6, epsilon_lambda = 1e-6) {
  n_samples <- nrow(data)
  n_features <- ncol(data)
  
  initial_assignments <- sample(1:k, n_samples, replace = TRUE)
  pi_k <- as.numeric(table(factor(initial_assignments, levels = 1:k))) / n_samples
  pi_k[pi_k == 0] <- epsilon_lambda 
  pi_k <- pi_k / sum(pi_k)
  
  lambdas_jk <- matrix(epsilon_lambda, nrow = k, ncol = n_features)
  for (j in 1:k) {
    if (sum(initial_assignments == j) > 0) {
      cluster_data <- data[initial_assignments == j, , drop = FALSE]
      lambdas_jk[j, ] <- pmax(epsilon_lambda, colMeans(cluster_data, na.rm = TRUE))
    } else { 
      prior_mean_lambda = gamma_prior_shape / gamma_prior_rate
      lambdas_jk[j, ] <- pmax(epsilon_lambda, rep(prior_mean_lambda, n_features) + rnorm(n_features, 0, 0.1*prior_mean_lambda))
    }
  }
  lambdas_jk[is.na(lambdas_jk)] <- epsilon_lambda 
  lambdas_jk[lambdas_jk <= 0] <- epsilon_lambda 
  
  log_likelihood_old <- -Inf
  responsibilities_ij <- matrix(0, nrow = n_samples, ncol = k)
  
  for (iter in 1:max_iter) {
    log_prob_data_given_cluster <- matrix(0, nrow = n_samples, ncol = k)
    for (j in 1:k) {
      log_prob_data_given_cluster[, j] <- rowSums(
        sapply(1:n_features, function(f) dpois(data[, f], lambdas_jk[j, f], log = TRUE)),
        na.rm = TRUE 
      )
    }
    
    weighted_log_prob <- sweep(log_prob_data_given_cluster, 2, log(pi_k), "+")
    
    for (i in 1:n_samples) {
      log_sum_exp_row_i <- log_sum_exp(weighted_log_prob[i, ])
      if(is.finite(log_sum_exp_row_i)){ 
        responsibilities_ij[i, ] <- exp(weighted_log_prob[i, ] - log_sum_exp_row_i)
      } else { 
        responsibilities_ij[i, ] <- 1/k 
      }
    }
    row_sums_resp <- rowSums(responsibilities_ij, na.rm = TRUE)
    problematic_rows <- (row_sums_resp < epsilon_lambda | is.na(row_sums_resp) | is.nan(row_sums_resp))
    if(any(problematic_rows)){
      responsibilities_ij[problematic_rows, ] <- 1/k 
    }
    responsibilities_ij <- responsibilities_ij / rowSums(responsibilities_ij, na.rm = TRUE) 
    responsibilities_ij[is.nan(responsibilities_ij)] <- 1/k
    
    N_k <- colSums(responsibilities_ij, na.rm = TRUE) 
    pi_k_new <- N_k / n_samples
    pi_k_new[N_k < epsilon_lambda] <- epsilon_lambda 
    pi_k_new <- pi_k_new / sum(pi_k_new)
    
    lambdas_jk_new <- matrix(epsilon_lambda, nrow = k, ncol = n_features)
    for (j in 1:k) {
      if (N_k[j] > epsilon_lambda) { 
        weighted_sum_data_f <- colSums(responsibilities_ij[, j] * data, na.rm = TRUE) 
        lambdas_jk_new[j, ] <- (gamma_prior_shape + weighted_sum_data_f) / (gamma_prior_rate + N_k[j])
        lambdas_jk_new[j, ] <- pmax(epsilon_lambda, lambdas_jk_new[j, ]) 
      } else { 
        prior_mean_lambda = gamma_prior_shape / gamma_prior_rate
        lambdas_jk_new[j, ] <- pmax(epsilon_lambda, rep(prior_mean_lambda, n_features) + rnorm(n_features, 0, 0.01*prior_mean_lambda))
      }
    }
    lambdas_jk_new[is.na(lambdas_jk_new)] <- epsilon_lambda 
    lambdas_jk_new[lambdas_jk_new <= 0] <- epsilon_lambda 
    
    pi_k <- pi_k_new
    lambdas_jk <- lambdas_jk_new
    
    current_log_likelihood <- 0
    log_prob_data_given_cluster_updated <- matrix(0, nrow = n_samples, ncol = k)
    for (j in 1:k) {
      log_prob_data_given_cluster_updated[, j] <- rowSums(
        sapply(1:n_features, function(f) dpois(data[, f], lambdas_jk[j, f], log = TRUE)),
        na.rm = TRUE
      )
    }
    for (i in 1:n_samples) {
      log_likelihood_val_i <- log_sum_exp(log(pi_k) + log_prob_data_given_cluster_updated[i, ])
      if(is.finite(log_likelihood_val_i)){
        current_log_likelihood <- current_log_likelihood + log_likelihood_val_i
      } else {
        current_log_likelihood <- current_log_likelihood - 1e10 
      }
    }
    
    if (is.na(current_log_likelihood) || !is.finite(current_log_likelihood)) { 
      current_log_likelihood <- -Inf 
      break 
    }
    if (abs(current_log_likelihood - log_likelihood_old) / (abs(log_likelihood_old) + 1e-3) < tol && iter > 1) { 
      break
    }
    log_likelihood_old <- current_log_likelihood
  }
  
  cluster_assignments <- apply(responsibilities_ij, 1, which.max)
  if(any(tabulate(cluster_assignments, nbins=k) == 0) && k > 1){ 
    current_log_likelihood <- -Inf 
  }
  
  return(list(
    lambdas = lambdas_jk, pi = pi_k, responsibilities = responsibilities_ij,
    log_likelihood = current_log_likelihood, cluster_assignments = cluster_assignments, iterations = iter
  ))
}

#' Poisson EM with multiple initializations and Gamma prior
poisson_em_multistart <- function(data, k, gamma_prior_shape, gamma_prior_rate, n_inits = 5, ...) {
  best_log_likelihood <- -Inf
  best_result <- NULL
  
  if (k <= 0) stop("Number of clusters k must be at least 1.")
  if (k == 1) { 
    n_samples <- nrow(data)
    n_features <- ncol(data)
    sum_data_f <- colSums(data, na.rm = TRUE) 
    lambdas_1k_f <- (gamma_prior_shape + sum_data_f) / (gamma_prior_rate + n_samples)
    lambdas_jk <- matrix(pmax(1e-6, lambdas_1k_f), nrow = 1, ncol = n_features)
    lambdas_jk[is.na(lambdas_jk)] <- 1e-6 
    pi_k <- 1
    responsibilities_ij <- matrix(1, nrow = n_samples, ncol = 1)
    log_likelihood <- sum(sapply(1:n_samples, function(sample_idx) {
      sum(dpois(data[sample_idx,], lambdas_jk[1,], log = TRUE), na.rm=TRUE) 
    }), na.rm=TRUE)
    if(is.na(log_likelihood) || !is.finite(log_likelihood)) log_likelihood <- -Inf 
    
    return(list(
      lambdas = lambdas_jk, pi = pi_k, responsibilities = responsibilities_ij,
      log_likelihood = log_likelihood, cluster_assignments = rep(1, n_samples), iterations = 1 
    ))
  }
  
  for (i_init in 1:n_inits) {
    current_result <- tryCatch({
      poisson_em_single_run(data, k, gamma_prior_shape, gamma_prior_rate, ...)
    }, error = function(e) {
      list(log_likelihood = -Inf, cluster_assignments = rep(1, nrow(data))) 
    })
    
    if (!is.null(current_result$log_likelihood) && 
        !is.na(current_result$log_likelihood) &&
        is.finite(current_result$log_likelihood) && 
        current_result$log_likelihood > best_log_likelihood) {
      best_log_likelihood <- current_result$log_likelihood
      best_result <- current_result
    }
  }
  if (is.null(best_result)) { 
    return(list( 
      lambdas = matrix(gamma_prior_shape/gamma_prior_rate, nrow = k, ncol = ncol(data)),
      pi = rep(1/k, k), responsibilities = matrix(1/k, nrow = nrow(data), ncol = k),
      log_likelihood = -Inf, cluster_assignments = rep(1, nrow(data)), iterations = 0
    ))
  }
  return(best_result)
}

# --------------------------------------------------------------------------
# 2. Data Generation Function (MODIFIED to sample true lambdas from Gamma)
# --------------------------------------------------------------------------
generate_poisson_mixture_data <- function(n_samples_per_cluster,
                                          n_features,
                                          num_noisy_features,
                                          true_k,
                                          gamma_shape_for_data_gen, # Hyperparameter for sampling true lambdas
                                          gamma_rate_for_data_gen,  # Hyperparameter for sampling true lambdas
                                          epsilon_lambda = 1e-6) { # Min value for sampled lambdas
  
  if (num_noisy_features < 0 || num_noisy_features > n_features) {
    stop("num_noisy_features must be between 0 and n_features.")
  }
  
  # Sample true lambdas for distinguishing features and for noise from Gamma distribution
  # These are the "true" parameters for THIS specific dataset being generated.
  lambdas_dist_pattern_sampled <- pmax(epsilon_lambda, 
                                       rgamma(true_k, shape = gamma_shape_for_data_gen, rate = gamma_rate_for_data_gen))
  lambda_noise_value_sampled <- pmax(epsilon_lambda, 
                                     rgamma(1, shape = gamma_shape_for_data_gen, rate = gamma_rate_for_data_gen))
  
  # Ensure the distinguishing lambdas are somewhat different, e.g., by sorting or adding small unique jitter if needed for identifiability,
  # though rgamma itself should provide variation. For very specific separation, one might enforce minimum differences.
  # For now, direct samples are used.
  
  n_total_samples <- n_samples_per_cluster * true_k
  true_labels <- rep(1:true_k, each = n_samples_per_cluster)
  data_matrix <- matrix(0, nrow = n_total_samples, ncol = n_features)
  num_dist_features <- n_features - num_noisy_features
  all_feature_indices <- 1:n_features
  if (num_noisy_features > 0) {
    noisy_feature_indices <- sample(all_feature_indices, num_noisy_features)
  } else {
    noisy_feature_indices <- c()
  }
  dist_feature_indices <- setdiff(all_feature_indices, noisy_feature_indices)
  
  for (i in 1:n_total_samples) {
    current_cluster <- true_labels[i]
    if (num_dist_features > 0 && length(dist_feature_indices) > 0) { 
      # Use the sampled lambda for the current_cluster's distinguishing features
      lambda_for_dist_feature <- lambdas_dist_pattern_sampled[current_cluster] 
      data_matrix[i, dist_feature_indices] <- rpois(length(dist_feature_indices), lambda_for_dist_feature)
    }
    if (num_noisy_features > 0 && length(noisy_feature_indices) > 0) { 
      # Use the sampled lambda for noisy features
      data_matrix[i, noisy_feature_indices] <- rpois(length(noisy_feature_indices), lambda_noise_value_sampled)
    }
  }
  data_matrix[is.na(data_matrix)] <- 0
  return(list(data = data_matrix, true_labels = true_labels, 
              true_dist_lambdas = lambdas_dist_pattern_sampled, 
              true_noise_lambda = lambda_noise_value_sampled)) # Also return true lambdas for potential inspection
}


# --------------------------------------------------------------------------
# 3. Simulation Parameters
# --------------------------------------------------------------------------
N_FEATURES <- 100 
N_SAMPLES_PER_CLUSTER <- 100 
TRUE_K <- 3
# LAMBDAS_DIST_PATTERN and LAMBDA_NOISE_VALUE are now generated inside the loop

GAMMA_PRIOR_SHAPE <- 3.0 # alpha_0 (Used for EM prior AND for data generation)
GAMMA_PRIOR_RATE <- 7.0  # beta_0  (Used for EM prior AND for data generation)
# With SHAPE=1, RATE=1, prior mean for lambdas is 1. 
# If you want generated lambdas to be larger on average, increase SHAPE or decrease RATE.
# For example, SHAPE=5, RATE=1 would mean data-generating lambdas average around 5.

N_EM_INITS <- 5 
FLEXMIX_NREP <- 5 
MAX_EM_ITER <- 100
EM_TOL <- 1e-5 
K_RANGE_BIC <- 2:5 
EPSILON_LAMBDA_GLOBAL <- 1e-6 # Global epsilon for flooring lambdas

results_list <- list()


# --------------------------------------------------------------------------
# 4. Simulation Loop
# --------------------------------------------------------------------------
cat("Starting simulation...\n")
for (num_noisy in 0:N_FEATURES) { 
  cat(sprintf("Processing: Dataset with %d noisy features out of %d\n", num_noisy, N_FEATURES))
  
  # Generate data, true lambdas are now sampled internally using GAMMA_PRIOR_SHAPE/RATE
  sim_data <- generate_poisson_mixture_data(
    n_samples_per_cluster = N_SAMPLES_PER_CLUSTER,
    n_features = N_FEATURES,
    num_noisy_features = num_noisy,
    true_k = TRUE_K,
    gamma_shape_for_data_gen = GAMMA_PRIOR_SHAPE, # Use global prior params for data gen
    gamma_rate_for_data_gen = GAMMA_PRIOR_RATE,
    epsilon_lambda = EPSILON_LAMBDA_GLOBAL
  )
  data_matrix <- sim_data$data
  true_labels <- sim_data$true_labels
  
  # My EM - Pass prior parameters
  time_myem_fixed_start <- Sys.time()
  myem_fixed_res <- poisson_em_multistart(data_matrix, k = TRUE_K, 
                                          gamma_prior_shape = GAMMA_PRIOR_SHAPE, 
                                          gamma_prior_rate = GAMMA_PRIOR_RATE,
                                          n_inits = N_EM_INITS, 
                                          max_iter = MAX_EM_ITER, tol = EM_TOL, epsilon_lambda = EPSILON_LAMBDA_GLOBAL)
  time_myem_fixed_end <- Sys.time()
  time_myem_fixed <- as.numeric(difftime(time_myem_fixed_end, time_myem_fixed_start, units = "secs"))
  ari_myem_fixed <- if (!is.null(myem_fixed_res) && !is.null(myem_fixed_res$cluster_assignments) && is.finite(myem_fixed_res$log_likelihood) && myem_fixed_res$log_likelihood > -Inf) {
    tryCatch(adjustedRandIndex(true_labels, myem_fixed_res$cluster_assignments), error = function(e) NA_real_)
  } else { NA_real_ }
  
  # Flexmix part remains unchanged
  time_flexmix_fixed_start <- Sys.time()
  flexmix_fixed_model <- tryCatch({
    flexmix(data_matrix ~ 1, k = TRUE_K, model = FLXMCmvpois(), 
            control = list(nrep = FLEXMIX_NREP, iter.max = MAX_EM_ITER, tol = EM_TOL, verbose=0))
  }, error = function(e) { NULL })
  time_flexmix_fixed_end <- Sys.time()
  time_flexmix_fixed <- as.numeric(difftime(time_flexmix_fixed_end, time_flexmix_fixed_start, units = "secs"))
  ari_flexmix_fixed <- if (!is.null(flexmix_fixed_model)) {
    tryCatch(adjustedRandIndex(true_labels, clusters(flexmix_fixed_model)), error = function(e) NA_real_)
  } else { NA_real_ }
  
  # My EM BIC - Pass prior parameters
  time_myem_bic_start <- Sys.time()
  best_bic_myem <- Inf
  best_k_myem <- NA_integer_
  final_myem_model_bic <- NULL
  
  for (k_test in K_RANGE_BIC) {
    if (k_test > nrow(data_matrix)) next 
    current_myem_res_bic <- poisson_em_multistart(data_matrix, k = k_test, 
                                                  gamma_prior_shape = GAMMA_PRIOR_SHAPE,
                                                  gamma_prior_rate = GAMMA_PRIOR_RATE,
                                                  n_inits = N_EM_INITS, 
                                                  max_iter = MAX_EM_ITER, tol = EM_TOL, epsilon_lambda = EPSILON_LAMBDA_GLOBAL)
    if (!is.null(current_myem_res_bic) && !is.na(current_myem_res_bic$log_likelihood) && is.finite(current_myem_res_bic$log_likelihood) && current_myem_res_bic$log_likelihood > -Inf ) {
      n_params <- (k_test - 1) + (k_test * N_FEATURES) 
      bic_val <- -2 * current_myem_res_bic$log_likelihood + n_params * log(nrow(data_matrix))
      if (!is.na(bic_val) && is.finite(bic_val) && bic_val < best_bic_myem) {
        best_bic_myem <- bic_val
        best_k_myem <- k_test
        final_myem_model_bic <- current_myem_res_bic
      }
    }
  }
  time_myem_bic_end <- Sys.time()
  time_myem_bic <- as.numeric(difftime(time_myem_bic_end, time_myem_bic_start, units = "secs"))
  k_selected_myem <- best_k_myem
  ari_myem_bic <- if (!is.null(final_myem_model_bic) && !is.null(final_myem_model_bic$cluster_assignments)) {
    tryCatch(adjustedRandIndex(true_labels, final_myem_model_bic$cluster_assignments), error = function(e) NA_real_)
  } else { NA_real_ }
  
  # Flexmix BIC part remains unchanged
  time_flexmix_bic_start <- Sys.time()
  flexmix_bic_model_step <- tryCatch({
    stepFlexmix(data_matrix ~ 1, k = K_RANGE_BIC, model = FLXMCmvpois(), 
                nrep = FLEXMIX_NREP, verbose = FALSE,
                control = list(iter.max = MAX_EM_ITER, tol = EM_TOL, verbose=0)) 
  }, error = function(e) { NULL })
  time_flexmix_bic_end <- Sys.time()
  time_flexmix_bic <- as.numeric(difftime(time_flexmix_bic_end, time_flexmix_bic_start, units = "secs"))
  
  k_selected_flexmix <- NA_integer_
  ari_flexmix_bic <- NA_real_
  best_flexmix_model_bic <- NULL
  
  if (!is.null(flexmix_bic_model_step)) {
    best_flexmix_model_bic <- tryCatch(getModel(flexmix_bic_model_step, "BIC"), error = function(e) NULL)
    if (!is.null(best_flexmix_model_bic)) {
      k_selected_flexmix <- best_flexmix_model_bic@k 
      cl_flex_bic <- tryCatch(clusters(best_flexmix_model_bic), error = function(e) NULL)
      if (!is.null(cl_flex_bic)) {
        ari_flexmix_bic <- tryCatch(adjustedRandIndex(true_labels, cl_flex_bic), error = function(e) NA_real_)
      }
    }
  }
  
  results_list[[length(results_list) + 1]] <- data.frame(
    NumNoisyFeatures = num_noisy,
    ARI_MyEM_FixedK = ari_myem_fixed, Time_MyEM_FixedK = time_myem_fixed,
    ARI_Flexmix_FixedK = ari_flexmix_fixed, Time_Flexmix_FixedK = time_flexmix_fixed,
    ARI_MyEM_BIC = ari_myem_bic, Time_MyEM_BIC = time_myem_bic, K_MyEM_BIC = k_selected_myem,
    ARI_Flexmix_BIC = ari_flexmix_bic, Time_Flexmix_BIC = time_flexmix_bic, K_Flexmix_BIC = k_selected_flexmix
    # Optionally, you could also store sim_data$true_dist_lambdas and sim_data$true_noise_lambda if you want to analyze their impact later
  )
}
cat("Simulation finished.\n")

results_df <- bind_rows(results_list)
results_df <- results_df %>% mutate(across(starts_with("ARI_"), as.numeric))

cat("\nSummary of results_df before plotting:\n")
print(summary(results_df)) 
cat("\n")


# --------------------------------------------------------------------------
# 5. Plotting Results (Unchanged titles, but plots reflect new data generation)
# --------------------------------------------------------------------------
cat("Generating plots...\n")

plot_ari_fixed_k <- ggplot(results_df, aes(x = NumNoisyFeatures)) +
  geom_line(aes(y = ARI_MyEM_FixedK, color = "My EM (Fixed K)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_MyEM_FixedK, color = "My EM (Fixed K)"), na.rm = TRUE) +
  geom_line(aes(y = ARI_Flexmix_FixedK, color = "Flexmix (Fixed K)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_Flexmix_FixedK, color = "Flexmix (Fixed K)"), na.rm = TRUE) +
  labs(title = paste("ARI vs. Noisy Features (True K =", TRUE_K, ", Params from Gamma Prior)"), 
       x = "Number of Noisy Features (i)", y = "Adjusted Rand Index (ARI)", color = "Algorithm") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (Fixed K)" = "blue", "Flexmix (Fixed K)" = "red"))
print(plot_ari_fixed_k)

plot_ari_bic <- ggplot(results_df, aes(x = NumNoisyFeatures)) +
  geom_line(aes(y = ARI_MyEM_BIC, color = "My EM (BIC)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_MyEM_BIC, color = "My EM (BIC)"), na.rm = TRUE) +
  geom_line(aes(y = ARI_Flexmix_BIC, color = "Flexmix (BIC)"), na.rm = TRUE) +
  geom_point(aes(y = ARI_Flexmix_BIC, color = "Flexmix (BIC)"), na.rm = TRUE) +
  labs(title = "ARI vs. Noisy Features (BIC Selected K, Params from Gamma Prior)", 
       x = "Number of Noisy Features (i)", y = "Adjusted Rand Index (ARI)", color = "Algorithm") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange"))
print(plot_ari_bic)

results_df_plot_k <- results_df %>%
  mutate(
    K_MyEM_BIC_plot = ifelse(is.na(K_MyEM_BIC), min(K_RANGE_BIC, na.rm=TRUE)-0.5, K_MyEM_BIC), 
    K_Flexmix_BIC_plot = ifelse(is.na(K_Flexmix_BIC), min(K_RANGE_BIC, na.rm=TRUE)-0.5, K_Flexmix_BIC)
  )
results_df_k_long <- results_df_plot_k %>%
  select(NumNoisyFeatures, K_MyEM_BIC_plot, K_Flexmix_BIC_plot) %>%
  rename(K_MyEM_BIC = K_MyEM_BIC_plot, K_Flexmix_BIC = K_Flexmix_BIC_plot) %>%
  pivot_longer(cols = c(K_MyEM_BIC, K_Flexmix_BIC), names_to = "Algorithm", values_to = "SelectedK") %>%
  mutate(Algorithm = case_when(
    Algorithm == "K_MyEM_BIC" ~ "My EM (BIC)",
    Algorithm == "K_Flexmix_BIC" ~ "Flexmix (BIC)", TRUE ~ Algorithm ))

plot_selected_k <- ggplot(results_df_k_long, aes(x = NumNoisyFeatures, y = SelectedK, color = Algorithm)) +
  geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  geom_hline(yintercept = TRUE_K, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(results_df$NumNoisyFeatures, na.rm=TRUE)*0.9, y = TRUE_K + 0.2, label = paste("True K =", TRUE_K), color="gray30") +
  labs(title = "Selected Number of Clusters (K) by BIC (Params from Gamma Prior)", 
       x = "Number of Noisy Features (i)", y = "Selected K", color = "Algorithm") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange")) +
  scale_y_continuous(breaks = (min(K_RANGE_BIC,na.rm=TRUE)-1):max(K_RANGE_BIC,na.rm=TRUE), 
                     labels = function(x) ifelse(x < min(K_RANGE_BIC, na.rm=TRUE), "NA/Fail", x))
print(plot_selected_k)

results_df_time_long <- results_df %>%
  select(NumNoisyFeatures, Time_MyEM_FixedK, Time_Flexmix_FixedK, Time_MyEM_BIC, Time_Flexmix_BIC) %>%
  pivot_longer(cols = -NumNoisyFeatures, names_to = "Algorithm_Metric", values_to = "Time") %>%
  mutate( Algorithm = case_when(
    grepl("MyEM_FixedK", Algorithm_Metric) ~ "My EM (Fixed K)",
    grepl("Flexmix_FixedK", Algorithm_Metric) ~ "Flexmix (Fixed K)",
    grepl("MyEM_BIC", Algorithm_Metric) ~ "My EM (BIC)",
    grepl("Flexmix_BIC", Algorithm_Metric) ~ "Flexmix (BIC)", TRUE ~ Algorithm_Metric ))

plot_time <- ggplot(results_df_time_long, aes(x = NumNoisyFeatures, y = Time, color = Algorithm)) +
  geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
  labs(title = "Execution Time vs. Number of Noisy Features (Params from Gamma Prior)", 
       x = "Number of Noisy Features (i)", y = "Time (seconds)", color = "Algorithm & Mode") +
  theme_minimal() +
  scale_color_manual(values = c("My EM (Fixed K)" = "blue", "Flexmix (Fixed K)" = "red",
                                "My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange")) 
print(plot_time)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_filename <- paste0("poisson_em_flexmix_comparison_results_", timestamp, ".csv")
# write.csv(results_df, results_filename, row.names = FALSE)
# cat(paste("Results saved to:", results_filename, "\n"))

cat("Script finished. Check plots.\n")

#---------------------------------------------------------------------------------------------------------------------------------------------------



































# 0. Libraries and Setup
# Ensure these packages are installed: install.packages(c("flexmix", "mclust", "ggplot2", "dplyr", "tidyr"))
library(flexmix)
library(mclust) # For adjustedRandIndex
library(ggplot2)
library(dplyr)
library(tidyr) # For pivot_longer

set.seed(123) # For reproducibility

# --------------------------------------------------------------------------
# 1. Custom Poisson EM Algorithm Implementation with Gamma Prior for Lambdas
#    (Functions: log_sum_exp, poisson_em_single_run, poisson_em_multistart - these remain the same as your last version)
# --------------------------------------------------------------------------
log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

poisson_em_single_run <- function(data, k, gamma_prior_shape, gamma_prior_rate, 
                                  max_iter = 100, tol = 1e-6, epsilon_lambda = 1e-6) {
  n_samples <- nrow(data)
  n_features <- ncol(data)
  
  initial_assignments <- sample(1:k, n_samples, replace = TRUE)
  pi_k <- as.numeric(table(factor(initial_assignments, levels = 1:k))) / n_samples
  pi_k[pi_k == 0] <- epsilon_lambda 
  pi_k <- pi_k / sum(pi_k)
  
  lambdas_jk <- matrix(epsilon_lambda, nrow = k, ncol = n_features)
  for (j in 1:k) {
    if (sum(initial_assignments == j) > 0) {
      cluster_data <- data[initial_assignments == j, , drop = FALSE]
      lambdas_jk[j, ] <- pmax(epsilon_lambda, colMeans(cluster_data, na.rm = TRUE))
    } else { 
      prior_mean_lambda = gamma_prior_shape / gamma_prior_rate
      lambdas_jk[j, ] <- pmax(epsilon_lambda, rep(prior_mean_lambda, n_features) + rnorm(n_features, 0, 0.1*prior_mean_lambda))
    }
  }
  lambdas_jk[is.na(lambdas_jk)] <- epsilon_lambda 
  lambdas_jk[lambdas_jk <= 0] <- epsilon_lambda 
  
  log_likelihood_old <- -Inf
  responsibilities_ij <- matrix(0, nrow = n_samples, ncol = k)
  
  for (iter in 1:max_iter) {
    log_prob_data_given_cluster <- matrix(0, nrow = n_samples, ncol = k)
    for (j in 1:k) {
      log_prob_data_given_cluster[, j] <- rowSums(
        sapply(1:n_features, function(f) dpois(data[, f], lambdas_jk[j, f], log = TRUE)),
        na.rm = TRUE 
      )
    }
    
    weighted_log_prob <- sweep(log_prob_data_given_cluster, 2, log(pi_k), "+")
    
    for (i in 1:n_samples) {
      log_sum_exp_row_i <- log_sum_exp(weighted_log_prob[i, ])
      if(is.finite(log_sum_exp_row_i)){ 
        responsibilities_ij[i, ] <- exp(weighted_log_prob[i, ] - log_sum_exp_row_i)
      } else { 
        responsibilities_ij[i, ] <- 1/k 
      }
    }
    row_sums_resp <- rowSums(responsibilities_ij, na.rm = TRUE)
    problematic_rows <- (row_sums_resp < epsilon_lambda | is.na(row_sums_resp) | is.nan(row_sums_resp))
    if(any(problematic_rows)){
      responsibilities_ij[problematic_rows, ] <- 1/k 
    }
    responsibilities_ij <- responsibilities_ij / rowSums(responsibilities_ij, na.rm = TRUE) 
    responsibilities_ij[is.nan(responsibilities_ij)] <- 1/k
    
    N_k <- colSums(responsibilities_ij, na.rm = TRUE) 
    pi_k_new <- N_k / n_samples
    pi_k_new[N_k < epsilon_lambda] <- epsilon_lambda 
    pi_k_new <- pi_k_new / sum(pi_k_new)
    
    lambdas_jk_new <- matrix(epsilon_lambda, nrow = k, ncol = n_features)
    for (j in 1:k) {
      if (N_k[j] > epsilon_lambda) { 
        weighted_sum_data_f <- colSums(responsibilities_ij[, j] * data, na.rm = TRUE) 
        lambdas_jk_new[j, ] <- (gamma_prior_shape + weighted_sum_data_f) / (gamma_prior_rate + N_k[j])
        lambdas_jk_new[j, ] <- pmax(epsilon_lambda, lambdas_jk_new[j, ]) 
      } else { 
        prior_mean_lambda = gamma_prior_shape / gamma_prior_rate
        lambdas_jk_new[j, ] <- pmax(epsilon_lambda, rep(prior_mean_lambda, n_features) + rnorm(n_features, 0, 0.01*prior_mean_lambda))
      }
    }
    lambdas_jk_new[is.na(lambdas_jk_new)] <- epsilon_lambda 
    lambdas_jk_new[lambdas_jk_new <= 0] <- epsilon_lambda 
    
    pi_k <- pi_k_new
    lambdas_jk <- lambdas_jk_new
    
    current_log_likelihood <- 0
    log_prob_data_given_cluster_updated <- matrix(0, nrow = n_samples, ncol = k)
    for (j in 1:k) {
      log_prob_data_given_cluster_updated[, j] <- rowSums(
        sapply(1:n_features, function(f) dpois(data[, f], lambdas_jk[j, f], log = TRUE)),
        na.rm = TRUE
      )
    }
    for (i in 1:n_samples) {
      log_likelihood_val_i <- log_sum_exp(log(pi_k) + log_prob_data_given_cluster_updated[i, ])
      if(is.finite(log_likelihood_val_i)){
        current_log_likelihood <- current_log_likelihood + log_likelihood_val_i
      } else {
        current_log_likelihood <- current_log_likelihood - 1e10 
      }
    }
    
    if (is.na(current_log_likelihood) || !is.finite(current_log_likelihood)) { 
      current_log_likelihood <- -Inf 
      break 
    }
    if (abs(current_log_likelihood - log_likelihood_old) / (abs(log_likelihood_old) + 1e-3) < tol && iter > 1) { 
      break
    }
    log_likelihood_old <- current_log_likelihood
  }
  
  cluster_assignments <- apply(responsibilities_ij, 1, which.max)
  if(any(tabulate(cluster_assignments, nbins=k) == 0) && k > 1){ 
    current_log_likelihood <- -Inf 
  }
  
  return(list(
    lambdas = lambdas_jk, pi = pi_k, responsibilities = responsibilities_ij,
    log_likelihood = current_log_likelihood, cluster_assignments = cluster_assignments, iterations = iter
  ))
}

poisson_em_multistart <- function(data, k, gamma_prior_shape, gamma_prior_rate, n_inits = 5, ...) {
  best_log_likelihood <- -Inf
  best_result <- NULL
  
  if (k <= 0) stop("Number of clusters k must be at least 1.")
  if (k == 1) { 
    n_samples <- nrow(data)
    n_features <- ncol(data)
    sum_data_f <- colSums(data, na.rm = TRUE) 
    lambdas_1k_f <- (gamma_prior_shape + sum_data_f) / (gamma_prior_rate + n_samples)
    lambdas_jk <- matrix(pmax(1e-6, lambdas_1k_f), nrow = 1, ncol = n_features)
    lambdas_jk[is.na(lambdas_jk)] <- 1e-6 
    pi_k <- 1
    responsibilities_ij <- matrix(1, nrow = n_samples, ncol = 1)
    log_likelihood <- sum(sapply(1:n_samples, function(sample_idx) {
      sum(dpois(data[sample_idx,], lambdas_jk[1,], log = TRUE), na.rm=TRUE) 
    }), na.rm=TRUE)
    if(is.na(log_likelihood) || !is.finite(log_likelihood)) log_likelihood <- -Inf 
    
    return(list(
      lambdas = lambdas_jk, pi = pi_k, responsibilities = responsibilities_ij,
      log_likelihood = log_likelihood, cluster_assignments = rep(1, n_samples), iterations = 1 
    ))
  }
  
  for (i_init in 1:n_inits) {
    current_result <- tryCatch({
      poisson_em_single_run(data, k, gamma_prior_shape, gamma_prior_rate, ...)
    }, error = function(e) {
      list(log_likelihood = -Inf, cluster_assignments = rep(1, nrow(data))) 
    })
    
    if (!is.null(current_result$log_likelihood) && 
        !is.na(current_result$log_likelihood) &&
        is.finite(current_result$log_likelihood) && 
        current_result$log_likelihood > best_log_likelihood) {
      best_log_likelihood <- current_result$log_likelihood
      best_result <- current_result
    }
  }
  if (is.null(best_result)) { 
    return(list( 
      lambdas = matrix(gamma_prior_shape/gamma_prior_rate, nrow = k, ncol = ncol(data)),
      pi = rep(1/k, k), responsibilities = matrix(1/k, nrow = nrow(data), ncol = k),
      log_likelihood = -Inf, cluster_assignments = rep(1, nrow(data)), iterations = 0
    ))
  }
  return(best_result)
}


# --------------------------------------------------------------------------
# 2. Data Generation Function (MODIFIED to use FIXED true lambdas)
# --------------------------------------------------------------------------
generate_poisson_mixture_data <- function(n_samples_per_cluster,
                                          n_features,
                                          num_noisy_features,
                                          true_k,
                                          fixed_lambdas_dist_pattern, # Use pre-sampled distinguishing lambdas
                                          fixed_lambda_noise_value) {  # Use pre-sampled noise lambda
  
  if (num_noisy_features < 0 || num_noisy_features > n_features) {
    stop("num_noisy_features must be between 0 and n_features.")
  }
  if (length(fixed_lambdas_dist_pattern) != true_k) {
    stop("Length of fixed_lambdas_dist_pattern must match true_k.")
  }
  
  n_total_samples <- n_samples_per_cluster * true_k
  true_labels <- rep(1:true_k, each = n_samples_per_cluster)
  data_matrix <- matrix(0, nrow = n_total_samples, ncol = n_features)
  num_dist_features <- n_features - num_noisy_features
  all_feature_indices <- 1:n_features
  if (num_noisy_features > 0) {
    noisy_feature_indices <- sample(all_feature_indices, num_noisy_features)
  } else {
    noisy_feature_indices <- c()
  }
  dist_feature_indices <- setdiff(all_feature_indices, noisy_feature_indices)
  
  for (i in 1:n_total_samples) {
    current_cluster <- true_labels[i]
    if (num_dist_features > 0 && length(dist_feature_indices) > 0) { 
      lambda_for_dist_feature <- fixed_lambdas_dist_pattern[current_cluster] 
      data_matrix[i, dist_feature_indices] <- rpois(length(dist_feature_indices), lambda_for_dist_feature)
    }
    if (num_noisy_features > 0 && length(noisy_feature_indices) > 0) { 
      data_matrix[i, noisy_feature_indices] <- rpois(length(noisy_feature_indices), fixed_lambda_noise_value)
    }
  }
  data_matrix[is.na(data_matrix)] <- 0 # Should not happen but good practice
  return(list(data = data_matrix, true_labels = true_labels))
}


# --------------------------------------------------------------------------
# 3. Simulation Parameters
# --------------------------------------------------------------------------
N_FEATURES <- 100 
N_SAMPLES_PER_CLUSTER <- 100 
TRUE_K <- 3
EPSILON_LAMBDA_GLOBAL <- 1e-6 # Global epsilon for flooring lambdas

# --- Alpha and Beta values for EM Prior to iterate over ---
ALPHA_VALUES_TO_TEST <- c(0.5, 1.0, 2.0, 5.0, 10.0) # Example values
BETA_VALUES_TO_TEST  <- c(0.5, 1.0, 2.0, 5.0, 10.0) # Example values
# This will result in 5x5 = 25 full simulation runs. Adjust as needed.

# --- FIXED True Lambdas for Data Generation (Sampled ONCE) ---
# These are sampled from a Gamma distribution with *fixed* hyperparameters,
# e.g., shape=2, rate=1 (mean=2) for distinguishing, shape=1, rate=1 (mean=1) for noise, or any other choice.
# This ensures the underlying data characteristics are the same for all alpha/beta prior tests.
DATA_GEN_DIST_LAMBDA_SHAPE <- 2.0 
DATA_GEN_DIST_LAMBDA_RATE  <- 0.5  # Mean = 4
DATA_GEN_NOISE_LAMBDA_SHAPE <- 1.0
DATA_GEN_NOISE_LAMBDA_RATE  <- 1.0  # Mean = 1

FIXED_LAMBDAS_DIST_PATTERN <- pmax(EPSILON_LAMBDA_GLOBAL, 
                                   rgamma(TRUE_K, shape = DATA_GEN_DIST_LAMBDA_SHAPE, rate = DATA_GEN_DIST_LAMBDA_RATE))
FIXED_LAMBDA_NOISE_VALUE   <- pmax(EPSILON_LAMBDA_GLOBAL, 
                                   rgamma(1, shape = DATA_GEN_NOISE_LAMBDA_SHAPE, rate = DATA_GEN_NOISE_LAMBDA_RATE))

cat("Using FIXED True Lambdas for Data Generation:\n")
cat("Distinguishing Lambdas Pattern:", FIXED_LAMBDAS_DIST_PATTERN, "\n")
cat("Noise Lambda Value:", FIXED_LAMBDA_NOISE_VALUE, "\n\n")


N_EM_INITS <- 5 
FLEXMIX_NREP <- 5 
MAX_EM_ITER <- 100
EM_TOL <- 1e-5 
K_RANGE_BIC <- 2:5 

# --- File for saving all results ---
RESULTS_CSV_FILENAME <- paste0("full_simulation_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
cat("Results will be saved to:", RESULTS_CSV_FILENAME, "\n")
# Initialize master list to store data frames from each alpha-beta pair
all_results_dfs_list <- list()


# --------------------------------------------------------------------------
# 4. Main Simulation Loop (Iterating over Alpha, Beta, and NumNoisy)
# --------------------------------------------------------------------------
cat("Starting full simulation...\n")

for (current_alpha_prior in ALPHA_VALUES_TO_TEST) {
  for (current_beta_prior in BETA_VALUES_TO_TEST) {
    
    cat(sprintf("\n\n----- RUNNING SIMULATION FOR EM PRIOR: ALPHA = %.2f, BETA = %.2f -----\n\n", 
                current_alpha_prior, current_beta_prior))
    
    current_alpha_beta_results_list <- list() # Results for this specific alpha-beta pair
    
    for (num_noisy in 0:N_FEATURES) { 
      cat(sprintf("Processing: Alpha=%.2f, Beta=%.2f, Noisy Features=%d/%d\n", 
                  current_alpha_prior, current_beta_prior, num_noisy, N_FEATURES))
      
      # Generate data using FIXED true lambdas
      sim_data <- generate_poisson_mixture_data(
        n_samples_per_cluster = N_SAMPLES_PER_CLUSTER,
        n_features = N_FEATURES,
        num_noisy_features = num_noisy,
        true_k = TRUE_K,
        fixed_lambdas_dist_pattern = FIXED_LAMBDAS_DIST_PATTERN,
        fixed_lambda_noise_value = FIXED_LAMBDA_NOISE_VALUE
      )
      data_matrix <- sim_data$data
      true_labels <- sim_data$true_labels
      
      # --- My EM (Fixed K) using current_alpha_prior and current_beta_prior ---
      time_myem_fixed_start <- Sys.time()
      myem_fixed_res <- poisson_em_multistart(data_matrix, k = TRUE_K, 
                                              gamma_prior_shape = current_alpha_prior, 
                                              gamma_prior_rate = current_beta_prior,
                                              n_inits = N_EM_INITS, 
                                              max_iter = MAX_EM_ITER, tol = EM_TOL, 
                                              epsilon_lambda = EPSILON_LAMBDA_GLOBAL)
      time_myem_fixed_end <- Sys.time()
      time_myem_fixed <- as.numeric(difftime(time_myem_fixed_end, time_myem_fixed_start, units = "secs"))
      ari_myem_fixed <- if (!is.null(myem_fixed_res) && !is.null(myem_fixed_res$cluster_assignments) && is.finite(myem_fixed_res$log_likelihood) && myem_fixed_res$log_likelihood > -Inf) {
        tryCatch(adjustedRandIndex(true_labels, myem_fixed_res$cluster_assignments), error = function(e) NA_real_)
      } else { NA_real_ }
      
      # --- Flexmix (Fixed K) - remains the same ---
      time_flexmix_fixed_start <- Sys.time()
      flexmix_fixed_model <- tryCatch({
        flexmix(data_matrix ~ 1, k = TRUE_K, model = FLXMCmvpois(), 
                control = list(nrep = FLEXMIX_NREP, iter.max = MAX_EM_ITER, tol = EM_TOL, verbose=0))
      }, error = function(e) { NULL })
      time_flexmix_fixed_end <- Sys.time()
      time_flexmix_fixed <- as.numeric(difftime(time_flexmix_fixed_end, time_flexmix_fixed_start, units = "secs"))
      ari_flexmix_fixed <- if (!is.null(flexmix_fixed_model)) {
        tryCatch(adjustedRandIndex(true_labels, clusters(flexmix_fixed_model)), error = function(e) NA_real_)
      } else { NA_real_ }
      
      # --- My EM (BIC) using current_alpha_prior and current_beta_prior ---
      time_myem_bic_start <- Sys.time()
      best_bic_myem <- Inf; best_k_myem <- NA_integer_; final_myem_model_bic <- NULL
      for (k_test in K_RANGE_BIC) {
        if (k_test > nrow(data_matrix)) next 
        current_myem_res_bic <- poisson_em_multistart(data_matrix, k = k_test, 
                                                      gamma_prior_shape = current_alpha_prior,
                                                      gamma_prior_rate = current_beta_prior,
                                                      n_inits = N_EM_INITS, max_iter = MAX_EM_ITER, 
                                                      tol = EM_TOL, epsilon_lambda = EPSILON_LAMBDA_GLOBAL)
        if (!is.null(current_myem_res_bic) && !is.na(current_myem_res_bic$log_likelihood) && is.finite(current_myem_res_bic$log_likelihood) && current_myem_res_bic$log_likelihood > -Inf ) {
          n_params <- (k_test - 1) + (k_test * N_FEATURES) 
          bic_val <- -2 * current_myem_res_bic$log_likelihood + n_params * log(nrow(data_matrix))
          if (!is.na(bic_val) && is.finite(bic_val) && bic_val < best_bic_myem) {
            best_bic_myem <- bic_val; best_k_myem <- k_test; final_myem_model_bic <- current_myem_res_bic
          }
        }
      }
      time_myem_bic_end <- Sys.time()
      time_myem_bic <- as.numeric(difftime(time_myem_bic_end, time_myem_bic_start, units = "secs"))
      k_selected_myem <- best_k_myem
      ari_myem_bic <- if (!is.null(final_myem_model_bic) && !is.null(final_myem_model_bic$cluster_assignments)) {
        tryCatch(adjustedRandIndex(true_labels, final_myem_model_bic$cluster_assignments), error = function(e) NA_real_)
      } else { NA_real_ }
      
      # --- Flexmix (BIC) - remains the same ---
      time_flexmix_bic_start <- Sys.time()
      flexmix_bic_model_step <- tryCatch({
        stepFlexmix(data_matrix ~ 1, k = K_RANGE_BIC, model = FLXMCmvpois(), nrep = FLEXMIX_NREP, verbose = FALSE,
                    control = list(iter.max = MAX_EM_ITER, tol = EM_TOL, verbose=0)) 
      }, error = function(e) { NULL })
      time_flexmix_bic_end <- Sys.time()
      time_flexmix_bic <- as.numeric(difftime(time_flexmix_bic_end, time_flexmix_bic_start, units = "secs"))
      k_selected_flexmix <- NA_integer_; ari_flexmix_bic <- NA_real_; best_flexmix_model_bic <- NULL
      if (!is.null(flexmix_bic_model_step)) {
        best_flexmix_model_bic <- tryCatch(getModel(flexmix_bic_model_step, "BIC"), error = function(e) NULL)
        if (!is.null(best_flexmix_model_bic)) {
          k_selected_flexmix <- best_flexmix_model_bic@k 
          cl_flex_bic <- tryCatch(clusters(best_flexmix_model_bic), error = function(e) NULL)
          if (!is.null(cl_flex_bic)) {
            ari_flexmix_bic <- tryCatch(adjustedRandIndex(true_labels, cl_flex_bic), error = function(e) NA_real_)
          }
        }
      }
      
      # Store results for this num_noisy iteration, including current alpha and beta for EM prior
      current_alpha_beta_results_list[[length(current_alpha_beta_results_list) + 1]] <- data.frame(
        Alpha_Prior_EM = current_alpha_prior, Beta_Prior_EM = current_beta_prior, # Add current prior params
        NumNoisyFeatures = num_noisy,
        ARI_MyEM_FixedK = ari_myem_fixed, Time_MyEM_FixedK = time_myem_fixed,
        ARI_Flexmix_FixedK = ari_flexmix_fixed, Time_Flexmix_FixedK = time_flexmix_fixed,
        ARI_MyEM_BIC = ari_myem_bic, Time_MyEM_BIC = time_myem_bic, K_MyEM_BIC = k_selected_myem,
        ARI_Flexmix_BIC = ari_flexmix_bic, Time_Flexmix_BIC = time_flexmix_bic, K_Flexmix_BIC = k_selected_flexmix
      )
    } # End of num_noisy loop
    
    # --- Process results for the current alpha-beta pair ---
    if (length(current_alpha_beta_results_list) > 0) {
      results_df_current_pair <- bind_rows(current_alpha_beta_results_list)
      results_df_current_pair <- results_df_current_pair %>% mutate(across(starts_with("ARI_"), as.numeric))
      
      # Append to master list (optional, if you want all in memory at end)
      all_results_dfs_list[[paste0("alpha",current_alpha_prior,"_beta",current_beta_prior)]] <- results_df_current_pair
      
      # Save to CSV (append)
      write.table(results_df_current_pair, file = RESULTS_CSV_FILENAME, 
                  append = TRUE, sep = ",", 
                  row.names = FALSE, 
                  col.names = !file.exists(RESULTS_CSV_FILENAME) || file.info(RESULTS_CSV_FILENAME)$size == 0)
      cat(sprintf("Results for Alpha=%.2f, Beta=%.2f appended to %s\n", 
                  current_alpha_prior, current_beta_prior, RESULTS_CSV_FILENAME))
      
      # --- Plotting for the current alpha-beta pair ---
      cat(sprintf("Generating plots for Alpha=%.2f, Beta=%.2f...\n", current_alpha_prior, current_beta_prior))
      
      plot_title_suffix <- sprintf("(EM Prior: Alpha=%.2f, Beta=%.2f | Data Lambdas Fixed)", 
                                   current_alpha_prior, current_beta_prior)
      
      plot_ari_fixed_k <- ggplot(results_df_current_pair, aes(x = NumNoisyFeatures)) +
        geom_line(aes(y = ARI_MyEM_FixedK, color = "My EM (Fixed K)"), na.rm = TRUE) +
        geom_point(aes(y = ARI_MyEM_FixedK, color = "My EM (Fixed K)"), na.rm = TRUE) +
        geom_line(aes(y = ARI_Flexmix_FixedK, color = "Flexmix (Fixed K)"), na.rm = TRUE) +
        geom_point(aes(y = ARI_Flexmix_FixedK, color = "Flexmix (Fixed K)"), na.rm = TRUE) +
        labs(title = paste("ARI vs. Noisy Features (True K =", TRUE_K, ")", plot_title_suffix), 
             x = "Number of Noisy Features (i)", y = "Adjusted Rand Index (ARI)", color = "Algorithm") +
        theme_minimal() + ylim(0,1) +
        scale_color_manual(values = c("My EM (Fixed K)" = "blue", "Flexmix (Fixed K)" = "red"))
      print(plot_ari_fixed_k)
      
      plot_ari_bic <- ggplot(results_df_current_pair, aes(x = NumNoisyFeatures)) +
        geom_line(aes(y = ARI_MyEM_BIC, color = "My EM (BIC)"), na.rm = TRUE) +
        geom_point(aes(y = ARI_MyEM_BIC, color = "My EM (BIC)"), na.rm = TRUE) +
        geom_line(aes(y = ARI_Flexmix_BIC, color = "Flexmix (BIC)"), na.rm = TRUE) +
        geom_point(aes(y = ARI_Flexmix_BIC, color = "Flexmix (BIC)"), na.rm = TRUE) +
        labs(title = paste("ARI vs. Noisy Features (BIC Selected K)", plot_title_suffix), 
             x = "Number of Noisy Features (i)", y = "Adjusted Rand Index (ARI)", color = "Algorithm") +
        theme_minimal() + ylim(0,1) +
        scale_color_manual(values = c("My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange"))
      print(plot_ari_bic)
      
      results_df_plot_k <- results_df_current_pair %>%
        mutate(
          K_MyEM_BIC_plot = ifelse(is.na(K_MyEM_BIC), min(K_RANGE_BIC, na.rm=TRUE)-0.5, K_MyEM_BIC), 
          K_Flexmix_BIC_plot = ifelse(is.na(K_Flexmix_BIC), min(K_RANGE_BIC, na.rm=TRUE)-0.5, K_Flexmix_BIC)
        )
      results_df_k_long <- results_df_plot_k %>%
        select(NumNoisyFeatures, K_MyEM_BIC_plot, K_Flexmix_BIC_plot) %>%
        rename(K_MyEM_BIC = K_MyEM_BIC_plot, K_Flexmix_BIC = K_Flexmix_BIC_plot) %>%
        pivot_longer(cols = c(K_MyEM_BIC, K_Flexmix_BIC), names_to = "Algorithm", values_to = "SelectedK") %>%
        mutate(Algorithm = case_when(
          Algorithm == "K_MyEM_BIC" ~ "My EM (BIC)",
          Algorithm == "K_Flexmix_BIC" ~ "Flexmix (BIC)", TRUE ~ Algorithm ))
      plot_selected_k <- ggplot(results_df_k_long, aes(x = NumNoisyFeatures, y = SelectedK, color = Algorithm)) +
        geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
        geom_hline(yintercept = TRUE_K, linetype = "dashed", color = "gray50") +
        annotate("text", x = max(results_df_current_pair$NumNoisyFeatures, na.rm=TRUE)*0.9, y = TRUE_K + 0.2, label = paste("True K =", TRUE_K), color="gray30") +
        labs(title = paste("Selected K by BIC", plot_title_suffix), 
             x = "Number of Noisy Features (i)", y = "Selected K", color = "Algorithm") +
        theme_minimal() +
        scale_color_manual(values = c("My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange")) +
        scale_y_continuous(breaks = (min(K_RANGE_BIC,na.rm=TRUE)-1):max(K_RANGE_BIC,na.rm=TRUE), 
                           labels = function(x) ifelse(x < min(K_RANGE_BIC, na.rm=TRUE), "NA/Fail", x))
      print(plot_selected_k)
      
      results_df_time_long <- results_df_current_pair %>%
        select(NumNoisyFeatures, Time_MyEM_FixedK, Time_Flexmix_FixedK, Time_MyEM_BIC, Time_Flexmix_BIC) %>%
        pivot_longer(cols = -NumNoisyFeatures, names_to = "Algorithm_Metric", values_to = "Time") %>%
        mutate( Algorithm = case_when(
          grepl("MyEM_FixedK", Algorithm_Metric) ~ "My EM (Fixed K)",
          grepl("Flexmix_FixedK", Algorithm_Metric) ~ "Flexmix (Fixed K)",
          grepl("MyEM_BIC", Algorithm_Metric) ~ "My EM (BIC)",
          grepl("Flexmix_BIC", Algorithm_Metric) ~ "Flexmix (BIC)", TRUE ~ Algorithm_Metric ))
      plot_time <- ggplot(results_df_time_long, aes(x = NumNoisyFeatures, y = Time, color = Algorithm)) +
        geom_line(na.rm = TRUE) + geom_point(na.rm = TRUE) +
        labs(title = paste("Execution Time vs. Noisy Features", plot_title_suffix), 
             x = "Number of Noisy Features (i)", y = "Time (seconds)", color = "Algorithm & Mode") +
        theme_minimal() +
        scale_color_manual(values = c("My EM (Fixed K)" = "blue", "Flexmix (Fixed K)" = "red",
                                      "My EM (BIC)" = "darkgreen", "Flexmix (BIC)" = "orange")) 
      print(plot_time)
    } else {
      cat(sprintf("No results generated for Alpha=%.2f, Beta=%.2f. Skipping save and plot.\n",
                  current_alpha_prior, current_beta_prior))
    }
  } # End of current_beta_prior loop
} # End of current_alpha_prior loop

cat("\n\nFull simulation across all Alpha/Beta pairs finished.\n")
cat("All results appended to:", RESULTS_CSV_FILENAME, "\n")
