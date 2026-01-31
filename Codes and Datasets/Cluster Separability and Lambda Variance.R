# Load libraries
library(ggplot2)
library(gridExtra)
library(reshape2)
library(MCMCpack)  # for Dirichlet

set.seed(42)

# Parameters
n <- 500     # number of samples
k <- 5       # number of clusters
d <- 10      # number of features
alpha <- 10   # Gamma shape
beta <- 0.02  # Gamma rate -> high mean = 400, high variance

# Step 1: Dirichlet mixing proportions
pi <- rdirichlet(1, rep(1, k))[1, ]

# Step 2: Cluster assignments
z <- sample(1:k, n, replace = TRUE, prob = pi)

# Step 3: Generate λ (k x d) from Gamma prior
lambda <- matrix(rgamma(k * d, shape = alpha, rate = beta), nrow = k, ncol = d)

# Step 4: Generate Poisson data
X <- matrix(0, nrow = n, ncol = d)
for (i in 1:n) {
  X[i, ] <- rpois(d, lambda[z[i], ])
}

# Create dataframes for plotting
X_unsorted <- X
X_sorted <- X[order(z), ]

# Melt for ggplot
df_unsorted <- melt(X_unsorted)
df_sorted <- melt(X_sorted)

colnames(df_unsorted) <- c("Sample", "Feature", "Count")
colnames(df_sorted) <- c("Sample", "Feature", "Count")

# Plot a: Unsorted
p1 <- ggplot(df_unsorted, aes(x = Feature, y = Sample, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis") +
  labs(title = "Heatmap of Poisson data (unsorted)",
       subtitle = "(a) Unsorted data (high variance)",
       x = "Feature", y = "Sample") +
  theme_minimal() +
  theme(plot.subtitle = element_text(hjust = 0.5, face = "italic"))

# Plot b: Sorted by cluster
p2 <- ggplot(df_sorted, aes(x = Feature, y = Sample, fill = Count)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis") +
  labs(title = "Heatmap of Poisson data (samples sorted by cluster)",
       subtitle = "(b) Sorted data (high variance)",
       x = "Feature", y = "Sample") +
  theme_minimal() +
  theme(plot.subtitle = element_text(hjust = 0.5, face = "italic"))

# Combine the two plots
grid.arrange(p1, p2, ncol = 2)
