# LAB 10

# I)
# Ex1
# Create a sequence from 1 to 100
x <- seq(1, 100)

# Sampling without replacement
s_no_replace <- sample(x, 10)
cat("Sampling WITHOUT replacement:\n", s_no_replace, "\n\n")

# Sampling with replacement
s_with_replace <- sample(x, 10, replace = TRUE)
cat("Sampling WITH replacement:\n", s_with_replace, "\n\n")

# Ex2
# install.packages("gtools")
library(gtools)
x<-c("A","B","C","D")
# -------- COMBINATIONS (Order does NOT matter) --------
cat("Combinations without repetition (Choose 2 elements from {1,2,3,4}):\n")
print(combinations(n = length(x), r = 3, v = x, repeats.allowed = FALSE))
# Output: All possible 2-element subsets of {1,2,3,4} where order does NOT matter

cat("\nCombinations with repetition (Choose 2 elements, repetition allowed):\n")
print(combinations(n = length(x), r = 3, v = x, repeats.allowed = TRUE))
# Output: All possible 2-element subsets allowing repeated elements

# -------- PERMUTATIONS (Order matters) --------
cat("\nPermutations without repetition (Choose 2 elements from {1,2,3,4}):\n")
print(permutations(n = length(x), r = 2, v = x, repeats.allowed = FALSE))
# Output: All possible 2-element sequences of {1,2,3,4} where order MATTERS

cat("\nPermutations with repetition (Choose 2 elements, repetition allowed):\n")
print(permutations(n = length(x), r = 2, v = x, repeats.allowed = TRUE))
# Output: All possible 2-element sequences allowing repeated elements


# II)
# Ex1
# Set parameters for the binomial distribution
n <- 10     # Number of trials
p1 <- 0.4   # Probability of success
p2 <- 0.7   # Second probability value for comparison
m <- 3      # Number of successes

# (a) Print probability value for given m, n, p1
prob_value <- dbinom(m, size = n, prob = p1)
cat("P(X = 3) for Binomial(n=10, p=0.4):", prob_value, "\n")

# (b) Print cumulative probability value for given m, n, p1
cum_prob_value <- pbinom(m, size = n, prob = p1)
cat("P(X ≤ 3) for Binomial(n=10, p=0.4):", cum_prob_value, "\n")

# (c) Find m for cumulative probability of 0.8
m_for_cum_08 <- qbinom(0.8, size = n, prob = p1)
cat("Smallest m such that P(X ≤ m) ≥ 0.8:", m_for_cum_08, "\n")

# (d) Random sampling: Generate 5 points from Binomial(n=10, p=0.4)
random_samples <- rbinom(5, size = n, prob = p1)
cat("5 Randomly Sampled Points from Binomial(n=10, p=0.4):", random_samples, "\n")

# (e) Plot PDF for p=0.4 and p=0.7
x_vals <- 0:n  # Possible values for X
pdf_p1 <- dbinom(x_vals, size = n, prob = p1)
pdf_p2 <- dbinom(x_vals, size = n, prob = p2)

# Plot the first PDF
plot(x_vals, pdf_p1, type = "l", col = "blue", lwd = 2, ylim = c(0, max(pdf_p1, pdf_p2)),
     xlab = "Number of Successes", ylab = "Probability", main = "Binomial PDF for p=0.4 and p=0.7")

# Add the second PDF to the same plot
lines(x_vals, pdf_p2, type = "l", col = "red", lwd = 2)
legend("topright", legend = c("p=0.4", "p=0.7"), col = c("blue", "red"), lwd = 2)

# Case 1: p = 0.4 (blue curve)
# This means heads (success) happens 40% of the time and tails (failure) happens 60% of the time.
# Since getting heads is less likely, you will mostly get fewer heads in 10 flips.
# The graph peaks around 4 heads (since 40% of 10 is 4), but there’s a longer tail stretching to the right because sometimes you might still get more heads.
# That’s why the graph is right-skewed—most values are on the left, but there’s a tail extending to the right.
# 
# Case 2: p = 0.7 (red curve)
# Now, heads (success) happens 70% of the time, and tails (failure) happens 30% of the time.
# Since getting heads is more likely, you will mostly get more heads in 10 flips.
# The graph peaks around 7 heads, but sometimes you still get fewer heads.
# That’s why the graph is left-skewed—most values are on the right, but there’s a tail extending to the left.
# 
# Intuitive Trick:
#   If success is rare (p < 0.5), most trials fail → Graph is right-skewed (more low values).
#   If success is common (p > 0.5), most trials succeed → Graph is left-skewed (more high values).

# (f) Generate 100 and 10,000 points and create frequency tables
set.seed(42)  # Set seed for reproducibility
sample_100 <- rbinom(100, size = n, prob = p1)
sample_10000 <- rbinom(10000, size = n, prob = p1)

freq_100 <- table(sample_100)
freq_10000 <- table(sample_10000)

# Plot frequency distributions
par(mfrow = c(2, 1))  # 2x1 grid layout

barplot(freq_100, col = "skyblue", main = "Frequency Distribution (n=100)",
        xlab = "Number of Successes", ylab = "Frequency", border = "blue")

barplot(freq_10000, col = "orange", main = "Frequency Distribution (n=10000)",
        xlab = "Number of Successes", ylab = "Frequency", border = "red")

par(mfrow = c(1, 1))  # Reset layout


# Ex2
# Set parameters for Hypergeometric Distribution
N <- 100   # Total population size
K <- 70    # Total number of "success" items in the population
n <- 12    # Sample size drawn
# p <- 0.3   # Not needed for hypergeometric, given as reference

# (a) Plot histogram of the hypergeometric probability density function
x_values <- 0:n  # Possible values of k (successes in sample)
prob_values <- dhyper(x_values, K, N - K, n)  # Compute PMF

barplot(prob_values, names.arg = x_values, col = "skyblue",
        main = "Hypergeometric Distribution (N=100, K=70, n=12)",
        xlab = "Number of Successes (k)", ylab = "Probability")

# Add text annotations for parameters
text(4, max(prob_values) * 0.9, paste("N =", N, "\nK =", K, "\nn =", n), col = "red")

# (b) Compute cumulative probability up to x = 10
cum_prob <- phyper(10, K, N - K, n)
print(paste("Cumulative probability up to x = 10:", round(cum_prob, 3)))

# (c) Find the x value for cumulative probability 0.9
x_value <- qhyper(0.9, K, N - K, n)
print(paste("X value corresponding to cumulative probability 0.9:", x_value))

# (d) Sample 5 random points from the hypergeometric distribution
set.seed(42)  # For reproducibility
samples <- rhyper(5, K, N - K, n)
print("Random samples:") 
print(signif(samples, 2)) # 2 significant digits


# Ex3

# (a)
# Set up plotting area (1x2 grid)
par(mfrow = c(1, 2))

# Define x values (trial numbers)
x_vals <- 1:10  # First 10 trials

# Compute geometric probabilities for p = 0.3
pdf_p1 <- dgeom(x_vals-1, prob = 0.3)  # P(X = x)

# Compute geometric probabilities for p = 0.8
pdf_p2 <- dgeom(x_vals-1, prob = 0.8)

# Plot for p = 0.3
barplot(pdf_p1, names.arg = x_vals, col = "blue",
        main = "Geometric PDF (p=0.3)",
        xlab = "Trial Number", ylab = "Probability")

# Plot for p = 0.8
barplot(pdf_p2, names.arg = x_vals, col = "red",
        main = "Geometric PDF (p=0.8)",
        xlab = "Trial Number", ylab = "Probability")

# Reset plotting layout
par(mfrow = c(1, 1))

# Observations:
# When p = 0.3, the probability mass is more spread out, meaning more trials are needed before the first success.
# When p = 0.8, the probability mass is more concentrated at smaller values, meaning success happens earlier in fewer trials.

# (b)
cum_prob_x4 <- pgeom(4 - 1, prob = 0.3)  # P(X ≤ 4)
print(round(cum_prob_x4, 3))  # Round to 3 decimal places

# (c)
m_value <- qgeom(0.2, prob = 0.3) + 1  # Find trial number where P(X ≤ m) = 0.2
print(m_value)

# (d)
set.seed(42)  # Ensure reproducibility
samples <- rgeom(6, prob = 0.4) + 1  # Generate 6 random samples
print(samples)


# Ex 4

# Set seed for reproducibility
set.seed(42)

# (4a) Compute and print the negative binomial probability density for y = 5, r = 3, p = 0.3
y <- 5
r <- 3
p <- 0.3
nbinom_pdf <- dnbinom(y, size = r, prob = p)
cat("Negative Binomial PDF for y =", y, ":", nbinom_pdf, "\n")

# (4b) Compute and print the cumulative negative binomial probability density up to y = 5
nbinom_cdf <- pnbinom(5, size = 3, prob = 0.3)
cat("Cumulative Probability up to y = 5:", nbinom_cdf, "\n")

# (4c) Find the y-value corresponding to a cumulative probability of 0.5 (median)
y_median <- qnbinom(0.5, size = 3, prob = 0.3)
cat("y value for cumulative probability of 0.5:", y_median, "\n")

# (4d) Print 4 random points sampled from this distribution with r = 3 and p = 0.3
nbinom_samples <- rnbinom(4, size = 3, prob = 0.3)
cat("Random samples from Negative Binomial Distribution:", nbinom_samples, "\n")

# (4e) Plot the negative binomial probability density function with r = 10, p = 0.3
y_vals <- 0:50  # Range of failures
nbinom_probs <- dnbinom(y_vals, size = 10, prob = 0.3)

barplot(nbinom_probs, names.arg = y_vals, col = "blue",
        main = "Negative Binomial Probability Distribution (r=10, p=0.3)",
        xlab = "Number of Failures", ylab = "Probability")

# (4f) Generate a frequency histogram of 10,000 random deviates from this distribution
nbinom_samples_large <- rnbinom(10000, size = 10, prob = 0.3)

hist(nbinom_samples_large, breaks = 30, col = "lightblue", border = "black",
     main = "Histogram of 10,000 Negative Binomial Deviates (r=10, p=0.3)",
     xlab = "Number of Failures", ylab = "Frequency", probability = TRUE)


# Ex 5

# Set seed for reproducibility
set.seed(42)

# (5a) Compute and print the Poisson probability given λ = 10 and m = 7
lambda <- 10
m <- 7
poisson_pmf <- dpois(m, lambda)
cat("Poisson PMF for m =", m, "and λ =", lambda, ":", poisson_pmf, "\n")

# (5b) Calculate and print the cumulative probability P(X ≤ m) for λ = 10 and m = 7
poisson_cdf <- ppois(m, lambda)
cat("Cumulative Probability P(X ≤", m, ") for λ =", lambda, ":", poisson_cdf, "\n")

# (5c) Compare a Binomial(n=1000, p=0.3) distribution with a Poisson(λ=np)
n <- 1000
p <- 0.3
lambda_binom <- n * p  # Expected mean of Poisson equivalent

# Compute probabilities
x_vals <- 1:n
binom_probs <- dbinom(x_vals, size = n, prob = p)
poisson_probs <- dpois(x_vals, lambda = lambda_binom)

par(mfrow=c(2,1))
# Plot the binomial distribution
barplot(binom_probs, names.arg = x_vals, col = "red", border = "black",
        main = "Binomial Distribution (n=1000, p=0.3)",
        xlab = "m", ylab = "Probability")

# Plot the Poisson distribution
barplot(poisson_probs, names.arg = x_vals, col = "blue", border = "black",
        main = "Poisson Approximation (λ=np)",
        xlab = "m", ylab = "Probability")
par(mfrow=c(1,1))

### overlay binomial and poisson 

# compare Binomial and Poisson Distributions
n <- 10000  # Large n for binomial
p <- 0.01
lambda_binom <- n * p  # λ = np

# Define m values (x values)
m_vals <- 0:n

# Generate Binomial and Poisson distributions
binom_pmf <- dbinom(m_vals, n, p)
poisson_pmf <- dpois(m_vals, lambda = lambda_binom)

# Plot Binomial PMF as vertical lines
plot(m_vals, binom_pmf, type = "h", col = "skyblue", lwd = 2,
     main = "Overlay: Binomial vs Poisson Distributions",
     xlab = "m (Number of Events)", ylab = "Probability",
     ylim = range(c(binom_pmf, poisson_pmf)),
     xlim=c(50,150))

# Overlay Poisson PMF as points
lines(m_vals, poisson_pmf, col = "violet", type = 'h', cex = 0.6)

# Add legend
legend("topright", legend = c("Binomial", "Poisson"),
       col = c("skyblue", "violet"), lwd = c(2, 2), bty = "n")


# (5d) Find the quantile value corresponding to cumulative probability of 0.22 for λ = 10
quantile_value <- qpois(0.22, lambda)
cat("Quantile value for cumulative probability of 0.22 with λ =", lambda, ":", quantile_value, "\n")

# (5e) Obtain 10,000 random sample points from a Poisson(λ = 9) and plot histogram
poisson_samples <- rpois(10000, lambda = 9)

hist(poisson_samples, breaks = 30, col = "lightblue", border = "black",
     main = "Histogram of 10,000 Poisson Deviates (λ=9)",
     xlab = "m", ylab = "Frequency", probability = TRUE)


# 6) Gaussian Distribution

# (a) Compute and print the unit normal PDF value for µ = 12 and σ = 2.

mu = 12
sigma = 2
x <- mu

pdf_value <- dnorm(x,mean = mu,sd=sigma)
print(paste("PDF value at X=",mu,":",round(pdf_value,4)))

# (b) Calculate and print the cumulative probability for Z = 2.0. Is this same as 1-
#     CPDF(Z=-2)?

z <- 2
cpdf_z <- pnorm(z)
cpdf_neg_z <- pnorm(-z)

print(paste("CPDF at z=",cpdf_z,":",round(cpdf_z,4)))
print(paste("CPDF at negative z=",cpdf_neg_z,":",round(cpdf_neg_z,4)))
print(paste("Checking if CPDF and neg CPDF is same or not :",cpdf_z == (1 - cpdf_neg_z)))

# (c) Plot a unit normal curve for the above parameters with X range of ±4σ and add a
# text box to the plot showing the parameter symbols and their values.


# (c) Plot a normal curve for μ = 12, σ = 2 with range ±4σ
x_vals <- seq(mu - 4*sigma, mu + 4*sigma, length.out=100)
y_vals <- dnorm(x_vals, mean=mu, sd=sigma)

plot(x_vals, y_vals, type="l", col="blue", lwd=2, main="Normal Distribution",
     xlab="X", ylab="Density")

text(mu, max(y_vals) * 0.9, labels=paste("μ =", mu, "\nσ =", sigma), col="red")

# (d) Compute the 75th quantile
q_75 <- qnorm(0.75, mean=mu, sd=sigma)
print(paste("75th Quantile:", round(q_75, 4)))

# (e) Generate 10,000 random deviates and plot histogram with normal curve overlay
set.seed(42)
samples <- rnorm(10000, mean=mu, sd=sigma)

hist(samples, breaks=50, probability=TRUE, col="lightgray", main="Histogram of Random Deviates")
curve(dnorm(x, mean=mu, sd=sigma), col="red", lwd=2, add=TRUE)

# (f) Histogram plot of a ‘normalised’ binomial distribution with μ = np = 10 and p = 0.5.
# Set parameters
# Binomial parameters
# for n value less - use less breaks, for higer n use higher breaks 
n <- 1000
p <- 0.5
mu_bin <- n * p
sigma_bin <- sqrt(n * p * (1 - p))

# Generate binomial samples
m <- rbinom(10000, size = n, prob = p)

# Compute W
W <- (m - mu_bin) / sigma_bin

# Plot histogram of W
hist(W, breaks = 15, probability = TRUE, col = "skyblue",
     main = "Normalized Binomial vs Standard Normal",
     xlab = "W", xlim = c(-4, 4))

# Overlay standard normal curve
curve(dnorm(x, mean = 0, sd = 1), from = -4, to = 4, add = TRUE, col = "red", lwd = 2)

legend("topright", legend = c("Normalized Binomial", "Standard Normal N(0,1)"),
       fill = c("skyblue", NA), border = NA, col = c("skyblue", "red"), lwd = 2,inset=0.05)


# (g) Plot Poisson PDFs for λ values 1, 10, 100, 1000 and overlay normal approximation
par(mfrow = c(2, 2))  # Set 2x2 grid

lambda_vals <- c(1, 10, 100, 1000)

for (lambda in lambda_vals) {
  # Integer m values around the mean ± 4σ
  m_vals <- seq(floor(lambda - 4 * sqrt(lambda)), ceiling(lambda + 4 * sqrt(lambda)))
  
  # Poisson PMF values
  pois_pmf <- dpois(m_vals, lambda)
  
  # Normal approximation at integer points
  norm_approx <- dnorm(m_vals, mean = lambda, sd = sqrt(lambda))
  
  # Plot Poisson PMF
  plot(m_vals, pois_pmf, type = "h", col = "blue", lwd = 2,
       main = paste("Poisson vs Normal (λ =", lambda, ")"),
       xlab = "m", ylab = "Probability", ylim = c(0, max(c(pois_pmf, norm_approx)) * 1.1))
  
  # Overlay Normal approximation as points
  points(m_vals, norm_approx, col = "red", pch = 16)
  
  # Add legend
  legend("topright", legend = c("Poisson PMF", "Normal PDF (approx)"),
         col = c("blue", "red"), lty = c(1, NA), pch = c(NA, 16), lwd = 2, bty = "n")
}

par(mfrow = c(1, 1))  # Reset plot layout


# (h) Correlated Normal Distributions using MASS
library(MASS)
tolerance=0.5
xy <- mvrnorm(1000, mu=c(50,60), Sigma=matrix(c(4,tolerance,tolerance,9),2,2))
print("Variance-Covariance Matrix:")
print(var(xy))

# Extract x and y
x <- xy[,1]
y <- xy[,2]

# Scatter plot of x and y
plot(x, y, main="Scatter Plot of Correlated Variables", xlab="X", ylab="Y", col="blue")
print(paste("Var(X) =", var(x)))
print(paste("Var(Y) =", var(y)))

# Check independence: sum of variances vs variance of sum
var_sum <- var(x) + var(y)
var_combined <- var(x + y)
print(paste("Sum of Individual Variances =", round(var_sum, 4)))
print(paste("Variance of (X+Y) =", round(var_combined, 4)))
print(paste("Are they independent?", round(var_sum, 4) == round(var_combined, 4)))

# Compute covariance using correlation coefficient
cov_computed <- cor(x, y) * sqrt(var(x) * var(y))
cov_reported <- var(xy)[1,2]

print(paste("Computed Covariance =", round(cov_computed, 4)))
print(paste("Reported Covariance =", round(cov_reported, 4)))
print(paste("Do they match?", round(cov_computed, 4) == round(cov_reported, 4)))


# Ex 7 - uniform distributions

# (a) Generate 5 uniform random numbers between 0 and 1:
set.seed(123)  # For reproducibility
runif(5)

# (b) Generate 5 random samples from a uniform distribution between 50 and 100:
runif(5, min = 50, max = 100)

# (c) Generate 10,000 uniform deviates and plot a histogram with x-limits 1 and 2:
samples <- runif(10000)  # Uniform [0, 1]

# Plot histogram with x-limits from 1 to 2
hist(samples, breaks = 50, col = "skyblue", xlab = "Value", main = "Histogram of Uniform(0,1) Samples", xlim = c(0, 1))


# Ex 8

# (a) What is the probability density corresponding to x = 3 and λ = 2?
x <- 3
lambda <- 2
density <- dexp(x, rate = lambda)
print(paste("Density at x = 3 and λ = 2:", density))

# (b) What is the quantile value corresponding to cumulative probability 0.995?
p <- 0.995
quantile_val <- qexp(p, rate = lambda)
print(paste("Quantile at p = 0.995 and λ = 2:", quantile_val))

# (c) Plot the exponential CDFs for λ = 2, 10, 100 on the same graph:
x_vals <- seq(0, 2, length.out = 500)
plot(x_vals, pexp(x_vals, rate = 2), type = "l", col = "blue", lwd = 2,
     ylab = "Cumulative Probability", xlab = "x", main = "Exponential CDFs")
lines(x_vals, pexp(x_vals, rate = 10), col = "red", lwd = 2)
lines(x_vals, pexp(x_vals, rate = 100), col = "green", lwd = 2)
legend("bottomright", legend = c("λ = 2", "λ = 10", "λ = 100"),
       col = c("blue", "red", "green"), lwd = 2)

# (d) Compute and print 4 random deviates from exponential distribution with λ = 3:
set.seed(42)  # Optional, for reproducibility
random_vals <- rexp(4, rate = 3)
print("4 Random Exponential Deviates (λ = 3):")
print(random_vals)


# Ex 9 - Gamma distribution

# Set up 1x2 plotting grid
par(mfrow = c(1, 2))

# (a) Plot PDFs for alpha = 1, 2, 3, 4 with theta = 4
x_vals <- seq(0, 50, length.out = 1000)
theta <- 4
colors <- c("black", "blue", "red", "magenta")

plot(x_vals, dgamma(x_vals, shape = 1, scale = theta), type = "l", col = colors[1],
     main = "Varying α (θ = 4)", xlab = "x", ylab = "Density", ylim = c(0, 0.15))

for (a in 2:4) {
  lines(x_vals, dgamma(x_vals, shape = a, scale = theta), col = colors[a])
}

legend("topright", legend = paste("α =", 1:4), col = colors, lwd = 2,inset = 0.05)

# (a) Plot PDFs for theta = 1, 2, 3, 4 with alpha = 4
alpha <- 4
plot(x_vals, dgamma(x_vals, shape = alpha, scale = 1), type = "l", col = colors[1],
     main = "Varying θ (α = 4)", xlab = "x", ylab = "Density", ylim = c(0, 0.4))

for (t in 2:4) {
  lines(x_vals, dgamma(x_vals, shape = alpha, scale = t), col = colors[t])
}

legend("topright", legend = paste("θ =", 1:4), col = colors, lwd = 2)


# The Gamma distribution models the waiting time for α events to occur in a Poisson process.
# α (shape parameter):
# Controls how many events are required before stopping.
# More events = More spread out distribution.
# θ (scale parameter):
# Controls the time per event.
# Larger θ = Slower process = More spread out distribution.

# Reset layout
par(mfrow = c(1, 1))

# (b) Probability density at x = 6, α = 4, θ = 1
density_val <- dgamma(6, shape = 4, scale = 1)
print(paste("PDF at x = 6, α = 4, θ = 1:", round(density_val, 5)))

# (c) Cumulative probability up to x = 6
cum_prob <- pgamma(6, shape = 4, scale = 1)
print(paste("Cumulative probability up to x = 6:", round(cum_prob, 5)))

# (d) Quantile value for 0.95 cumulative probability
quantile_val <- qgamma(0.95, shape = 4, scale = 1)
print(paste("x value at 0.95 cumulative probability:", round(quantile_val, 5)))

# (e) 10,000 random deviates and histogram
set.seed(42)
gamma_sample <- rgamma(10000, shape = 4, scale = 1)
hist(gamma_sample, breaks = 30, main = "Histogram of Gamma(α=4, θ=1)", xlab = "Value", col = "skyblue", border = "white")


# Ex 10 - Chi square distribution 

# (a) Plot χ² distribution for df = 2, 3, 5, 10
x_vals <- seq(0, 30, length.out = 1000)
df_vals <- c(2, 3, 5,10)
colors <- c("black", "blue", "red", "magenta")

plot(x_vals, dchisq(x_vals, df = df_vals[1]), type = "l", col = colors[1],
     main = "Chi-square PDFs", xlab = "x", ylab = "Density", ylim = c(0, 0.3))

for (i in 2:4) {
  lines(x_vals, dchisq(x_vals, df = df_vals[i]), col = colors[i])
}

legend("topright", legend = paste("df =", df_vals), col = colors, lwd = 2)

# Shows multiple chi-square PDFs with different degrees of freedom: 2, 3, 5, 10.
# You can observe:
# As df increases, the distribution becomes more symmetric and spreads out.
# For small df (e.g., 2), it is highly skewed right.

# (b) Probability density at x = 6, df = 5
pdf_val <- dchisq(6, df = 5)
print(paste("PDF at x = 6, df = 5:", round(pdf_val, 5)))

# (c) Cumulative probability up to x = 6, df = 10
cum_val <- pchisq(6, df = 10)
print(paste("Cumulative probability up to x = 6, df = 10:", round(cum_val, 5)))

# (d) 85th quantile for df = 6
quantile_85 <- qchisq(0.85, df = 6)
print(paste("85th percentile for df = 6:", round(quantile_85, 5)))

# (e) Histogram of 10,000 random deviates, df = 6
set.seed(42)
chisq_sample <- rchisq(10000, df = 6)
hist(chisq_sample, breaks = 30, col = "red", border = "white", 
     main = "Chi-square r = 6", xlab = "Value")
text(20, 1000, "r = 6", col = "black", cex = 1.5)

# (f) Z² = ((x - μ)^2) / σ² and χ²(1) overlay
mu <- 2
sigma <- 1
x_vals <- seq(-5, 10, length.out = 1000)
z_sq_vals <- ((x_vals - mu)^2) / sigma^2
chisq_pdf <- dchisq(z_sq_vals, df = 1)

plot(z_sq_vals, chisq_pdf, type = "l", col = "darkgreen", lwd = 2,
     main = expression(paste("Chi-square PDF with ", df == 1)),
     xlab = expression(Z^2), ylab = "Density")

# Shows a chi-squared PDF with df = 1, based on your transformation Z^2
# Very sharp peak near 0 and heavy right tail.
# This is expected for chi-squared with df = 1: it's highly skewed and concentrated near zero.

# Ex 11 - CLT
# Set seed for reproducibility
set.seed(42)

# (1) CLT with Uniform Distribution
# ----------------------------------

# Step 1: Generate 10,000 samples of 5 values each from U[0,10]
samples <- replicate(10000, runif(5, min = 0, max = 10))

# Step 2: Compute the mean of each sample
sample_means <- apply(samples, 2, mean)

# Compute the mean and SD of the sample means
mean_sample_means <- mean(sample_means)
sd_sample_means <- sd(sample_means)

# Standardize the sample means
Z_samples <- (sample_means - mean_sample_means) / sd_sample_means

# Step 3: Plot histogram of standardized means
hist_info <- hist(Z_samples,
                  breaks = 25,
                  col = "lightblue",
                  main = "Standardized Sample Means (Z-scores)",
                  xlab = "Z-score",
                  border = "white")

# Step 4: Generate normal curve based on standard normal distribution
x_seq <- seq(min(Z_samples), max(Z_samples), length.out = 200)
pdf_vals <- dnorm(x_seq, mean = 0, sd = 1)

# Step 5: Scale the normal PDF to match histogram
bin_width <- hist_info$breaks[2] - hist_info$breaks[1]
scaling_factor <- length(Z_samples) * bin_width
scaled_pdf <- pdf_vals * scaling_factor

# Step 6: Overlay normal curve
lines(x_seq, scaled_pdf, col = "red", lwd = 2)

# Step 7: Print stats
cat("Mean of sample means:", mean_sample_means, "\n")
cat("Standard deviation of sample means:", sd_sample_means, "\n")
cat("Bin width:", bin_width, "\n")


# (2) CLT with Dice Rolls
# ------------------------

# Step (i): Single dice throw (10,000 rolls)
a <- sample(1:6, size=10000, replace=TRUE)
hist(a, breaks=6, prob=TRUE, col="lightblue", main="Single Dice Roll", xlab="Dice Value", border="black")

# Step (ii): Two dice (sum of two rolls)
b <- sample(1:6, size=10000, replace=TRUE) + sample(1:6, size=10000, replace=TRUE)
hist(b, breaks=11, prob=TRUE, col="lightblue", main="Sum of Two Dice Rolls", xlab="Sum", border="black")

# Step (iii): Three dice (sum of three rolls)
c <- sample(1:6, size=10000, replace=TRUE) + sample(1:6, size=10000, replace=TRUE) + sample(1:6, size=10000, replace=TRUE)
hist(c, breaks=16, prob=TRUE, col="lightblue", main="Sum of Three Dice Rolls", xlab="Sum", border="black")

# Step (iv): Five dice (sum of five rolls)
d <- sample(1:6, size=10000, replace=TRUE) + sample(1:6, size=10000, replace=TRUE) +
  sample(1:6, size=10000, replace=TRUE) + sample(1:6, size=10000, replace=TRUE) +
  sample(1:6, size=10000, replace=TRUE)

# Compute mean and standard deviation for normal curve overlay
mean_d <- mean(d)
sd_d <- sd(d)

# Generate normal PDF with calculated mean and SD
x_vals <- seq(min(d), max(d), length=100)
normal_curve <- dnorm(x_vals, mean=mean_d, sd=sd_d)

# Scale the normal PDF
scaled_curve <- normal_curve * 10000 * 0.5  # Adjusting for histogram scaling

# Plot histogram with normal curve overlay
hist(d, breaks=30, prob=TRUE, col="lightblue", main="Sum of Five Dice Rolls", xlab="Sum", border="black")
lines(x_vals, scaled_curve, col="red", lwd=2)


# Ex 12 - ROC

# Step 1: Install and Load Required Package
# install.packages("pROC")  # Run this only if you don't have it
# Load required package
library(pROC)

# (1) Load white wine data
wine_data <- read.csv("/home/ibab/R/Lab10/winequality-white.csv", sep = ";")

# (1 continued) Create binary classification columns for thresholds 6 to 10
thresholds <- 6:10
for (t in thresholds) {
  col_name <- paste0("quality_", t)
  wine_data[[col_name]] <- ifelse(
    if (t == 10) wine_data$quality == t else wine_data$quality >= t, 1, 0)
}

# Setup for plotting multiple plots
par(mfrow = c(3, 2))  # adjust layout if needed

for (t in thresholds) {
  label_col <- paste0("quality_", t)
  labels <- wine_data[[label_col]]
  
  # Check if response has both classes (0 and 1)
  if (length(unique(labels)) < 2) {
    cat("Skipping threshold", t, "- not enough class diversity (needs both 0 and 1)\n")
    next
  }
  
  main_title <- paste("ROC Curve -", ifelse(t == 10, "Quality == 10", paste("Quality >=", t)))
  
  plot.roc(labels, wine_data$alcohol,
           main = main_title,
           legacy.axes = TRUE,
           ci = TRUE,
           print.auc = TRUE,
           identity.lwd = 2,
           print.thres = TRUE)
}
