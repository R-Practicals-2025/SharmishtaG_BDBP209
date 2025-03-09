# Ex1
plot(2,4,pch=0,col="magenta")


# Ex2
x <- seq(-pi, pi, 0.1)  # Generate x values from -π to π
y1 <- sin(x)            # Compute sine values
y2 <- cos(x)            # Compute cosine values

# Plot sine function with blue star points
plot(x, y1,type='o', col="blue", pch=8, xlab="x", ylab="y = sin(x) / y = cos(x)", 
     main="Sine vs Cos Function")

# Add cosine function with red cross points and connecting lines
#points(x, y2, col="red", pch=4)
lines(x, y2, col="red",type='o',pch=4)

# Add legend to the plot
legend("topright", inset=0.05, lty=1, legend=c("Sine", "Cos"), 
       col=c("blue", "red"), pch=c(8, 4), title="Sine/Cos Function")


#Ex3

# Define the number of programs and their frequencies
programs <- c(1, 2, 3, 4, 5, 6, 7, 8)
frequency <- c(62, 47, 39, 39, 58, 37, 4, 11)

# Convert to probability
total_families <- sum(frequency)
probability <- frequency / total_families

# Create the bar plot
bp <- barplot(probability, names.arg = programs, col = "gray90", border = "black",
              xlab = expression(italic(x) ~ "(number of assistance programs)"), 
              ylab = "Probability", 
              main = "", # No title in image
              ylim = c(0, 0.25), space = 0.5, axes=FALSE)

# Add X-axis with inward ticks
axis(1, at = bp, labels = programs, line = 0, col = "black", tck = 0.02)

# Add Y-axis with inward ticks
axis(2, at=seq(0,0.25,0.05), tck = 0.02, las=1)

# Explicitly add the X-axis line at y=0
abline(h=-0.0002, col="black", lwd=3, lty=1)


## Add grid lines
#grid(nx = NA, ny = NULL, lty = "dotted", col = "darkgray")

## Add frequency values above bars
#text(x = bp, y = frequency+2, labels = frequency, cex = 0.8)


#Ex 4
# Set up a 2x3 plotting grid
par(mfrow=c(2,3))

# (i) x vs cos(x) with red color and lines
x <- seq(-pi, pi, length.out=100)
plot(x, cos(x), type='l', col='red', lwd=2, main="Cos(x)", xlab="x", ylab="cos(x)")

# (ii) x vs (x^2 / 3) + 4.2 with violet color, points and lines, linewidth 2 and linetype 1
x <- seq(-5, 5, length.out=100)
y <- (x^2 / 3) + 4.2
plot(x, y, type='o', col='purple', lwd=2, lty=1, pch=16, main="Quadratic", xlab="x", ylab="(x^2 / 3) + 4.2")

# (iii) and (iv)

# Custom function to compute binomial probability using factorial
binomial_prob <- function(k, n, p) {
  binomial_coeff <- factorial(n) / (factorial(k) * factorial(n - k))
  return(binomial_coeff * (p^k) * ((1 - p)^(n - k)))
}

# Set parameters for first plot (p=0.3)
n <- 12
p <- 0.3
x_vals <- 0:n

# Compute probabilities
probs <- sapply(x_vals, function(k) binomial_prob(k, n, p))

# Plot histogram for Binomial distribution (p=0.3)
barplot(probs, names.arg = x_vals, col = "gray", border = "black",
        main = paste("Binomial Dist (n=", n, ", p=", p, ")"), 
        xlab = "Number of Successes", ylab = "Probability")

# Set parameters for second plot (p=0.8)
p <- 0.8
probs <- sapply(x_vals, function(k) binomial_prob(k, n, p))

# Plot histogram for Binomial distribution (p=0.8)
barplot(probs, names.arg = x_vals, col = "gray", border = "black",
        main = paste("Binomial Dist (n=", n, ", p=", p, ")"), 
        xlab = "Number of Successes", ylab = "Probability")


# (v) Histogram plot using type='h' option
x <- seq(1, 10, 0.5)
y <- 50*x / (x + 2)
colors <- rep(c("blue", "orange"), length.out=length(x))
plot(x, y, type='h', col=colors, lwd=2, main="Histogram", xlab="x", ylab="50x / (x+2)")

# (vi) x vs log(x) with orange color and ‘step’ linetype
x <- seq(1, 10, length.out=100)  # Avoid log(0)
plot(x, log(x), type='s', col='orange', lwd=2, main="Log(x) Step", xlab="x", ylab="log(x)")

par(mfrow = c(1, 1))


# Ex 5
x<-c(1,3,5,7,9,11)
y<-c(2,7,5,10,8,10)
plot(x,y,lty="dashed",col="pink",type="o",xlab = "Time",ylab = "Performance",main = "This is a graph",col.main="blue",lwd=3)

labels=c(1,3,5,7,9,11)
text(x+0.3,y,labels,col="red")
legend("topleft", inset=0.05, lty="dashed", legend="Per curve", 
       col="pink")

# Ex 6
# Custom function to compute factorial using logarithm to avoid overflow
log_factorial_manual <- function(x) {
  if (x == 0) return(0)  # log(0!) = log(1) = 0
  return(sum(log(1:x)))  # Sum of logs to compute log(x!)
}

# Custom function to compute binomial coefficient using logs
log_choose_manual <- function(n, k) {
  if (k > n || k < 0) return(-Inf)  # Handle invalid cases
  return(log_factorial_manual(n) - log_factorial_manual(k) - log_factorial_manual(n - k))
}

# Custom function to compute hypergeometric probability using logs
hypergeometric_pmf <- function(k, N, K, n) {
  if (k < 0 || k > min(K, n)) return(0)  # Invalid values return 0
  log_numerator <- log_choose_manual(K, k) + log_choose_manual(N - K, n - k)
  log_denominator <- log_choose_manual(N, n)
  return(exp(log_numerator - log_denominator))  # Convert log probability back to normal scale
}

# Set parameters
N <- 500  # Population size
K <- 50   # Number of success states
n <- 30   # Sample size

# Define range for k (number of observed successes)
k_values <- 0:min(K, n)  # k ranges from 0 to min(K, n)

# Compute probabilities using manual function
probabilities <- sapply(k_values, function(k) hypergeometric_pmf(k, N, K, n))

# Normalize probabilities to avoid floating-point issues
probabilities <- probabilities / sum(probabilities)

# Plot the bar chart for hypergeometric distribution
barplot(probabilities, names.arg = k_values, col = "blue",
        main = "Hypergeometric Distribution (N=500, K=50, n=30)",
        xlab = "Number of Successes (k)", ylab = "Probability",
        border = "black", ylim = c(0, 0.25), las = 1)

# Add grid lines for clarity
grid()


# Ex 7: Convergence of Hypergeometric to Binomial
# Custom function for log-factorial (to avoid large numbers)
log_factorial_manual <- function(x) {
  if (x == 0) return(0)  # log(0!) = log(1) = 0
  return(sum(log(1:x)))  # Sum of logs to compute log(x!)
}

binomial_prob <- function(k, n, p) {
  if (k > n || k < 0) return(0)  # Handle invalid cases
  
  # Compute binomial coefficient using logs
  binomial_coeff <- exp(log_factorial_manual(n) - log_factorial_manual(k) - log_factorial_manual(n - k))
  
  # Compute probability
  return(binomial_coeff * (p^k) * ((1 - p)^(n - k)))
}

# Set parameters
N <- 500  # Population size
K <- 50   # Number of successes in population
p <- K / N  # Approximate probability of success

# Define increasing sample sizes (9 integer values)
n_values <- round(seq(10, 250, length.out = 9))  

# Set up 3×3 grid for multiple plots
par(mfrow = c(3, 3))

for (n in n_values) {
  k_values <- 0:min(K, n)  # Possible number of successes
  
  # Compute probabilities using defined functions
  hyper_probs <- sapply(k_values, function(k) hypergeometric_pmf(k, N, K, n))
  binomial_probs <- sapply(k_values, function(k) binomial_prob(k, n, p))
  
  # Normalize probabilities for better visualization
  hyper_probs <- hyper_probs / sum(hyper_probs, na.rm = TRUE)
  binomial_probs <- binomial_probs / sum(binomial_probs, na.rm = TRUE)
  
  # Define offset to separate the two stick plots
  offset <- 0.2  # Shift amount
  
  # Plot Hypergeometric PMF (stick plot, shifted left)
  plot(k_values - offset, hyper_probs, type = "h", lwd = 2, col = "blue",
       main = paste("n =", n), xlab = "Successes (k)", ylab = "Probability",
       ylim = c(0, max(hyper_probs, binomial_probs, na.rm = TRUE)),
       xlim = c(min(k_values) - 1, max(k_values) + 1))  # Add padding
  
  # Add Binomial PMF (stick plot, shifted right)
  lines(k_values + offset, binomial_probs, type = "h", lwd = 2, col = "red")
  
  # Add legend
  legend("topleft", legend = c("Hypergeometric", "Binomial"), 
         col = c("blue", "red"), lwd = 2, inset = c(0.5,0.05),cex = 0.8)  
}

# Reset plot layout
par(mfrow = c(1, 1))

#Ex 8
poisson_pmf <- function(k, lambda) {
  return((lambda^k * exp(-lambda)) / factorial(k))
}

k_values <- 0:60  # Range of k values
lambda_values <- c(3, 20, 45)

# Compute Poisson probabilities
poisson_probs <- lapply(lambda_values, function(lambda) sapply(k_values, function(k) poisson_pmf(k, lambda)))

# Plot Poisson distributions
plot(k_values, poisson_probs[[1]], type = "h", col = "blue", lwd = 2,
     main = "Poisson Distributions (λ = 3, 20, 45)", xlab = "k", ylab = "Probability", ylim = c(0, 0.25))
lines(k_values, poisson_probs[[2]], type = "h", col = "red", lwd = 2)
lines(k_values, poisson_probs[[3]], type = "h", col = "green", lwd = 2)

# Add legend
legend("topright", legend = c("λ = 3", "λ = 20", "λ = 45"), col = c("blue", "red", "green"), lty = 1, lwd = 2)

#Ex 9
# Load dataset
data <- read.csv('/home/ibab/R/Lab9/SOCR-HeightWeight.csv')

# Print column names and dataset dimensions
print(colnames(data))
print(dim(data))

# Check structure to confirm data types
str(data)

# Convert columns to numeric if needed
height <- data$Height.Inches.
weight <- data$Weight.Pounds.

# Remove NA values if conversion introduced any
height <- height[!is.na(height)]
weight <- weight[!is.na(weight)]

# Compute mean and standard deviation for heights
mean_height <- mean(height, na.rm = TRUE)
sd_height <- sd(height, na.rm = TRUE)

# Compute mean and standard deviation for weights
mean_weight <- mean(weight, na.rm = TRUE)
sd_weight <- sd(weight, na.rm = TRUE)

cat("Mean Height:", mean_height, "\nStandard Deviation Height:", sd_height, "\n")
cat("Mean Weight:", mean_weight, "\nStandard Deviation Weight:", sd_weight, "\n")

# (i) Histogram for Heights
xname <- "Height (Inches)"

hist(height, 
     breaks = "Sturges", 
     prob = TRUE,  # Normalize to probabilities
     col = "lightgray", 
     border = "black",
     main = paste("Histogram of", xname),  
     xlim = range(height), 
     xlab = xname, 
     ylab = "Probability Density")

grid()  # Add grid for better visualization

# (ii) Histogram for Weights
xname2 <- "Weight (Pounds)"

hist(weight, 
     breaks = "Sturges", 
     prob = TRUE,  # Normalize to probabilities
     col = "lightgray", 
     border = "black",
     main = paste("Histogram of", xname2),  
     xlim = range(weight), 
     xlab = xname2, 
     ylab = "Probability Density")

grid()  # Add grid for better visualization


# (iii) 
zcalc <- function(x, mu, sd) {
  z = (x - mu) / sd
  return(z)
} 

gaussian_pdf <- function(z) {
  return((1 / sqrt(2 * pi)) * exp(-0.5 * z^2))
}

# Compute Z-scores
z_height <- zcalc(height, mean(height), sd(height))
z_weight <- zcalc(weight, mean(weight), sd(weight))

# Sort Z-scores for smooth plotting
z_height_sorted <- sort(z_height)
z_weight_sorted <- sort(z_weight)

# Compute Gaussian PDF values
pdf_height <- gaussian_pdf(z_height_sorted)
pdf_weight <- gaussian_pdf(z_weight_sorted)

# Define consistent axis limits
ylim_max <- max(pdf_height, pdf_weight)

# Plot Gaussian curve for heights
plot(z_height_sorted, pdf_height, type = "l", col = "red", lwd = 2,
     main = "Gaussian Curve (Z-transformed Heights & Weights)", 
     xlab = "Z-score", ylab = "Probability Density", ylim = c(0, ylim_max),lty=2)

# Add Gaussian curve for weights
lines(z_weight_sorted, pdf_weight, col = "blue", lwd = 2,lty=3)

# Add legend
legend("topright", legend = c("Height", "Weight"), col = c("red", "blue"), lwd = 2,lty=c(2,3))


# (iv)
# Custom function to plot histograms with varying bin sizes
plot_histogram <- function(data, bins, title) {
  hist(data, breaks = bins, col = "lightblue", main = title,
       xlab = "Value", ylab = "Frequency", border = "black")
}

# Effect of decreasing bin size in histogram
par(mfrow = c(1, 3))  # Arrange plots side by side
plot_histogram(height, 10, "Bins = 10")
plot_histogram(height, 30, "Bins = 30")
plot_histogram(height, 50, "Bins = 50")
par(mfrow = c(1, 1))  # Reset layout


# Ex 10

# Function to compute the PDF of Uniform Distribution
uniform_pdf <- function(x, a, b) {
  if (x < a || x > b) {
    return(0)  # If x is outside [a, b], probability is 0
  } else {
    return(1 / (b - a))  # Uniform density function
  }
}

# Function to compute the CPDF (Cumulative Probability Density Function)
uniform_cpdf <- function(x, a, b) {
  if (x < a) {
    return(0)  # If x is less than a, probability is 0
  } else if (x > b) {
    return(1)  # If x is greater than b, probability is 1
  } else {
    return((x - a) / (b - a))  # Cumulative probability within the range
  }
}

# Function to plot the PDF of the Uniform Distribution with vertical shading
plot_uniform_pdf <- function(a, b, shade_up_to) {
  x_vals <- seq(a - 1, b + 1, length.out = 100)
  pdf_vals <- sapply(x_vals, function(x) uniform_pdf(x, a, b))
  
  plot(x_vals, pdf_vals, type = "l", col = "blue", lwd = 2,
       main = paste("Uniform PDF: U(", a, ",", b, ")"),
       xlab = "x", ylab = "Density", ylim = c(0, 1.2))
  
  # Shade the region under the curve using vertical line segments
  shade_x_vals <- seq(a, shade_up_to, length.out = 100)
  for (x in shade_x_vals) {
    lines(c(x, x), c(0, uniform_pdf(x, a, b)), col = rgb(0, 0, 1, 0.3), lwd = 1)  # Vertical shading
  }
}

# Function to plot the CPDF of the Uniform Distribution with vertical shading
plot_uniform_cpdf <- function(a, b, shade_up_to) {
  x_vals <- seq(a - 1, b + 1, length.out = 100)
  cpdf_vals <- sapply(x_vals, function(x) uniform_cpdf(x, a, b))
  
  plot(x_vals, cpdf_vals, type = "l", col = "red", lwd = 2,
       main = paste("Uniform CPDF: U(", a, ",", b, ")"),
       xlab = "x", ylab = "Cumulative Probability")
  
  # Shade the region under the curve using vertical line segments
  shade_x_vals <- seq(a, shade_up_to, length.out = 100)
  for (x in shade_x_vals) {
    lines(c(x, x), c(0, uniform_cpdf(x, a, b)), col = rgb(1, 0, 0, 0.3), lwd = 1)  # Vertical shading
  }
}

# Execute plotting functions
par(mfrow = c(1, 2))  # Arrange two plots side by side
plot_uniform_pdf(1, 2, 1.5)  # Shade region under PDF up to x = 1.5
plot_uniform_cpdf(1, 2, 1.5)  # Shade region under CPDF up to x = 1.5



# Ex 11,12,13

exponential_pdf <- function(x, lambda) {
  result <- numeric(length(x))
  for (i in seq_along(x)) {
    if (x[i] < 0) {
      result[i] <- 0
    } else {
      result[i] <- lambda * exp(-lambda * x[i])
    }
  }
  return(result)
}

exponential_cpdf <- function(x, lambda) {
  result <- numeric(length(x))
  for (i in seq_along(x)) {
    if (x[i] < 0) {
      result[i] <- 0
    } else {
      result[i] <- 1 - exp(-lambda * x[i])
    }
  }
  return(result)
}

gamma_pdf <- function(x, alpha, theta) {
  result <- numeric(length(x))
  for (i in seq_along(x)) {
    if (x[i] < 0) {
      result[i] <- 0
    } else {
      result[i] <- (x[i]^(alpha - 1) * exp(-x[i] / theta)) / (gamma(alpha) * theta^alpha)
    }
  }
  return(result)
}

gamma_cpdf <- function(x, alpha, theta) {
  return(pgamma(x, shape = alpha, scale = theta))
}

chisq_pdf <- function(x, df) {
  result <- numeric(length(x))
  for (i in seq_along(x)) {
    if (x[i] < 0) {
      result[i] <- 0
    } else {
      result[i] <- (x[i]^(df/2 - 1) * exp(-x[i]/2)) / (2^(df/2) * gamma(df/2))
    }
  }
  return(result)
}

chisq_cpdf <- function(x, df) {
  return(pchisq(x, df))
}

plot_exponential <- function(lambda, shade_up_to) {
  x_vals <- seq(0, 5, length.out = 100)
  pdf_vals <- sapply(x_vals, function(x) exponential_pdf(x, lambda))
  cpdf_vals <- sapply(x_vals, function(x) exponential_cpdf(x, lambda))
  
  plot(x_vals, pdf_vals, type = "l", col = "blue", lwd = 2,
       main = paste("Exponential PDF (λ =", lambda, ")"),
       xlab = "x", ylab = "Density")
  
  shade_x_vals <- seq(0, shade_up_to, length.out = 100)
  for (x in shade_x_vals) {
    y_val <- exponential_pdf(x, lambda)  
    if (y_val > 0) {  
      segments(x0 = x, y0 = 0, x1 = x, y1 = y_val, 
               col = rgb(0, 0, 1, 0.6), lwd = 1)
    }
  }
  
  plot(x_vals, cpdf_vals, type = "l", col = "red", lwd = 2,
       main = paste("Exponential CPDF (λ =", lambda, ")"),
       xlab = "x", ylab = "Cumulative Probability")
}

plot_gamma <- function(alpha, theta, shade_up_to) {
  x_vals <- seq(0, 40, length.out = 100)
  pdf_vals <- gamma_pdf(x_vals, alpha, theta)
  cpdf_vals <- gamma_cpdf(x_vals, alpha, theta)
  
  plot(x_vals, pdf_vals, type = "l", col = "blue", lwd = 2,
       main = paste("Gamma PDF (α =", alpha, ", θ =", theta, ")"),
       xlab = "x", ylab = "Density")
  
  shade_x_vals <- seq(0, shade_up_to, length.out = 50)
  for (x in shade_x_vals) {
    y_val <- gamma_pdf(x, alpha, theta)  
    if (y_val > 0) {  
      segments(x0 = x, y0 = 0, x1 = x, y1 = y_val, 
               col = rgb(0, 0, 1, 0.6), lwd = 1)
    }
  }
  
  plot(x_vals, cpdf_vals, type = "l", col = "red", lwd = 2,
       main = paste("Gamma CPDF (α =", alpha, ", θ =", theta, ")"),
       xlab = "x", ylab = "Cumulative Probability")
}

plot_chisq <- function(df, shade_up_to) {
  x_vals <- seq(0, 50, length.out = 100)
  pdf_vals <- sapply(x_vals, function(x) chisq_pdf(x, df))
  cpdf_vals <- sapply(x_vals, function(x) chisq_cpdf(x, df))
  
  plot(x_vals, pdf_vals, type = "l", col = "blue", lwd = 2,
       main = paste("Chi-square PDF (df =", df, ")"),
       xlab = "x", ylab = "Density")
  
  shade_x_vals <- seq(0, shade_up_to, length.out = 200)
  for (x in shade_x_vals) {
    y_val <- chisq_pdf(x, df)  
    if (y_val > 0) {
      segments(x0 = x, y0 = 0, x1 = x, y1 = y_val, 
               col = rgb(0, 0, 1, 0.6), lwd = 1)
    }
  }
  
  plot(x_vals, cpdf_vals, type = "l", col = "red", lwd = 2,
       main = paste("Chi-square CPDF (df =", df, ")"),
       xlab = "x", ylab = "Cumulative Probability")
}

# Set up 3x2 grid layout
par(mfrow = c(3, 2))

# Execute functions with shading applied using segments()
plot_exponential(10, 2.8)  # λ = 10, shade up to x = 2.8
plot_gamma(5, 3, 10)       # α = 5, θ = 3, shade up to x = 10
plot_chisq(20, 1.0)        # df = 20, shade up to x = 1.0

# Reset layout
par(mfrow = c(1, 1))
