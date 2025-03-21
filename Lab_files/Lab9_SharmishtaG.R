# NOTES

# Plots
# barplot()	Bar graph	Categorical data with bars with heights/lengths proportional to the values
# boxplot()	Box-and-whisker plot	Shows distribution of data points through quartiles
# dotchart()	Cleveland dotplot	An alternative to a bar graph; plots a dot for each observation on a scale
# hist()	Histogram	Distribution of numerical data
# plot()	Scatterplot	Plotted dots represent values for two numeric variables on the horizontal and vertical axes
# pie()	Pie chart	Displays data as a percentage of a whole

# Redirecting Output to a File in R
# svg("mygraph.svg")	Recommended: Resize without fuzziness or pixelation
# pdf("mygraph.pdf")	Recommended: Resize without fuzziness or pixelation
# win.metafile("mygraph.wmf")	Ideal for graphs used in Microsoft Word or PowerPoint
# png("mygraph.png")	Recommended for webpages
# postscript("mygraph.eps")	Suitable for embedding an image in documents or sending to printers
# tiff("mygraph.tiff")	Large file size
# jpeg("mygraph.jpg")	Generally for photographs/pixelated images
# bmp("mygraph.bmp")	Seldom used outside of the Windows operating system

# pch=0	Square
# pch=1	Circle
# pch=2	Triangle (point up)
# pch=3	Plus (+)
# pch=4	Cross (x)
# pch=5	Diamond
# pch=6	Triangle (point down)
# pch=7	Square cross
# pch=8	Star
# pch=9	Diamond plus
# pch=10	Circle plus
# pch=11	Triangles up and down
# pch=12	Square plus
# pch=13	Circle cross
# pch=14	Square and triangle down
# pch=15	Filled square (blue)
# pch=16	Filled circle (blue)
# pch=17	Filled triangle (point up, blue)
# pch=18	Filled diamond (blue)
# pch=19	Solid circle (blue)
# pch=20	Bullet (smaller circle)
# pch=21	Filled circle (red)
# pch=22	Filled square (red)
# pch=23	Filled diamond (red)
# pch=24	Filled triangle (point up, red)
# pch=25	Filled triangle (point down, red)

# Understanding cex in R
# cex = 1 → Default size
# cex = 1.5 → 150% of the default size
# cex = 0.5 → 50% of the default size
# Special Variants of cex
# Parameter	Effect
# cex.axis	Scales the axis labels
# cex.lab	Scales the axis titles (labels)
# cex.main	Scales the main plot title
# cex.sub	Scales the subtitle

# type="p"	Plots only points
# type="l"	Plots only lines
# type="b"	Plots both points and lines
# type="o"	Plots points overlaid by lines
# type="h"	Plots vertical lines like a histogram
# type="s"	Plots step-wise (staircase) lines
# type="n"	No plotting, just marks axes

#col="blue"

# lty = 0	"blank" (no line)
# lty = 1	"solid" (default)
# lty = 2	"dashed"
# lty = 3	"dotted"
# lty = 4	"dotdash" (dot-dash pattern)
# lty = 5	"longdash" (long dashes)
# lty = 6	"twodash" (two short dashes)

# col.main	Sets the color of the main title (uses same values as col).
# font.main	Sets the font style of the main title:
# 1 → Plain
# 2 → Bold
# 3 → Italic
# 4 → Bold Italic
# cex.main	Scales the size of the main title.

# Basic Plot with Title Customization
plot(xval, yval, pch=18, cex=1.8, col="purple", type="o", lty=6, lwd=1,
     main="This is the main title", col.main="blue", font.main=2, cex.main=1.5)

# Adding Axis Labels
plot(xval, yval, pch=18, cex=1.8, col="purple", type="o", lty=6, lwd=1,
     main="This is the main title", col.main="blue", font.main=1, cex.main=1.5,
     xlab = "Concentration (mmol/L)", ylab = "Velocity (mmol/L/sec)",
     font.lab = 2, col.lab="brown", cex.lab = 1.2)

# Adjusting Plot Ranges
# Range is set using xlim and ylim parameters
plot(Xvalue, Yvalue, col="blue", type="o", xlim=c(1,20), ylim=c(1,150))



# Text inside the plot

# Create data
Xvalue = c(1,2,3,4,5,6,7,8,9,10)
Yvalue = c(12, 23, 36, 48, 53, 64, 78, 89, 91, 110)

# Create a basic plot with circles as points
plot(Xvalue, Yvalue, type="o")

# Add text inside the plot at a specific coordinate
# The text will be centered at position (3, 100)
plot(Xvalue, Yvalue, text(3,100,"This is text inside plot", col="brown", font=2))



# Adding Labels to Data Points

Xvalue = c(1,2,3,4,5,6,7,8,9,10)
Yvalue = c(12, 23, 36, 48, 53, 64, 78, 89, 91, 110)

# These letters are labels for data points
cch = c("a","b","c","d","e","f","g","h","i","j")

# Plotting the graph with circles
plot(Xvalue, Yvalue, type="o")

# Adding text labels to each point
# The labels are placed slightly to the right of each point (Xvalue+0.3)
text(Xvalue+0.3, Yvalue, cch, col="blue", font=2)

# Adding a second text string inside the plot
text(5, 90, "This is text", col="red")


# Using Logarithmic Scales

# Create your data (Xval and Yval should be defined before this code)
plot(Xval, Yval, log="y", type="o", col="blue", lwd=2)



# Overlaying Multiple Graphs

x <- seq(2, 10, 0.1)
y1 <- x^2
y2 <- x^3

# Plot the first function (squares) with lines, color red
plot(x, y1, type='l', col='red')

# Add the second function (cubes) with lines, color green
lines(x, y2, col='green')

# Add a legend
legend('bottomright', inset=0.05, c("Squares", "Cubes"), lty=1, 
       col=c('red', 'dark green'), title="Graph function")



# Adding Horizontal and Vertical Lines

x <- seq(2, 10, 0.1)
y1 <- x^2
y2 <- x^3

# Plot the first function (squares)
plot(x, y1, type='l', col='red')

# Add the second function (cubes)
lines(x, y2, col='green')

# Add legend
legend('bottomright', inset=0.05, c("Squares", "Cubes"), lty=1, 
       col=c('red', 'dark green'), title="Graph function")

# Add a horizontal line at y=4 and vertical line at x=5
abline(a=4, b=0, col='blue')  # Horizontal line: y = a + b*x (where b=0 makes it horizontal)

# Add vertical lines at multiple x positions
abline(h=c(4, 6, 8), col="dark green", lty=2)  # Horizontal lines at y=4, y=6, and y=8

# Add horizontal lines at multiple y positions  
abline(v=c(4, 6, 8), col="dark green", lty=1)  # Vertical lines at x=4, x=6, and x=8



-------------------------------------------------------------------------------------------------
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
axis(2, at=seq(0,0.25,0.05), tck = 0.02)

# Explicitly add the X-axis line at y=0
abline(h=0.000, col="black", lwd=3, lty=1)


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
         col = c("blue", "red"), lwd = 2, inset = 0.05,cex = 0.8)  
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
data <- read.csv('/Users/sharmishtaganesh/Desktop/Biostat_Lab/Lab_files/SOCR-HeightWeight.csv')

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

par(mfrow=c(1,2))
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
par(mfrow=c(1,1))

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



------------------------------------WORKSHEET------------------------------------------------------
  
#1)
  
# (a) Read the data file "expression.csv" into a data frame and print column names
  expression_data <- read.csv("expression.csv", header = TRUE)
print(colnames(expression_data))

# (b) Compute control mean, treatment mean, and RatioCol
control_columns <- expression_data[, 2:9]  # Controls assumed in columns 2-9
treatment_columns <- expression_data[, 10:17]  # Treatments assumed in columns 10-17

control_mean <- rowMeans(control_columns, na.rm = TRUE)
treatment_mean <- rowMeans(treatment_columns, na.rm = TRUE)

# Avoid division by zero
RatioCol <- ifelse(treatment_mean != 0, control_mean / treatment_mean, NA)
expression_data$RatioCol <- RatioCol  # Add column to dataframe

# Print first 10 values of RatioCol
print(RatioCol[1:10])

# (c) Create 2x2 plot layout and generate the required plots
par(mfrow = c(2, 2))  # Divide canvas into 2x2 layout

# (i) Histogram of first control column (e.g., Control1)
hist(expression_data[, 2], main = "Histogram of Control1", col = "lightblue", 
     xlab = "Expression Levels", breaks = 30)

# (ii) Histogram of a treatment column (e.g., Treatment3)
hist(expression_data[, 12], main = "Histogram of Treatment3", col = "lightgreen", 
     xlab = "Expression Levels", breaks = 30)

# (iii) Scatter plot between Control1 and Control5
plot(expression_data[, 2], expression_data[, 6], main = "Control1 vs Control5",
     xlab = "Control1", ylab = "Control5", col = "blue", pch = 19)

# (iv) Scatter plot between Control6 and Treatment6
plot(expression_data[, 7], expression_data[, 15], main = "Control6 vs Treatment6",
     xlab = "Control6", ylab = "Treatment6", col = "red", pch = 19)

# Reset plotting layout
par(mfrow = c(1, 1))


# 2)
# (i) Load the data
expression_data <- read.csv("expression.csv", header = TRUE)

# (ii) Create subset 'subi' where both control1 and control2 values are > 50
subi <- subset(expression_data, expression_data[, 2] > 50 & expression_data[, 3] > 50)

# (iii) Create subset 'subrow' for specific genes
gene_names <- c("gene-4", "gene-20", "gene-37", "gene-100")  # Gene names to extract
subrow <- expression_data[expression_data[, 1] %in% gene_names, ]  # Assuming first column contains gene names

# (iv) Randomly sample 200 rows from the expression data
set.seed(123)  # Ensuring reproducibility
samp <- expression_data[sample(nrow(expression_data), 200), ]

# (v) Create a boxplot comparing selected control and treatment columns
selected_cols <- c(5, 6, 7, 10, 14, 17)  # Assuming these correspond to control4, control5, control6, treatment1, treatment5, treatment8

boxplot(expression_data[, selected_cols], 
        main = "Boxplot of Selected Controls and Treatments",
        names = colnames(expression_data)[selected_cols],  # Add proper labels
        col = c("lightblue", "lightblue", "lightblue", "lightgreen", "lightgreen", "lightgreen"),
        ylab = "Expression Levels")

# 3)
# (i) Generate 10,000 random deviates from a Gaussian distribution with mean=20 and sd=4
set.seed(123)  # Ensuring reproducibility
data_gaussian <- rnorm(10000, mean = 20, sd = 4)
x <- seq(min(data_gaussian), max(data_gaussian), length.out = 100)
# Plot histogram
hist(data_gaussian, main = "Histogram of Gaussian Distribution (mean=20, sd=4)",
     xlab = "Values", col = "lightblue", border = "black", freq = FALSE)

# Overlay normal density curve
curve(dnorm(x, mean = 20, sd = 4), col = "red", lwd = 2, add = TRUE)


# (ii) Compute the Full Width at Half Maximum (FWHM)
# Formula: FWHM = 2.355 * standard deviation (for a normal distribution)
FWHM <- 2.355 * 4
print(paste("Full Width at Half Maximum (FWHM):", FWHM))


# (iii) Compute the area under the curve for Gaussian P(x, μ=30, σ=6) between x=35 and x=50
mu <- 30  # Mean
sigma <- 6  # Standard deviation
x1 <- 35
x2 <- 50

# Compute cumulative probabilities
p1 <- pnorm(x1, mean = mu, sd = sigma)
p2 <- pnorm(x2, mean = mu, sd = sigma)

# Area under the curve between x1 and x2
area_under_curve <- p2 - p1
print(paste("Area under the curve between 35 and 50:", area_under_curve))


# 4) i)

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
probs1 <- sapply(x_vals, function(k) binomial_prob(k, n, p))

# Set parameters for second plot (p=0.8)
n<-13
p <- 0.7
probs2 <- sapply(x_vals, function(k) binomial_prob(k, n, p))

# Plot the first binomial distribution (n=12, p=0.3) in blue
plot(x_vals, probs1, type = "b", col = "blue", pch = 16, lwd = 2,
     xlab = "x", ylab = "Probability", main = "Binomial Distributions",
     ylim = c(0, max(c(probs1, probs2))))

# Add the second binomial distribution (n=13, p=0.7) in red
lines(x_vals, probs2, type = "b", col = "red", pch = 16, lwd = 2)

# Add a legend
legend("topright", legend = c("n=12, p=0.3", "n=13, p=0.7"), col = c("blue", "red"), pch = 16, lwd = 2)


# 4) ii)

# Set seed for reproducibility
set.seed(123)

# (i) Generate 5000 random deviates from a Poisson distribution with λ = 8
poisson_data <- rpois(5000, lambda = 8)

# Compute mean and standard deviation of Poisson distribution
mean_poisson <- mean(poisson_data)
sd_poisson <- sd(poisson_data)

# (ii) Plot Poisson frequency histogram
hist(poisson_data, breaks = 30, probability = TRUE, col = "lightblue", border = "black",
     main = "Poisson Distribution with Gaussian Approximation",
     xlab = "Values", ylim = c(0, 0.2))

# Define x values for the Gaussian curve
x_vals <- seq(min(poisson_data), max(poisson_data), length.out = 100)

# Compute Gaussian probabilities
gaussian_probs <- dnorm(x_vals, mean = mean_poisson, sd = sd_poisson)

# Overlay Gaussian curve in red
lines(x_vals, gaussian_probs, col = "red", lwd = 2)

# Add a legend
legend("topright", legend = c("Poisson (λ=8)", "Gaussian Approximation"), col = c("black", "red"), lwd = 2)


#[alternative]

# Compute mean and standard deviation of the Poisson distribution
# Set seed for reproducibility
set.seed(123)

# (i) Generate 5000 random deviates from a Poisson distribution with λ = 8
poisson_data <- rpois(5000, lambda = 8)

# Custom function to compute Poisson PMF
poisson_pmf <- function(k, lambda) {
  return((lambda^k * exp(-lambda)) / factorial(k))
}

# Custom function to compute Z-score
zcalc <- function(x, mu, sd) {
  z = (x - mu) / sd
  return(z)
} 

# Custom function to compute Gaussian PDF
gaussian_pdf <- function(z) {
  return((1 / sqrt(2 * pi)) * exp(-0.5 * z^2))
}

# Compute mean and standard deviation of Poisson distribution
mean_poisson <- mean(poisson_data)
sd_poisson <- sd(poisson_data)

# Define a range of k values for Poisson probabilities
k_values <- 0:20  # Reasonable range based on lambda = 8

# Compute Poisson probabilities
poisson_probs <- sapply(k_values, function(k) poisson_pmf(k, lambda = 8))

# Compute Gaussian approximation using the same mean and standard deviation
z_values <- sapply(k_values, function(k) zcalc(k, mean_poisson, sd_poisson))
gaussian_probs <- sapply(z_values, gaussian_pdf)

# Scale Gaussian to match Poisson distribution's height
gaussian_probs <- gaussian_probs * max(poisson_probs) / max(gaussian_probs)

# (ii) Plot Poisson histogram
hist(poisson_data, breaks = 30, probability = TRUE, col = "lightblue", border = "black",
     main = "Poisson Distribution with Gaussian Approximation",
     xlab = "Values", ylim = c(0, max(poisson_probs) * 1.2))

# Overlay Poisson PMF as blue points
points(k_values, poisson_probs, col = "blue", pch = 19, cex = 1.2)

# Overlay Gaussian approximation in red
lines(k_values, gaussian_probs, col = "red", lwd = 2)

# Add legend
legend("topright", legend = c("Poisson (λ=8)", "Gaussian Approximation"), col = c("blue", "red"), lwd = 2, pch = c(19, NA))


# 3) 
# (i) [alternative] - custom
# Set seed for reproducibility
set.seed(123)

# (i) Generate 10,000 random deviates using Box-Muller Transform
generate_gaussian_data_manual <- function(n, mean, sd) {
  u1 <- runif(n)  # Generate uniform random numbers
  u2 <- runif(n)
  
  # Box-Muller transform to generate standard normal deviates
  z <- sqrt(-2 * log(u1)) * cos(2 * pi * u2)
  
  # Scale to desired mean and standard deviation
  data <- mean + sd * z
  return(data)
}

# Generate 10,000 Gaussian deviates with mean=20, sd=4
data_gaussian <- generate_gaussian_data_manual(10000, mean = 20, sd = 4)

# Compute Z-scores
compute_z_score <- function(x, mean, sd) {
  return((x - mean) / sd)
}

z_scores <- sapply(data_gaussian, function(x) compute_z_score(x, mean = 20, sd = 4))


# (ii) Compute Gaussian PDF Using Formula
gaussian_pdf_manual <- function(x, mean, sd) {
  coef <- 1 / (sd * sqrt(2 * pi))  # Prefactor
  exponent <- exp(-((x - mean)^2) / (2 * sd^2))  # Exponential part
  return(coef * exponent)
}

# (iii) Compute Probabilities Using Trapezoidal Integration for Gaussian CDF
compute_gaussian_cdf_manual <- function(x, mean, sd, steps = 1000) {
  x_vals <- seq(-10 * sd + mean, x, length.out = steps)  # Approximate range for integration
  y_vals <- sapply(x_vals, function(x) gaussian_pdf_manual(x, mean, sd))
  
  # Trapezoidal approximation of the integral
  area <- sum((y_vals[-1] + y_vals[-length(y_vals)]) / 2 * (x_vals[2] - x_vals[1]))
  return(area)
}

# Compute Probability P(X ≤ 35) for N(30,6)
probability_35 <- compute_gaussian_cdf_manual(35, mean = 30, sd = 6)
print(paste("P(X ≤ 35) for N(30,6):", probability_35))

# Compute Probability P(35 ≤ X ≤ 50) by subtraction
probability_50 <- compute_gaussian_cdf_manual(50, mean = 30, sd = 6)
area_under_curve <- probability_50 - probability_35
print(paste("Area under the curve between 35 and 50:", area_under_curve))


# (iv) Plot Histogram and Gaussian Curve
hist(data_gaussian, breaks = 30, probability = TRUE, col = "lightblue", border = "black",
     main = "Gaussian Distribution (Manual Calculation)", xlab = "Values")

# Overlay Gaussian PDF
x_vals <- seq(min(data_gaussian), max(data_gaussian), length.out = 100)
pdf_vals <- sapply(x_vals, function(x) gaussian_pdf_manual(x, mean = 20, sd = 4))
lines(x_vals, pdf_vals, col = "red", lwd = 2)

# Add legend
legend("topright", legend = c("Histogram", "Gaussian PDF"), col = c("black", "red"), lwd = 2)


