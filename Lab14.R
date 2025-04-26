# EX I

# Question 1
means <- c(20.34, 19.49, 25.68)
stderr <- c(0.83, 1.51, 1.39)

# Create the barplot and store midpoints of bars
bar_positions <- barplot(
  means,
  names.arg = c("A", "B", "C"),
  col = "grey",
  ylim = c(0, max(means + stderr) + 2),
  main = "Errors on bar plot"
)

# Add error bars using arrows()
arrows(
  x0 = bar_positions, y0 = means + stderr,
  x1 = bar_positions, y1 = means - stderr,
  angle = 90, code = 3, length = 0.06, col = "red"
)
# This line draws vertical error bars centered on each bar, extending from the 
# mean + stderr to mean - stderr, with horizontal caps at both ends, all in red. 
# Perfect for showing uncertainty in measurements!


# Question 2
x = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)
y = c(5, 9, 13, 17, 20, 24, 26, 27, 28, 27)
errors = c(0.5, 0.9, 1.4, 1.5, 2.0, 2.2, 2.3, 2.5, 2.9, 3.0)

# Create the barplot and store midpoints of bars
plot(x,y,pch=16,xlab ="concentration",ylab = "optical activity",main = "Error bars on data points",ylim = c(min(y - errors) - 1, max(y + errors) + 1))

# Add vertical error bars using arrows()
arrows(
  x0 = x, y0 = y + errors,
  x1 = x, y1 = y - errors,
  angle = 90, code = 3, length = 0.05, col = "blue"
)


# Question 3
x = c(10,20,30,40,50,60,70,80,90,100)
y = c(95, 220, 279, 424, 499, 540, 720, 880, 950, 1200)
cov(x,y)
cor(x,y)
cor(longley) #multivariate data


# EX II - One sample T test

# Question 1

# Z-test with Q-Q plot for normality check
one_sample_Ztest_pvalue <- function(x, muzero, alpha, null) {
  n <- length(x)
  x_bar <- mean(x)
  sigma <- sd(x)
  
  cat("Sample size (n):", n, "\n")
  cat("Sample mean (x̄):", round(x_bar, 3), "\n")
  cat("Population standard deviation (σ):", sigma, "\n")
  cat("Hypothesized population mean (μ₀):", muzero, "\n")
  
  # Q-Q plot for normality check
  qqnorm(x)
  qqline(x, col = "red")
  cat("Q-Q plot generated. Check if the points lie along the red line to assess normality.\n")
  
  # Calculate Z statistic
  Z <- (x_bar - muzero) / (sigma / sqrt(n))
  cat("Computed Z-statistic:", round(Z, 4), "\n")
  
  # Critical Z value calculation for different test types
  if (null == "equal") {
    critical_value <- qnorm(1 - alpha / 2)  # Two-tailed critical value
    cat("Test Type: Two-tailed (H₀: μ = μ₀ vs H₁: μ ≠ μ₀)\n")
    
  } else if (null == "less_than_or_equal") {
    critical_value <- qnorm(1 - alpha)  # Right-tailed test critical value
    cat("Test Type: Right-tailed (H₀: μ ≤ μ₀ vs H₁: μ > μ₀)\n")
    
  } else if (null == "more_than_or_equal") {
    critical_value <- qnorm(alpha)  # Left-tailed test critical value
    cat("Test Type: Left-tailed (H₀: μ ≥ μ₀ vs H₁: μ < μ₀)\n")
    
  } else {
    stop("Invalid null hypothesis. Use 'equal', 'less_than_or_equal', or 'more_than_or_equal'.")
  }
  
  # Decision based on p-value
  if (null == "equal") {
    p_value <- 2 * (1 - pnorm(abs(Z)))  # Two-tailed test p-value
  } else if (null == "less_than_or_equal") {
    p_value <- 1 - pnorm(Z)  # Right-tailed test p-value
  } else if (null == "more_than_or_equal") {
    p_value <- pnorm(Z)  # Left-tailed test p-value
  }
  
  if (p_value < alpha) {
    conclusion <- "Reject H₀ (p-value < α or α/2(for 2 tailed) ⇒ significant)"
  } else {
    conclusion <- "Fail to reject H₀ (p-value ≥ α or α/2(for 2 tailed) ⇒ not significant)"
  }
  
  # Output the results
  cat("Computed p-value:", round(p_value, 5), "\n")
  cat("Significance level (α):", alpha, "\n")
  cat("Critical value for Z (", null, "test):", round(critical_value, 4), "\n")
  cat("Conclusion:", conclusion, "\n")
  
  return(list(Z_value = Z, p_value = p_value, conclusion = conclusion, critical_value = critical_value))
}

# Dataset
x <- c(96.0, 104.0, 99.1, 97.6, 99.4, 92.8, 105.6, 97.2,
       96.8, 92.1, 100.6, 101.5, 100.7, 97.3, 99.6, 105.9)

# Parameters
mu0 <- 100
alpha <- 0.05
null <- "equal"  # Two-tailed test

# Run the test
result <- one_sample_Ztest_pvalue(x, mu0, alpha, null)


# Question 2

# Function to perform a one-sample t-test
one_sample_t_test <- function(x, muzero, alpha, null) {
  n <- length(x)  # Sample size
  df <- n - 1  # Degrees of freedom
  x_bar <- mean(x)  # Sample mean
  s <- sd(x)  # Sample standard deviation
  
  # Print sample statistics
  cat("Sample size (n):", n, "\n")
  cat("Sample mean (x̄):", round(x_bar, 3), "\n")
  cat("Sample standard deviation (s):", round(s, 3), "\n")
  cat("Hypothesized mean (μ₀):", muzero, "\n")
  
  # Calculate the t-statistic
  t_stat <- (x_bar - muzero) / (s / sqrt(n))
  cat("Computed t-statistic:", round(t_stat, 4), "\n")
  
  # Critical t value calculation for different test types
  if (null == "equal") {
    critical_value <- qt(1 - alpha / 2, df)  # Two-tailed test critical value
    cat("Test Type: Two-tailed (H₀: μ = μ₀ vs H₁: μ ≠ μ₀)\n")
    
  } else if (null == "less_than_or_equal") {
    critical_value <- qt(1 - alpha, df)  # Right-tailed test critical value
    cat("Test Type: Right-tailed (H₀: μ ≤ μ₀ vs H₁: μ > μ₀)\n")
    
  } else if (null == "more_than_or_equal") {
    critical_value <- qt(alpha, df)  # Left-tailed test critical value
    cat("Test Type: Left-tailed (H₀: μ ≥ μ₀ vs H₁: μ < μ₀)\n")
    
  } else {
    stop("Invalid null hypothesis. Use 'equal', 'less_than_or_equal', or 'more_than_or_equal'.")
  }
  
  # Decision based on p-value
  if (null == "equal") {
    p_value <- 2 * (1 - pt(abs(t_stat), df))  # Two-tailed test p-value
  } else if (null == "less_than_or_equal") {
    p_value <- 1 - pt(t_stat, df)  # Right-tailed test p-value
  } else if (null == "more_than_or_equal") {
    p_value <- pt(t_stat, df)  # Left-tailed test p-value
  }
  
  # Conclusion based on the p-value and alpha
  if (p_value < alpha) {
    conclusion <- "Reject H₀ (p-value < α or α/2(for 2 tailed) ⇒ significant)"
  } else {
    conclusion <- "Fail to reject H₀ (p-value ≥ α or α/2(for 2 tailed) ⇒ not significant)"
  }
  
  # Output the results
  cat("Computed p-value:", round(p_value, 5), "\n")
  cat("Significance level (α):", alpha, "\n")
  cat("Critical value for t (", null, "test):", round(critical_value, 4), "\n")
  cat("Conclusion:", conclusion, "\n")
  
  return(list(t_value = t_stat, p_value = p_value, conclusion = conclusion, critical_value = critical_value))
}

# Dataset
x <- c(96.0, 104.0, 99.1, 97.6, 99.4, 92.8, 105.6, 97.2,
       96.8, 92.1, 100.6, 101.5, 100.7, 97.3, 99.6, 105.9)

# Null hypothesis parameters
muzero <- 100
alpha <- 0.05
null <- "equal"  # Two-tailed test

# Perform the one-sample t-test
result <- one_sample_t_test(x, muzero, alpha, null)


# Question 3

# Given data
x <- 710        # Number of successes (e.g., users clicked on an ad)
n <- 2600       # Total number of trials (e.g., users shown the ad)
p0 <- 0.25      # Hypothesized population proportion
alpha <- 0.05   # Significance level
alt <- "greater"  # We want to test if true proportion > 0.25

cat("\nExact Binomial Test\n")
binom_result <- binom.test(x = x, n = n, p = p0, alternative = alt)

# Print the results
print(binom_result)

# Explanation:
# binom.test uses the exact binomial distribution.
# It's most accurate for small or moderate n, but can still be used for larger n.

# 2. Approximate test using prop.test() WITH continuity correction
cat("\nProportion Test WITH Yates' Continuity Correction\n")
prop_result_correct <- prop.test(x = x, n = n, p = p0, alternative = alt, correct = TRUE)

# Print the results
print(prop_result_correct)

# Explanation:
# prop.test() uses normal approximation (z-test for proportions).
# With correct = TRUE, it applies Yates' continuity correction.
# This correction makes the test more conservative, especially when n is small(30-40).

# 3. Proportion test WITHOUT continuity correction
cat("\nProportion Test WITHOUT Continuity Correction\n")
prop_result_nocorrect <- prop.test(x = x, n = n, p = p0, alternative = alt, correct = FALSE)

# Print the results
print(prop_result_nocorrect)

# Explanation:
# When correct = FALSE, it skips Yates’ correction.
# Useful when n is large (like here: n = 2600), as correction may be overly conservative.

# Inference and Comparison

cat("\n Inference \n")
cat("P-value from binom.test:", binom_result$p.value, "\n")
cat("P-value from prop.test with correction:", prop_result_correct$p.value, "\n")
cat("P-value from prop.test without correction:", prop_result_nocorrect$p.value, "\n")
# No much difference in p value with/without correction because n is large (sample size)

# Confidence Intervals from intrinsic functions
cat("\nConfidence Interval (binom.test): ", round(binom_result$conf.int[1], 4), "to", round(binom_result$conf.int[2], 4), "\n")
cat("Confidence Interval (prop.test with correction): ", round(prop_result_correct$conf.int[1], 4), "to", round(prop_result_correct$conf.int[2], 4), "\n")
cat("Confidence Interval (prop.test without correction): ", round(prop_result_nocorrect$conf.int[1], 4), "to", round(prop_result_nocorrect$conf.int[2], 4), "\n")

# Decision rule
if (binom_result$p.value < alpha) {
  cat("\nConclusion from binom.test: Reject H₀ — evidence suggests true p > 0.25\n")
} else {
  cat("\nConclusion from binom.test: Fail to reject H₀ — insufficient evidence\n")
}

if (prop_result_correct$p.value < alpha) {
  cat("Conclusion from prop.test (corrected): Reject H₀ — evidence suggests true p > 0.25\n")
} else {
  cat("Conclusion from prop.test (corrected): Fail to reject H₀ — not significant\n")
}

if (prop_result_nocorrect$p.value < alpha) {
  cat("Conclusion from prop.test (uncorrected): Reject H₀ — stronger evidence (less conservative)\n")
} else {
  cat("Conclusion from prop.test (uncorrected): Fail to reject H₀ — not significant\n")
}

# Summary of Conceptual Differences-
# cat("\nConcept Summary\n")
# cat("• binom.test is based on exact binomial distribution (no approximation).\n")
# cat("• prop.test is a large-sample z-test for proportions using normal approximation.\n")
# cat("• Yates’ continuity correction in prop.test adjusts for discreteness of binomial vs continuous normal (used when n is small).\n")
# cat("• In large n scenarios, correction can make the test overly conservative.\n")


# Question 4 

one_sample_variance_test <- function(x, test_sigma, alpha, tail_type = "two-tailed") {
  # Step 1: Calculate sample statistics
  n <- length(x)  # Sample size
  sample_variance <- var(x)  # Sample variance
  
  # Step 2: Compute the chi-squared test statistic
  chi_squared_stat <- (n - 1) * sample_variance / (test_sigma^2)
  
  # Step 3: Determine critical values based on the tail type
  if (tail_type == "two-tailed") {
    # Two-tailed test: Calculate both the lower and upper critical values
    lower_critical_value <- qchisq(alpha / 2, df = n - 1)  # Lower bound (left tail)
    upper_critical_value <- qchisq(1 - alpha / 2, df = n - 1)  # Upper bound (right tail)
    conclusion <- if (chi_squared_stat < lower_critical_value || chi_squared_stat > upper_critical_value) {
      "Reject H0 (Variance differs from hypothesized value)"
    } else {
      "Fail to reject H0 (Variance is equal to hypothesized value)"
    }
  } else if (tail_type == "right-tailed") {
    # Right-tailed test: Calculate only the upper critical value
    upper_critical_value <- qchisq(1 - alpha, df = n - 1)  
    conclusion <- if (chi_squared_stat > upper_critical_value) {
      "Reject H0 (Variance is greater than hypothesized value)"
    } else {
      "Fail to reject H0 (Variance is less than or equal to hypothesized value)"
    }
    # Only calculate lower critical value for two-tailed or left-tailed tests
    lower_critical_value <- NA  # No need for lower critical value here
  } else if (tail_type == "left-tailed") {
    # Left-tailed test: Calculate only the lower critical value
    lower_critical_value <- qchisq(alpha, df = n - 1)
    conclusion <- if (chi_squared_stat < lower_critical_value) {
      "Reject H0 (Variance is less than hypothesized value)"
    } else {
      "Fail to reject H0 (Variance is greater than or equal to hypothesized value)"
    }
    # No need for upper critical value for left-tailed test
    upper_critical_value <- NA
  } else {
    stop("Invalid tail_type. Choose from 'left-tailed', 'right-tailed', or 'two-tailed'.")
  }
  
  # Output the results
  cat("Sample size (n):", n, "\n")
  cat("Sample variance (s^2):", round(sample_variance, 4), "\n")
  cat("Hypothesized variance (σ^2):", test_sigma^2, "\n")
  cat("Chi-squared statistic:", round(chi_squared_stat, 4), "\n")
  if (!is.na(lower_critical_value)) {
    cat("Lower critical value:", round(lower_critical_value, 4), "\n")
  }
  if (!is.na(upper_critical_value)) {
    cat("Upper critical value:", round(upper_critical_value, 4), "\n")
  }
  cat("Conclusion:", conclusion, "\n")
  
  # Return the test results
  return(list(chi_squared_stat = chi_squared_stat, 
              lower_critical_value = lower_critical_value, 
              upper_critical_value = upper_critical_value, 
              conclusion = conclusion))
}

# Data set
x <- c(142.8, 135.0, 157.5, 148.4, 135.9, 153.4, 149.0, 130.2,
       156.0, 189.7, 151.6, 156.5, 123.8, 152.9, 118.4, 145.8)

# Hypothesized variance (σ^2) and significance level (alpha)
test_sigma <- 29  # Hypothesized standard deviation (σ)
alpha <- 0.05  # Significance level

# Perform the one-sample variance test (two-tailed)
result_two_tailed <- one_sample_variance_test(x, test_sigma, alpha, tail_type = "two-tailed")

# Perform the one-sample variance test (right-tailed)
result_right_tailed <- one_sample_variance_test(x, test_sigma, alpha, tail_type = "right-tailed")

# Perform the one-sample variance test (left-tailed)
result_left_tailed <- one_sample_variance_test(x, test_sigma, alpha, tail_type = "left-tailed")

# The asymmetry of the Chi-squared distribution causes the critical values for different types of tests 
# (right-tailed, left-tailed, and two-tailed) to differ. This skewness explains why the rejection regions are 
# positioned differently for each type of hypothesis test, and why conclusions can change depending on the 
# test direction (left, right, or two-tailed). For right-tailed tests, we focus on the right tail, while for 
# left-tailed tests, we focus on the left tail. For a two-tailed test, both tails are considered.


# Question 5
# Given dataset
x <- c(176.9, 158.3, 152.1, 158.8, 172.4, 169.8, 159.7, 162.7,
       156.6, 174.5, 184.4, 165.2, 147.8, 177.8, 160.1, 161.5)

# Perform the Wilcoxon signed-rank test
wilcoxon_result <- wilcox.test(x, mu = 160, alternative = "less", conf.int = TRUE, conf.level = 0.95)

# Print the results summary
cat("Wilcoxon Signed-Rank Test Results:\n")
cat("W statistic:", wilcoxon_result$statistic, "\n")
cat("P-value:", wilcoxon_result$p.value, "\n")
cat("Confidence Interval for the median: ", 
    round(wilcoxon_result$conf.int[1], 4), " to ", 
    round(wilcoxon_result$conf.int[2], 4), "\n")

# Interpretation of the results
if (wilcoxon_result$p.value < 0.05) {
  cat("Conclusion: Reject H0 (The median is less than 160). There is enough evidence to suggest that the median is less than 160.\n")
} else {
  cat("Conclusion: Fail to reject H0 (The median is not less than 160). There is insufficient evidence to suggest that the median is less than 160.\n")
}


# Exercise III - Two sample tests


two_sample_Z_test <- function(x1, x2, sigma_x1, sigma_x2, alpha, null_hypothesis = "equal") {
  # Calculate sample means and sizes
  x1_bar <- mean(x1)
  x2_bar <- mean(x2)
  n1 <- length(x1)
  n2 <- length(x2)
  
  # Calculate the Z-statistic based on the chosen null hypothesis
  if (null_hypothesis == "equal") {
    null_value <- 0
  } else if (null_hypothesis == "greater") {
    null_value <- 0  # Tests if mu1 >= mu2
  } else if (null_hypothesis == "less") {
    null_value <- 0  # Tests if mu1 <= mu2
  } else {
    stop("Invalid null hypothesis. Choose 'equal', 'greater', or 'less'.")
  }
  
  Z <- (x1_bar - x2_bar - null_value) / sqrt((sigma_x1^2 / n1) + (sigma_x2^2 / n2))
  
  # Find the critical Z-value based on the significance level and type of test
  if (null_hypothesis == "equal") {
    # Two-tailed test
    Z_critical <- qnorm(1 - alpha / 2)
  } else if (null_hypothesis == "greater") {
    # Right-tailed test
    Z_critical <- qnorm(1 - alpha)
  } else if (null_hypothesis == "less") {
    # Left-tailed test
    Z_critical <- qnorm(alpha)
  }
  
  # Calculate p-value for the test
  if (null_hypothesis == "equal") {
    p_value <- 2 * (1 - pnorm(abs(Z)))  # Two-tailed p-value
  } else if (null_hypothesis == "greater") {
    p_value <- 1 - pnorm(Z)  # Right-tailed p-value
  } else if (null_hypothesis == "less") {
    p_value <- pnorm(Z)  # Left-tailed p-value
  }
  
  # Conclusion based on the Z-statistic and critical value
  if (abs(Z) > Z_critical) {
    conclusion <- paste("Reject H0: There is enough evidence to reject the null hypothesis (", null_hypothesis, ").", sep = "")
  } else {
    conclusion <- paste("Fail to reject H0: There is not enough evidence to reject the null hypothesis (", null_hypothesis, ").", sep = "")
  }
  
  # Output results
  return(list(
    Z_statistic = Z,
    Z_critical = Z_critical,
    p_value = p_value,
    Conclusion = conclusion
  ))
}

x1 = c( 258.0, 271.5, 189.1, 216.5, 237.2, 222.0, 231.3, 181.7, 220.0, 179.3, 238.1, 217.7,
        246.2, 241.5, 233.8, 222.3, 199.2, 167.9, 216.2, 240.4, 235.3, 187.0, 233.7, 214.7,
        174.6, 246.3, 185.7, 207.0, 244.3, 237.7, 245.2, 228.3, 201.8, 218.3, 242.7, 213.8,
        231.9, 257.3, 208.4, 250.7, 198.3, 206.7, 259.7, 253.3, 200.3, 196.6, 210.6, 257.6,
        173.5, 267.5, 167.2, 227.1, 172.1, 197.6, 256.9, 203.7, 195.1, 237.4, 210.2, 208.8,
        218.0, 205.1, 241.1, 216.8, 223.6, 191.0, 225.9, 215.1, 233.1, 243.0)


x2 = c( 221.0, 213.0, 199.3, 211.2, 225.2, 229.1, 253.9, 194.6, 243.0, 221.9, 230.9, 221.1,
        206.7, 217.2, 215.8, 203.0, 234.0, 196.3, 235.8, 234.3, 244.7, 248.8, 200.5, 232.0,
        233.3, 220.6, 289.2, 244.9, 230.8, 182.9, 199.3, 263.2, 220.6, 266.7, 258.0, 243.9,
        178.1, 200.7, 270.2, 224.4, 222.4, 234.6, 296.7, 202.3, 277.9, 204.3, 221.1, 257.0,
        243.4, 239.4, 230.0, 263.5, 241.3, 216.6, 227.9, 230.1, 230.5, 188.6, 289.3, 234.4,
        267.5, 256.0, 246.5, 210.5, 270.6, 295.5, 195.8, 235.3, 245.4, 245.4)

sigma_x1 <- 24.6  # Standard deviation for sample 1
sigma_x2 <- 27.8  # Standard deviation for sample 2
alpha <- 0.05     # Significance level

# Perform the two-sample Z test for different null hypotheses

# Test for equal means (two-tailed)
result_equal <- two_sample_Z_test(x1, x2, sigma_x1, sigma_x2, alpha, null_hypothesis = "equal")

# Test for greater means (right-tailed)
result_greater <- two_sample_Z_test(x1, x2, sigma_x1, sigma_x2, alpha, null_hypothesis = "greater")

# Test for less means (left-tailed)
result_less <- two_sample_Z_test(x1, x2, sigma_x1, sigma_x2, alpha, null_hypothesis = "less")

# # Print the results for each test
# cat("Result for Two-Tailed Test (Equal Means): \n")
# cat("Z-statistic:", result_equal$Z_statistic, "\n")
# cat("Critical Z-value:", result_equal$Z_critical, "\n")
# cat("P-value:", result_equal$p_value, "\n")
# cat("Conclusion:", result_equal$Conclusion, "\n\n")
# 
# cat("Result for Right-Tailed Test (Greater Means): \n")
# cat("Z-statistic:", result_greater$Z_statistic, "\n")
# cat("Critical Z-value:", result_greater$Z_critical, "\n")
# cat("P-value:", result_greater$p_value, "\n")
# cat("Conclusion:", result_greater$Conclusion, "\n\n")

cat("Result for Left-Tailed Test (Less Means): \n")
cat("Z-statistic:", result_less$Z_statistic, "\n")
cat("Critical Z-value:", result_less$Z_critical, "\n")
cat("P-value:", result_less$p_value, "\n")
cat("Conclusion:", result_less$Conclusion, "\n")


# Question 2

perform_t_test <- function(test_type, x, y = NULL, paired = FALSE, conf.level = 0.95) {
  alpha <- 0.05
  cat("Significance level (α):", alpha, "\n")
  
  if (test_type == "welch") {
    result <- t.test(x, y, alternative = "two.sided", var.equal = FALSE)
    cat("Welch's Two Sample t-test\n")
  } else if (test_type == "paired") {
    result <- t.test(x, y, alternative = "two.sided", paired = TRUE, conf.level = conf.level)
    cat("Paired Sample t-test\n")
  } else {
    cat("Invalid test type. Use 'welch' or 'paired'.\n")
    return(NULL)
  }
  
  cat("p-value:", result$p.value, "\n")
  cat("95% Confidence Interval:", result$conf.int[1], "to", result$conf.int[2], "\n")
  
  if (result$p.value < alpha) {
    cat("Conclusion: Reject H₀ → Significant difference between groups.\n")
  } else {
    cat("Conclusion: Fail to reject H₀ → No significant difference between groups.\n")
  }
}

# Welch's t-test
Xvar <- c(4.95, 5.37, 4.70, 4.96, 4.72, 5.17, 5.28, 5.12, 5.26, 5.48)
Yvar <- c(4.65, 4.86, 4.57, 4.56, 4.96, 4.63, 5.04, 4.92, 5.37, 4.58, 4.26, 4.40)
perform_t_test("welch", Xvar, Yvar)

# Paired t-test
data_before <- c(95, 106, 79, 71, 90, 79, 71, 77, 103, 103, 92, 63, 82, 76)
data_after <- c(97, 116, 82, 81, 82, 86, 107, 86, 94, 91, 85, 98, 91, 87)
perform_t_test("paired", data_before, data_after)

# Interpretation:
# > Welch Two Sample t-test: There is a significant difference between the two groups (Xvar and Yvar), 
# as indicated by a p-value of 0.006749 and a confidence interval that does not contain 0.
# > Paired t-test: There is no significant difference between the means of the two sets of data 
# (data_before and data_after), as indicated by a p-value of 0.101 and a confidence interval that includes 0.
# > These results suggest that while the two independent groups (Xvar and Yvar) differ significantly, the 
# paired data (data_before and data_after) does not show a significant difference after treatment.

# Question 3
two_sample_proportion_test <- function(test_type, x, n = NULL) {
  alpha <- 0.05
  cat("Significance level (α):", alpha, "\n")
  
  if (test_type == "proportion") {
    result <- prop.test(x, n, alternative = "two.sided", correct = TRUE)
    cat("p-value:", result$p.value, "\n")
    
    if (result$p.value < alpha) {
      cat("Conclusion: Reject H₀ → Significant difference in proportions.\n")
    } else {
      cat("Conclusion: Fail to reject H₀ → No significant difference in proportions.\n")
    }
    
  } else if (test_type == "fisher") {
    result <- fisher.test(x, alternative = "two.sided", conf.int = TRUE)
    cat("p-value:", result$p.value, "\n")
    
    if (result$p.value < alpha) {
      cat("Conclusion: Reject H₀ → Significant association between the groups.\n")
    } else {
      cat("Conclusion: Fail to reject H₀ → No significant association between the groups.\n")
    }
    
  } else {
    cat("Invalid test type. Use 'proportion' or 'fisher'.\n")
  }
}

# (a) Two-Sample Proportion Test
two_sample_proportion_test("proportion", x = c(520, 550), n = c(600, 600))

# (b) Fisher’s Exact Test
data_matrix <- matrix(c(11, 17, 42, 39), nrow = 2, byrow = TRUE)
two_sample_proportion_test("fisher", x = data_matrix)


# Question 4
two_sample_variance_test <- function(x, y, alpha = 0.05) {
  # Sample sizes
  n1 <- length(x)
  n2 <- length(y)
  
  # Sample variances
  var1 <- var(x)
  var2 <- var(y)
  
  # F statistic: larger variance / smaller variance
  F_stat <- max(var1, var2) / min(var1, var2)
  
  # Degrees of freedom
  df1 <- n1 - 1
  df2 <- n2 - 1
  
  # p-value (two-tailed)
  p_value <- 2 * min(
    pf(F_stat, df1, df2, lower.tail = FALSE),
    pf(F_stat, df1, df2, lower.tail = TRUE)
  )
  
  # Print results
  cat("F-statistic:", F_stat, "\n")
  cat("p-value:", p_value, "\n")
  cat("Significance level (alpha):", alpha, "\n")
  
  # Decision
  if (p_value < alpha) {
    cat("Conclusion: Reject H₀ → Variances are significantly different.\n")
  } else {
    cat("Conclusion: Fail to reject H₀ → No significant difference in variances.\n")
  }
}
# Provided datasets
x <- c(1067.7, 984.3, 998.8, 1025.9, 1060.9, 959.1, 1013.8,
       1047.0, 987.8, 1051.0, 885.2, 1049.5, 1098.2, 1001.5,
       1011.1, 991.6)

y <- c(957.6, 981.8, 1096.5, 984.4, 1074.3, 929.4, 1056.0,
       1012.3, 1040.7, 1099.5, 1006.1, 1064.3, 865.6, 944.4,
       1091.8, 952.1)

# Perform the variance test with alpha = 0.05
two_sample_variance_test(x, y, alpha = 0.05)


# Question 5

wilcoxon_signed_rank_test <- function(pre, post, alpha = 0.05) {
  result <- wilcox.test(pre, post, 
                        alternative = "greater", 
                        paired = TRUE, 
                        conf.level = 1 - alpha,
                        exact = FALSE)  # Use normal approximation for p-value when ties are present
  
  cat("Wilcoxon Signed-Rank Test Result:\n")
  cat("Test Statistic (V):", result$statistic, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Significance level (alpha):", alpha, "\n")
  
  if (result$p.value < alpha) {
    cat("Conclusion: Reject H₀ → Post_therapy values are significantly lower than Pre_therapy.\n")
  } else {
    cat("Conclusion: Fail to reject H₀ → No significant difference detected.\n")
  }
}

# Example usage:
Pre_therapy <- c(74, 72, 62, 58, 59, 65, 54, 63, 80, 66, 65, 64, 79, 60)
Post_therapy <- c(79, 55, 53, 53, 74, 55, 64, 55, 39, 44, 37, 68, 54, 54)

wilcoxon_signed_rank_test(Pre_therapy, Post_therapy)


# Question 6

wilcoxon_rank_sum_test <- function(group1, group2, alpha = 0.05) {
  result <- wilcox.test(group2, group1,  # placebo vs drug
                        alternative = "less", 
                        paired = FALSE,
                        conf.level = 1 - alpha,
                        exact = FALSE)  # handles ties if present
  
  cat("Wilcoxon Rank Sum Test (Mann-Whitney U Test):\n")
  cat("Test Statistic (W):", result$statistic, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Significance level (alpha):", alpha, "\n")
  
  if (result$p.value < alpha) {
    cat("Conclusion: Reject H₀ → Placebo group has significantly smaller mean than Drug group.\n")
  } else {
    cat("Conclusion: Fail to reject H₀ → No significant difference in means between groups.\n")
  }
}

# Sample Data
drug <- c(31.7, 75.0, 101.1, 60.5, 62.8, 59.3, 58.9, 91.3, 99.1, 52.0, 39.1)
placebo <- c(59.3, 72.7, 100.5, 64.7, 69.0, 72.7, 69.6, 97.4, 100.6, 65.1, 65.7)

# Function Call
wilcoxon_rank_sum_test(drug, placebo)


# Question 7

kruskal_wallis_test <- function(x, y, alpha = 0.05) {
  result <- kruskal.test(x ~ y)  # Kruskal-Wallis test
  
  cat("Kruskal-Wallis Test Result:\n")
  cat("Test Statistic (H):", result$statistic, "\n")
  cat("p-value:", result$p.value, "\n")
  cat("Significance level (alpha):", alpha, "\n")
  
  if (result$p.value < alpha) {
    cat("Conclusion: Reject H₀ → Significant difference between the groups.\n")
  } else {
    cat("Conclusion: Fail to reject H₀ → No significant difference between the groups.\n")
  }
}

# Data: Group-wise values
group_1 <- c(220, 214, 203, 184, 186, 200, 165)
group_2 <- c(262, 193, 225, 200, 164, 266, 179)
group_3 <- c(272, 192, 190, 208, 231, 235, 141)
group_4 <- c(190, 255, 247, 278, 230, 269, 289)

# Create the corresponding labels for each group
x_values <- c(group_1, group_2, group_3, group_4)
y_labels <- rep(1:4, times = c(length(group_1), length(group_2), length(group_3), length(group_4)))

# Perform the Kruskal-Wallis test
kruskal_wallis_test(x_values, y_labels)


# Question 8

chi_square_gof_test <- function(observed, expected, alpha = 0.05) {
  # Calculate the Chi-square statistic
  chi_square_statistic <- sum((observed - expected)^2 / expected)
  
  # Degrees of freedom = number of categories - 1
  df <- length(observed) - 1
  
  # Find the critical value using the qchisq function
  critical_value <- qchisq(1 - alpha, df)
  
  # Calculate the p-value using the pchisq function
  p_value <- 1 - pchisq(chi_square_statistic, df)
  
  # Print the results
  cat("Chi-Square Goodness of Fit Test Result:\n")
  cat("Chi-Square Statistic:", chi_square_statistic, "\n")
  cat("Critical Value:", critical_value, "\n")
  cat("p-value:", p_value, "\n")
  
  # Conclusion based on p-value
  if (p_value < alpha) {
    cat("Conclusion: Reject H₀ → Significant difference between observed and expected frequencies.\n")
  } else {
    cat("Conclusion: Fail to reject H₀ → No significant difference between observed and expected frequencies.\n")
  }
}

# Data: Observed and Expected frequencies
observed <- c(32, 82, 77, 49)
expected <- c(40, 80, 80, 40)

# Perform the Chi-Square GoF Test
chi_square_gof_test(observed, expected)


# Exercise 4 

# Question 1

# a)

# Load necessary libraries
library(ggplot2)
# Load dplyr package
library(dplyr)
# Read the titanic data
titanicData <- read.csv("/home/ibab/R/Lab14/titanic.csv")
# Ensure passenger_class is treated as a factor
titanicData$passenger_class <- as.factor(titanicData$passenger_class)
print(levels(titanicData$passenger_class))  # To check levels are 1st, 2nd, 3rd
# Set up the plotting area to have 3 rows and 1 column
par(mfrow = c(3,1))  # 3 rows, 1 column
# Plot histogram for 1st class
hist(titanicData$age[titanicData$passenger_class == "1st"],
     breaks=30,
     main="Age Distribution - 1st Class",
     xlab="Age",
     col="lightblue",
     border="black")
# Plot histogram for 2nd class
hist(titanicData$age[titanicData$passenger_class == "2nd"],
     breaks=30,
     main="Age Distribution - 2nd Class",
     xlab="Age",
     col="lightgreen",
     border="black")

# Plot histogram for 3rd class
hist(titanicData$age[titanicData$passenger_class == "3rd"],
     breaks=30,
     main="Age Distribution - 3rd Class",
     xlab="Age",
     col="lightpink",
     border="black")

# Reset plotting area back to default (optional)
par(mfrow = c(1,1))

# Hard to infer anything about normaility or equality of variance


# b)

# Group the data by passenger class and calculate mean and standard deviation
titanic_by_passenger_class <- group_by(titanicData, passenger_class)

# Summarize mean and standard deviation of age by group
summary_stats <- summarise(titanic_by_passenger_class, 
                           group_mean = mean(age, na.rm=TRUE),
                           group_sd = sd(age, na.rm=TRUE))

# Print the summary statistics
print(summary_stats)
cat("The standard deviations across passenger classes (14.9, 13.0, 11.3) are reasonably similar.\n")
cat("Thus, the assumption of equal variance needed for ANOVA is satisfied.\n")


# c)
# Fit the ANOVA model
lmresults <- lm(age ~ passenger_class, data=titanicData)
# Perform ANOVA
anova_results <- anova(lmresults)
# Print the ANOVA results
print(anova_results)

# Conclusion Based on the Example:
# The F-value (75.903) is large, indicating that there is a significant difference between the passenger classes in terms of age.
# The p-value (< 2.2e-16) is extremely small, confirming that the difference is statistically significant.
# Thus, we reject the null hypothesis and conclude that there are significant differences in age among the passenger classes.

# '*' (p ≤ 0.001)**: Very highly significant, strong evidence against the null hypothesis.
# '' (0.001 < p ≤ 0.01)**: Significant, strong evidence against the null hypothesis.
# '*' (0.01 < p ≤ 0.05): Marginally significant, some evidence against the null hypothesis.
# '.' (0.05 < p ≤ 0.1): Borderline significant, weak evidence against the null hypothesis.
# ' ' (p > 0.1): Not significant, not enough evidence to reject the null hypothesis.

# d)
# Perform Tukey's HSD test
tukey_results <- TukeyHSD(aov(lmresults))

# Print the results
print(tukey_results)

# Interpretation of Tukey's Results:
# The p-values for all comparisons (2nd vs 1st, 3rd vs 1st, 3rd vs 2nd) are very small (0.0000000 and 0.0116695), which means that all pairwise comparisons are statistically significant. This confirms that the mean ages of the passengers in different classes are significantly different from each other.
# 2nd vs 1st: Mean age of 2nd class is significantly lower than 1st class.
# 3rd vs 1st: Mean age of 3rd class is significantly lower than 1st class.
# 3rd vs 2nd: Mean age of 3rd class is significantly lower than 2nd class


# e)
# Perform Kruskal-Wallis test
kruskal_results <- kruskal.test(age ~ passenger_class, data=titanicData)

# Print the result
print(kruskal_results)

# The Kruskal-Wallis test is a non-parametric alternative to ANOVA that doesn't assume normality. It tests whether there are significant differences between the medians of the groups.
# Interpretation of Kruskal-Wallis Results:
# Kruskal-Wallis chi-squared = 116.08, df = 2, p-value < 2.2e-16:
# The p-value is extremely small (< 2.2e-16), indicating that there is a significant difference between the passenger classes in terms of age.
# Since the p-value is much smaller than the significance level (usually 0.05), we reject the null hypothesis that the median ages of the groups are equal. This suggests that the ages of passengers from different classes are significantly different.


# Conclusion: The Chi-Square Goodness of Fit test suggests that there is no significant 
# difference between the observed and expected frequencies, while both the Tukey HSD 
# and Kruskal-Wallis tests suggest significant differences in the mean or median ages of passengers in different classes.


# Question 2

# a)

# Load the data
cuckooData <- read.csv("/home/ibab/R/Lab14/cuckooeggs.csv")
str(cuckooData)
dim(cuckooData)
head(cuckooData)
len_of_unique<-length(unique(cuckooData$host_species))
# Plot histograms for each host species
par(mfrow = c(len_of_unique/2, 2))  # Set layout

for (species in unique(cuckooData$host_species)) {
  hist(cuckooData$egg_length[cuckooData$host_species == species],
       main = species, xlab = "Egg Length", col = "lightblue",breaks=10)
}

par(mfrow = c(1,1)) 

# (b)

# Group the data by host species
group_cuckoo_byspecies <- group_by(cuckooData, host_species)

# Summarize mean and standard deviation by species
summary_stats <- summarise(group_cuckoo_byspecies, 
                           group_mean = mean(egg_length, na.rm = TRUE),
                           group_sd = sd(egg_length, na.rm = TRUE))

# Print the summary statistics
print(summary_stats)

# Comment about variances
cat("The standard deviations across host species classes are similar.\n")
cat("Thus, the assumption of equal variance needed for ANOVA is satisfied.\n")

# (c)
# Fit ANOVA model
anova_model <- aov(egg_length ~ host_species, data = cuckooData)
summary(anova_model)
# Conclusion:
# If p-value < 0.05 → Means are different across species.

# (d)
# Perform Tukey's HSD (Honestly Significant Difference) test
tukey_result <- TukeyHSD(anova_model)
# Print results
print(tukey_result)


# Significant differences (p < 0.05):
# Meadow Pipit vs Hedge Sparrow (p = 0.043) → Significant
# Wren vs Hedge Sparrow (p ≈ 0.0000006) → Highly Significant
# Tree Pipit vs Meadow Pipit (p ≈ 0.047) → Significant
# Wren vs Meadow Pipit (p ≈ 0.000486) → Highly Significant
# Wren vs Pied Wagtail (p ≈ 0.0000070) → Highly Significant
# Wren vs Robin (p ≈ 0.000318) → Highly Significant
# Wren vs Tree Pipit (p ≈ 0.0000006) → Highly Significant

# Not significant (p > 0.05):
# Pied Wagtail vs Hedge Sparrow
# Robin vs Hedge Sparrow
# Tree Pipit vs Hedge Sparrow
# Pied Wagtail vs Meadow Pipit
# Robin vs Meadow Pipit
# Robin vs Pied Wagtail
# Tree Pipit vs Pied Wagtail
# Tree Pipit vs Robin

# Question 3
# Read the data
malariaData <- read.csv("/home/ibab/R/Lab14/malaria vs maize.csv")

# Quick look
print(head(malariaData))

# (a) Multiple histograms for original incidence rates
# Subset data
low <- malariaData$incidence_rate_per_ten_thousand[malariaData$maize_yield == "Low"]
medium <- malariaData$incidence_rate_per_ten_thousand[malariaData$maize_yield == "Medium"]
high <- malariaData$incidence_rate_per_ten_thousand[malariaData$maize_yield == "High"]

# Plot histograms
par(mfrow = c(1, 3))  # 1 row, 3 columns
hist(low, 
     main = "Low Maize Yield", 
     xlab = "Malaria Incidence Rate", 
     col = "lightblue", 
     breaks = 10)

hist(medium, 
     main = "Medium Maize Yield", 
     xlab = "Malaria Incidence Rate", 
     col = "lightgreen", 
     breaks = 10)

hist(high, 
     main = "High Maize Yield", 
     xlab = "Malaria Incidence Rate", 
     col = "lightpink", 
     breaks = 10)

# Reset layout
par(mfrow = c(1,1))

# (b) Calculate standard deviations
sd_low <- sd(low)
sd_medium <- sd(medium)
sd_high <- sd(high)

cat("Standard deviations (original data):\n")
cat("Low:", sd_low, "\n")
cat("Medium:", sd_medium, "\n")
cat("High:", sd_high, "\n")

# (b) Check ANOVA assumptions
# If SDs are very different, variance homogeneity assumption is violated!

# (c) Log-transform incidence rates
malariaData$log_incidence <- log(malariaData$incidence_rate_per_ten_thousand)

# Subset log data
low_log <- malariaData$log_incidence[malariaData$maize_yield == "Low"]
medium_log <- malariaData$log_incidence[malariaData$maize_yield == "Medium"]
high_log <- malariaData$log_incidence[malariaData$maize_yield == "High"]

# Plot histograms for log-transformed data
par(mfrow = c(1, 3))
hist(low_log, 
     main = "Low Maize Yield (Log)", 
     xlab = "Log Malaria Incidence Rate", 
     col = "lightblue", 
     breaks = 10)

hist(medium_log, 
     main = "Medium Maize Yield (Log)", 
     xlab = "Log Malaria Incidence Rate", 
     col = "lightgreen", 
     breaks = 10)

hist(high_log, 
     main = "High Maize Yield (Log)", 
     xlab = "Log Malaria Incidence Rate", 
     col = "lightpink", 
     breaks = 10)

# Reset layout
par(mfrow = c(1,1))

# (c) Calculate standard deviations of log data
sd_low_log <- sd(low_log)
sd_medium_log <- sd(medium_log)
sd_high_log <- sd(high_log)

cat("Standard deviations (log-transformed data):\n")
cat("Low (log):", sd_low_log, "\n")
cat("Medium (log):", sd_medium_log, "\n")
cat("High (log):", sd_high_log, "\n")

# (d) Test association
# ANOVA on log-transformed data
anova_result <- aov(log_incidence ~ maize_yield, data = malariaData)
summary(anova_result)


# Question 4

# (a)

# Load the data
circadianData <- read.csv("/home/ibab/R/Lab14/circadian mutant health.csv")
# Subset data for each genotype
tim01_data <- circadianData$days_to_death[circadianData$genotype == "tim01"]
rescued_data <- circadianData$days_to_death[circadianData$genotype == "tim01 (rescued)"]
wildtype_data <- circadianData$days_to_death[circadianData$genotype == "wild type"]

# Set up 1 row 3 columns layout for plots
par(mfrow = c(1, 3))

# Plot histograms
hist(tim01_data, 
     main = "tim01 Mutant", 
     xlab = "Days to Death", 
     col = "red", 
     breaks = 10, 
     xlim = c(0, 22))

hist(rescued_data, 
     main = "tim01 Rescued", 
     xlab = "Days to Death", 
     col = "blue", 
     breaks = 10, 
     xlim = c(0, 22))

hist(wildtype_data, 
     main = "Wild Type", 
     xlab = "Days to Death", 
     col = "green", 
     breaks = 10, 
     xlim = c(0, 22))

# Reset layout
par(mfrow = c(1,1))


# (b)
kruskal.test(days_to_death ~ genotype, data = circadianData)
# If the p-value is < 0.05, it means there is a statistically significant difference in median survival times between at least two groups.
