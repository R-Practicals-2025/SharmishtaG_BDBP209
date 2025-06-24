# BDBP209: Biostatistics and R Programming Laboratory

This repository contains lab exercises, notes, and code for the BDBP209 R Programming Laboratory conducted at IBAB, Bengaluru. The course introduces essential R programming skills and foundational statistical analysis methods through hands-on exercises using biological and real-world datasets.

## Course Overview

The laboratory is designed to help students:
- Learn the basics of R programming
- Explore data types, data structures, and file handling
- Perform data visualization using R base and ggplot2
- Conduct statistical analyses using built-in functions
- Work with real datasets to perform regression and clustering
- Understand distributions, hypothesis testing, and ANOVA

## Labs Covered

### Fundamentals and Data Handling
- Lab 1: R environment, operators, objects, and functions
- Lab 2 & 3: Vectors, sequences, missing values, matrices
- Lab 5 & 6: Data frames, factors, apply family, file operations
- Lab 7: Matrix algebra, set operations, loops, and conditional statements
- Lab 8: String manipulation and DNA palindrome detection

### Visualization and Distributions
- Lab 9: Base R plots, histograms, density plots, and barplots
- Lab 10 & 11: Probability distributions (binomial, Poisson, Gaussian, etc.), sampling, Central Limit Theorem, ROC curves

### Statistical Analysis
- Lab 14:
  - Error bars, covariance, and correlation
  - One-sample and two-sample tests (Z-test, t-test, Wilcoxon)
  - Proportion tests, variance tests, Kruskal-Wallis test
  - Chi-square Goodness-of-Fit test
  - One-way ANOVA

### Regression and Clustering
- Lab 15:
  - Linear regression (simple, multiple, polynomial)
  - Non-linear regression using `nls()`
  - Hierarchical clustering with gene expression data
  - Dendrograms and heatmaps

## Datasets Used

- winequality-red.csv and winequality-white.csv
- trees (built-in dataset)
- spellman-wide.csv (yeast gene expression)
- synthetic and user-defined datasets

## Requirements

Install required R packages before starting:

```r
install.packages(c("ggplot2", "dplyr", "readr", "tidyr", "stringr", "gtools", "lattice", "dendextend", "RColorBrewer"))
