# LAB 7 

# Ex 1
# (i)
amat <- matrix(seq(10, 120, by = 10), nrow = 3, ncol = 4, byrow = TRUE)
print(amat)
amat2 <- matrix(seq(10, 120, by = 10), nrow = 3, ncol = 4, byrow = FALSE)
print(amat2)
identical(amat2, t(amat))
#Both matrices are transposes of each other
#Difference between amat and amat2 is that with byrow=FALSE the elements are distributed column-wise, first column is filled, then second column is filled and so on.

# (ii)
rownames(amat)<-c("R1","R2","R3")
dimnames(amat)<-list(rownames(amat),c("C1", "C2", "C3", "C4"))
print(amat)

# (iii)
A <- matrix(c(2,5,7,3, 1,8,9,10, 1,12,5,10, 4,17,15,11), nrow=4, ncol=4, byrow=TRUE)
print("Matrix A:")
print(A)
B <- matrix(c(12,5,3,17, 1,18,9,10, 1,12,5,10, 4,15,15,4), nrow=4, ncol=4, byrow=TRUE)
print("Matrix B:")
print(B)
element_wise_product <- A * B
print("Element-wise Multiplication Result:")
print(element_wise_product)
matrix_product <- A %*% B
print("Matrix-Matrix Multiplication Result:")
print(matrix_product)

# (iv)
X <- c(5, 6, 8, 9)
print(class(X))
Y <- c(8, 10, 12, 5)
print(class(Y))
outer_product <- X %o% Y
print("Outer Product:")
print(outer_product)
inner_product <- sum(X * Y)
print("Inner Product:")
print(inner_product)
# Alternative - for 1 D vector its the same as inner product
inner_product <- X %*% Y
print(inner_product)
print(class(outer_product))
print(class(inner_product))

# (v)
X <- c(5, 6, 8, 9)
diag_matrix <- diag(X)
print("Diagonal Matrix:")
print(diag_matrix)

# (vi)
A <- matrix(c(2,5,7,3, 1,8,9,10, 1,12,5,10, 4,17,15,11), nrow=4, ncol=4, byrow=TRUE)
diagonal_A <- diag(A)
print("Diagonal Elements of A:")
print(diagonal_A)

# (vii)
identity_matrix <- diag(6)
print("6x6 Identity Matrix:")
print(identity_matrix)

# (viii)
A <- matrix(c(3,4,-2, 4,-5,1, 10,-6,5), nrow=3, ncol=3)
print("Matrix A:")
print(A)

# (ix)
B <- matrix(c(5,-3,13), nrow=3, ncol=1)
print("Matrix B:")
print(B)

# (x)
X <- solve(A, B) # solve(A, B) computes X=A ^−1*B, solving for the unknown vector 𝑋
print("Solution X:")
print(X)
print("Type of X:")
print(class(X))
print(typeof(X)) #how is it stored in memory
# x<-1l Ensures its an integer
# difference between attribute,mode,typeof,class

# Summary of Differences
# attributes() → Metadata (e.g., dim, names).
# mode() → High-level data type (e.g., numeric, character).
# typeof() → Low-level internal data structure (e.g., integer, double).
# class() → Object’s behavior and structure (e.g., matrix, data.frame).
# 
# Quick Analogy for Clarity
# Imagine an Excel sheet:
# attributes() → The sheet’s labels, column names, or dimensions.
# mode() → The data type inside cells (e.g., numbers or text).
# typeof() → The technical format (e.g., integer vs double).
# class() → Defines the sheet’s structure (e.g., table format).


# (xi)
Ainv <- solve(A)
print("Inverse of A:")
print(Ainv)

# Check if Ainv * A is an identity matrix
identity_check <- Ainv %*% A
print("Ainv * A (Should be Identity Matrix):")
print(identity_check)
all.equal(identity_check, diag(3))


# (xii)
results <- eigen(A)
print("Eigenvalues of A:")
print(results$values)

print("Eigenvectors of A:")
print(results$vectors)

print("Type of results:")
print(class(results))

eigenvector2 <- results$vectors[,2]  # Select the second eigenvector
multiplication_result <- A %*% eigenvector2
print("A * Second Eigenvector:")
print(multiplication_result)

# Compare with λ * v (theoretical property of eigenvectors)
lambda2 <- results$values[2]  # Second eigenvalue
expected_result <- lambda2 * eigenvector2
print("λ * Second Eigenvector:")
print(expected_result)

# For an eigenvector𝑣with eigenvalue 𝜆 we should get:
# 𝐴𝑣=𝜆𝑣 - satisfies eigen function
# Av=λv is the fundamental property of eigenvalues and eigenvectors.
# The matrix A stretches (or scales) the eigenvector by the corresponding eigenvalue.
# Eigenvectors are directional invariants under the transformation defined by A.

# Ex2
brain_cancer_data=read.csv("/Users/sharmishtaganesh/Desktop/Biostat_Lab/Lab_files/BrainCancer.csv", header=TRUE)
print(colnames(brain_cancer_data))
# Print first few rows to verify
print(head(brain_cancer_data))

# (i) Create a new column: (GTV^2 + time)
brain_cancer_data$new_column <- (brain_cancer_data$gtv^2) + brain_cancer_data$time

# (ii) Print row and column names
print("Row Names:")
print(rownames(brain_cancer_data))
print("Column Names:")
print(colnames(brain_cancer_data))

# (iii) Change row names to 'Row-1', 'Row-2', ...
rownames(brain_cancer_data) <- paste("Row", 1:nrow(brain_cancer_data), sep="-")

# (iv) Remove the 'ki' column
brain_cancer_data$ki <- NULL

# Print modified data
print(head(brain_cancer_data))


# Ex3 
# (i) Install the readxl package (only needed once)
install.packages("readxl")

# (ii) Load the readxl package
library(readxl)

# (iii) Read the Excel file (Replace <path_to_file.xlsx> with the actual file path)
data <- read_excel("<path_to_pone.0148733.s001.xlsx>", sheet = 1)

# (iv) Print column names and dimensions
print("Column Names:")
print(colnames(data))

print("Data Dimensions (Rows, Columns):")
print(dim(data))


# Ex4
# (i) Define vectors A and B
A <- c("a", "b", "c", "d", "e")
B <- c("d", "e", "f", "g")

# Print the vectors
print("Vector A:")
print(A)
print("Vector B:")
print(B)

# (ii) Union of A and B
union_AB <- union(A, B)
print("Union of A and B:")
print(union_AB)

# (iii) Intersection of A and B
intersection_AB <- intersect(A, B)
print("Intersection of A and B:")
print(intersection_AB)

# (iv) Difference of A - B
difference_A_B <- setdiff(A, B)
print("Difference A - B:")
print(difference_A_B)

# Difference of B - A
difference_B_A <- setdiff(B, A)
print("Difference B - A:")
print(difference_B_A)

# (v) Check if (A - B) ∪ (A ∩ B) ∪ (B - A) equals union(A, B)
combined_sets <- union(difference_A_B, union(intersection_AB, difference_B_A))
print("Is (A - B) ∪ (A ∩ B) ∪ (B - A) equal to union(A, B)?")
print(setequal(combined_sets, union_AB))

# different ways to perform intersection
# (vi) List elements of B present in A using two different methods
print("Elements of B present in A (Method 1 - Using %in% operator):")
print(B[B %in% A])

print("Elements of B present in A (Method 2 - Using intersect() function):")
print(intersect(B, A))

# (vii) List elements of A present in B
print("Elements of A present in B:")
print(A[A %in% B])

print(setdiff(setA,setdiff(setA,setB)))


# Ex5

# (i) Create vector and filter elements
vec <- c(8,10,12,7,14,16,2,4,9,19,20,3,6)

# (a) Values greater than 12
print("Values greater than 12:")
print(vec[vec > 12])

# (b) Values greater than 10 and less than 20
print("Values greater than 10 and less than 20:")
print(vec[vec > 10 & vec < 20])


# (ii) Create array A with NA values
A <- c(2,7,29,32,41,11,15,NA,NA,55,32,NA,42,109)

# Remove NA and keep values less than 100
A_filtered <- A[!is.na(A) & A < 100]
print("Array after removing NA and values less than 100:")
print(A_filtered)

# (iii) Replace NA with 0
A[is.na(A)] <- 0
print("Array after replacing NA with 0:")
print(A)


# (iv) Create vectors for gene names and gender
genes <- paste("gene", 1:7, sep="-")
gender <- c("M", "M", "F", "M", "F", "F", "M")


# (v) Enter experiment data as vectors
result1 <- c(12.3, 11.5, 13.6, 15.4, 9.4, 8.1, 10.0)
result2 <- c(22.1, 25.7, 32.5, 42.5, 12.6, 15.5, 17.6)
result3 <- c(15.5, 13.4, 11.5, 21.7, 14.5, 16.5, 12.1)
result4 <- c(14.4, 16.6, 45.0, 11.0, 9.7, 10.0, 12.5)
result5 <- c(12.2, 15.5, 17.4, 19.4, 10.2, 9.8, 9.0)
result6 <- c(13.3, 14.5, 21.6, 17.9, 15.6, 14.4, 12.0)
result7 <- c(11.0, 10.0, 12.2, 14.3, 23.3, 19.8, 13.4)


# (vi) Create dataframe
datframe <- data.frame(genes, gender, result1, result2, result3, result4, result5, result6, result7)

# (vii) Rename columns
colnames(datframe) <- c("GeneName", "Gender", "expt1", "expt2", "expt3", "expt4", "expt5", "expt6", "expt7")

# Print the dataframe
print("DataFrame:")
print(datframe)


# (viii) Subset with expt2 values greater than 20
subset_expt2_gt20 <- subset(datframe, expt2 > 20)
print("Subset with expt2 > 20:")
print(subset_expt2_gt20)


# (ix) Subset with only Female gender
subset_female <- subset(datframe, Gender == "F")
print("Subset with only Female gender:")
print(subset_female)


# (x) Subset with Male gender and expt2 < 30
subset_male_expt2_lt30 <- subset(datframe, Gender == "M" & expt2 < 30)
print("Subset with Male gender and expt2 < 30:")
print(subset_male_expt2_lt30)


# Ex6
# (i) Determine the quadrant of an angle
find_quadrant <- function(angle) {
  angle <- angle %% 360  # Convert to within 0-360 range
  if (angle == 0 | angle == 360) {
    print("Positive X-axis")
  } else if (angle == 90) {
    print("Positive Y-axis")
  } else if (angle == 180) {
    print("Negative X-axis")
  } else if (angle == 270) {
    print("Negative Y-axis")
  } else if (angle > 0 & angle < 90) {
    print("First quadrant")
  } else if (angle > 90 & angle < 180) {
    print("Second quadrant")
  } else if (angle > 180 & angle < 270) {
    print("Third quadrant")
  } else {
    print("Fourth quadrant")
  }
}
find_quadrant(45)  

# (ii) Arrange three numbers in decreasing order without using built-in sorting
sort_descending <- function(a, b, c) {
  if (a >= b & a >= c) {
    if (b >= c) {
      print(c(a, b, c))
    } else {
      print(c(a, c, b))
    }
  } else if (b >= a & b >= c) {
    if (a >= c) {
      print(c(b, a, c))
    } else {
      print(c(b, c, a))
    }
  } else {
    if (a >= b) {
      print(c(c, a, b))
    } else {
      print(c(c, b, a))
    }
  }
}

sort_descending(12, 45, 7)  # Expected Output: c(45, 12, 7)


# (iii) Calculate journey ticket cost based on distance and age
calculate_ticket_cost <- function(distance, age) {
  if (distance <= 100) {
    cost <- 100
  } else if (distance <= 1000) {
    cost <- 100 + (distance - 100) * 1.5
  } else {
    cost <- 100 + (900 * 1.5) + (distance - 1000) * 2
  }
  
  # Applying discount based on age
  if (age > 60) {
    cost <- cost * 0.75  # 25% discount for senior citizens
  } else if (age < 6) {
    cost <- cost * 0.5   # 50% discount for children under 6
  }
  
  print(paste("Ticket cost:", cost, "Rs"))
}

calculate_ticket_cost(150, 65)  # Senior citizen, 150 km
calculate_ticket_cost(1200, 4)  # Child under 6, 1200 km
calculate_ticket_cost(500, 30)  # Regular adult, 500 km


# Ex7

# (i)
v<-c(9,-1,4,5,6,-88,-7,55)
k=1
for (i in v){
  if (i<0){
    v[k]=0
  }
  k=k+1
}
print(v)

# (ii)
stirling_factorial <- function(n) {
  if (n < 0) {
    stop("Factorial is not defined for negative numbers.")
  }
  if (n == 0) {
    return(1)  # 0! is 1
  }
  
  pi_term <- sqrt(2 * pi * n)
  power_term <- (n / exp(1))^n
  correction <- 1 + 1/(12*n) + 1/(288*n^2) - 139/(51840*n^3) - 571/(2488320*n^4)
  
  approx_factorial <- pi_term * power_term * correction
  return(approx_factorial)
}

# Test cases
print(stirling_factorial(5))   # Expected output: Close to 120 (5!)
print(stirling_factorial(10))  # Expected output: Close to 3.6288e+06 (10!)
print(stirling_factorial(0))   # Expected output: 1
print(stirling_factorial(20))  # Expected output: Close to 2.4329e+18 (20!)


# (iii)
sumOfdigits <- function(x) {
  sum <- 0  # Initialize sum
  
  while (x > 0) {  # Loop until x becomes 0
    rem <- x %% 10  # Get last digit
    x <- x %/% 10   # Remove last digit
    sum <- sum + rem  # Add to sum
  }
  
  return(sum)
}

# Test cases
print(sumOfdigits(10))   # Expected output: 1 (1+0)
print(sumOfdigits(1234)) # Expected output: 10 (1+2+3+4)
print(sumOfdigits(999))  # Expected output: 27 (9+9+9)
print(sumOfdigits(0))    # Expected output: 0


