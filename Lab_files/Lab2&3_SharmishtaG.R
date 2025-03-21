# Lab2: 10 Jan 2025

############### Q1 ###############

# Using round() to round numbers to 3 decimal places
round(123456789, digits = 3)      # Large number: R doesn't display decimals for such large integers
round(12.1343, digits = 3)        # Q1.1: Rounds to 12.134
round(123.12344, digits = 3)      # Q1.2: Rounds to 123.123
round(1234.12344, digits = 3)     # Q1.3: Rounds to 1234.123
x <- round(12345.12344, digits = 3)  # Q1.4: Rounds to 12345.123
print(x)                          # Default printing: limited to digits displayed
options(digits = 15)              # Q1.5: Changing the number of digits displayed in R's global settings
print(x)                          # Shows 12345.123 without numerical rounding noise

# Numerical noise example
y <- round(12345.12344, digits = 3)  
print(y)                          # Same result as x
options(digits = 20)              # Displaying 20 digits for clarity
print(y)                          # More digits are now visible to clarify the rounding process

# Using formatC to control precision and formatting explicitly
formatC(round(12345.12344, digits = 3), format = "f", digits = 3)  # Q1.6: Formats with exactly 3 decimal places
formatC(round(12345.12344, digits = 3), digits = 3)                # Q1.6: Removes excess formatting

# Exploring print behavior
print(1234.12344)                # Q1.7: Default behavior: prints with R's global options digits
print(1234.723, digits = 3)      # Q1.8: Rounds to nearest integer since digits <= length before the decimal
print(1234.723, digits = 5)      # Q1.8: Prints 1234.72300 as more precision allows decimal retention

# Demonstrating the effect of rounding large numbers
round(123456788.123, digits = 3)  # Q1.9: Rounds to 123456788.123
print(round(123456788.123, digits = 2), digits = 20)  # Q1.10: Prints 123456788.12 with extra significant digits
print(round(123456789.1234, digits = 4), digits = 20)  # Q1.10: Displays with greater detail using digits option

# Exploring the paste() function
paste("Hello World", sep = ",")          # Q1.11: Prints "Hello World" (sep is unnecessary here for one argument)
paste("Hello", "World", sep = ",")       # Q1.11: Concatenates "Hello,World"
paste(1:10)                              # Q1.12: Pastes numbers 1 to 10 as individual elements
paste(1:10)[4]                           # Q1.12: Extracts the 4th element of the pasted result, i.e., "4"
as.numeric(paste(1:10))                  # Q1.14: Converts the pasted result into numeric, i.e., numbers 1 to 10
paste(1:10, collapse = ".")              # Q1.15: Collapses all pasted elements into one string: "1.2.3...10"

# Combining strings and objects of different lengths
print(paste("Hello", 1:10, sep = "-"))   # Q1.16: Concatenates "Hello" with each number from 1 to 10 with "-" separator


############### Q2 ###############

# Generating simple sequences
print(0:10)                           # (i) Sequence from 0 to 10
print(15:5)                           # (ii) Sequence from 15 to 5 (descending)
print(seq(0, 1.5, 0.1))               # (iii) Sequence from 0 to 1.5 with step size 0.1
print(seq(6, 4, -0.2))                # (iv) Sequence from 6 to 4 with step size -0.2

# Define N vector
N <- c(55, 76, 92, 103, 84, 88, 121, 91, 65, 77, 99)  

# Generating sequences with 'seq'
print(seq(from = 0.04, by = 0.01, length = 11))   # (v) Sequence from 0.04 in steps of 0.01, length 11
print(seq(0.04, by = 0.01, along = N))           # (vi) Matching length of N
print(seq(from = 0.04, to = 0.14, along = N))    # (vii) Sequence with same start & end, matches the above?

# Using `sequence` for unequal lengths
print(sequence(c(4, 3, 4, 4, 4, 5)))            # (viii) Unequal sequences

# Generating repeats with `rep`
print(rep(9, 5))                                # (ix) Repeat 9 five times
print(rep(1:4, 2))                              # (ix) Repeat entire sequence twice
print(rep(1:4, each = 2))                       # (ix) Repeat each element twice
print(rep(1:4, each = 2, times = 3))            # (ix) Repeat each element twice, repeated three times
print(rep(1:4, 1:4))                            # (ix) Repeat each element based on its position
print(rep(1:4, c(4, 1, 4, 2)))                  # (x) Repeat elements based on specified counts
print(rep(c("cat", "dog", "goldfish", "rat"), c(2, 3, 2, 1)))  # (x) Repeat strings based on counts

# Using `seq` with varying step sizes
print(seq(-1, 1, by = 0.1))                     # (xi) Sequence from -1 to 1 with step size 0.1
print(seq(-1, 1, length = 7))                   # (xii) Sequence with 7 elements from -1 to 1

# Generate sequence without using `seq`
numbers <- (-1 + (0:20) * 0.1)                  # (xiii) Generate sequence from -1 to 1 with interval 0.1
print(numbers)                                  # Print the generated sequence


############### Q3: Missing values, infinity, and NaN handling ###############

# (i) 3 divided by 0
result1 <- 3 / 0  # Dividing a positive number by 0 gives Inf
print(result1)

# (ii) exp(-Inf)
result2 <- exp(-Inf)  # exp(-Inf) evaluates to 0
print(result2)

# (iii) exp(Inf)
result3 <- exp(Inf)  # exp(Inf) evaluates to Inf
print(result3)

# (iv) (0:3) raised to the power of Inf and c(0:3)**Inf
result4 <- (0:3) ** Inf  # This results in NaN
result5 <- c(0:3) ** Inf # This also results in NaN
print(result4)
print(result5)

# (v) 0 divided by 0
result6 <- 0 / 0  # Dividing 0 by 0 gives NaN
print(result6)

# (vi) Inf - Inf
result7 <- Inf - Inf  # Subtracting Inf from Inf gives NaN
print(result7)

# (vii) is.finite(10)
result8 <- is.finite(10)  # Checks if 10 is finite
print(result8)

# (viii) is.infinite(10)
result9 <- is.infinite(10)  # Checks if 10 is infinite
print(result9)

# (ix) y <- c(4, NA, 7)
y <- c(4, NA, 7)

# Checking if elements in y are "NA" (not working properly)
result10 <- y == "NA"  # Incorrect comparison as "NA" is not the same as NA
print(result10)

# (x) Checking if elements in y are NA
result11 <- is.na(y)  # Correct method to check for NA in vector y
print(result11)

# (xi) Strip y of the NA entry
result12 <- y[!is.na(y)]  # Removes NA entries from the vector y
print(result12)

# (xii) Creating a data frame and handling missing values
c1 <- c(1, 2, 3, NA)
c2 <- c(5, 6, NA, 8)
c3 <- c(9, NA, 11, 12)
c4 <- c(NA, 14, 15, 16)

# Creating a full data frame
full.frame <- data.frame(c1, c2, c3, c4)
print("Full Data Frame:")
print(full.frame)

# Removing rows with NA in c1
reduced.frame <- full.frame[!is.na(full.frame$c1),]
print("Reduced Data Frame:")
print(reduced.frame)

# Vector with missing values
x <- c(1, 2, 3, NA, 5)

# Calculate the mean with missing values (returns NA)
mean(x)  # This will give NA because of the NA present in the vector.

# Calculate the mean ignoring NA values
mean(x, na.rm = TRUE)  # This will return 2.75 as it ignores the NA value.


# (xiii) Generate a vector with NAs and identify their indices
v <- c(1:6, NA, NA, 9:12)
print("Vector with NAs:")
print(v)

# Get indices of NA using seq(along=v)
result13 <- seq(along=v)[is.na(v)]  # This shows the positions where NA is present
print(result13)

# OR Use the which function to get indices of NA
result14 <- which(is.na(v))  # This also gives the positions of NA
print(result14)

# Checking for NA presence directly in v
result15 <- is.na(v)  # Boolean check for NA in v
print(result15)


############### Q4 ###############

# (i) Vector assignment and queries
vec <- c(4, 7, 6, 5, 6, 7)
# Checking class and getting basic statistics
class(vec)  # Check the type of the vector
length(vec)  # Number of elements in the vector
min(vec)  # Minimum value in the vector
max(vec)  # Maximum value in the vector

# (ii) Vector creation using keyboard input
vec1 <- scan()  # Create vector using user input
# The user is expected to enter elements in the terminal/console and press Enter to finish the input.

# (iii) Vector subscripts (Indexing starts from 1)
vec[4]  # Access the 4th element of the vector

# (iv) Extracting multiple elements using two ways
ind <- c(2, 3, 6)  # Create index vector
vec[ind]  # Access elements at indices 2, 3, and 6
vec[c(2, 3, 6)]  # Same as above using a different notation

# (v) Drop elements using minus symbol (removes the 1st element from vec)
vec[-1]  # Drops the first element (access all except the 1st)

# (vi) What does vec[-length(vec)] do?
vec[-length(vec)]  # Drops the last element of vec

# (vii) Write a function called trim to remove the largest two and smallest two values
# from a vector called vec. Use this function on the above vec vector and print the result.

trim <- function(v) {
  sort(v)[-c(1, 2, (length(v) - 1), length(v))]  # Sorts and removes smallest two and largest two values
}

trim(vec)  # Call trim function with vec

# (viii) Use sequences to extract elements: vec[1:3], vec[seq(2,length(vec),2)]
vec[1:3]  # Get the first 3 elements
vec[seq(2, length(vec), 2)]  # Get every 2nd element, starting from the 2nd element

# Another way to extract every 2nd element, skipping odd indices
vec[-seq(1, length(vec), 2)]  # Removes all elements at odd indices (leaving even indices)
vec[1:length(vec) %% 2 == 0]  # Extract elements at even indices (without skipping any)

# (ix) Working with logical subscripts: Find the sum of values less than 5
x <- 0:10  # Create a vector of 0 through 10
x[x < 5]  # Pick elements that are less than 5
sum(x[x < 5])  # Sum of those elements

# Using subset() to filter values less than 5
sum(subset(x, x < 5))  # Filters values less than 5 and computes their sum

# Using 'which()' to get indices where elements are less than 5 and then summing
sum(x[which(x < 5)])  # Filters elements where values are less than 5 using 'which()'

# (x) Sum of the three largest values in a vector
x1 <- c(2, 64, 1, 43, 22, 33, 66)
sum(sort(x1, decreasing = TRUE)[1:3])  # Sort in descending order and sum the largest 3 values
# OR alternatively
sum(sort(x1)[(length(x1) - 2):length(x1)])  # Another approach by sorting in ascending and accessing the last 3 elements

# (xi) Finding index of vector corresponding to maximum or minimum value
which.max(x1)  # Index of the maximum value in x1
which.min(x1)  # Index of the minimum value in x1

# (xii) Combining vectors as columns or rows
cbind(1:10, 10:1)  # Combine two vectors as columns (must have same length)
rbind(1:10, 10:1)  # Combine two vectors as rows

# (xiii) Basic operations with vectors:
X <- c(1:10)  # Create vector X with values 1 through 10
X  # Show vector X

Y <- c(1:10 * 5)  # Create vector Y (elements are 1:10 multiplied by 5)
Y  # Show vector Y

X * Y  # Element-wise multiplication of X and Y
X + Y  # Element-wise addition of X and Y
X - Y  # Element-wise subtraction of X and Y
X / Y  # Element-wise division of X and Y
X ^ Y  # Element-wise exponentiation of X raised to powers of Y
log(X)  # Logarithm of each element in X
exp(Y)  # Exponent of each element in Y


############### Q5 ###############

# Matrices: 2D arrays where all elements are numbers.
# They have a fixed number of rows and columns.

# Arrays: Multidimensional structures that extend the matrix concept to more than 2 dimensions.
# Arrays can hold numbers and allow you to specify multiple dimensions (e.g., rows, columns, and depth).
# For example, in a 3D array with dimensions M[i, j, k]:
# - i refers to the number of rows.
# - j refers to the number of columns.
# - k refers to the number of 2D matrices (layers) stacked in the 3rd dimension.

# Dataframes: 2D structures that allow a mix of different data types across columns.
# For example, one column can contain numeric data, another can contain character strings, 
# and another can contain logical values.

# Lists: Flexible data structures that can contain mixed data types in different components.
# Lists can hold vectors, matrices, arrays, data frames, or even other lists.

# Example:
# Arrays (3D): Represented as M[i, j, k], where elements are accessed with three subscripts.
# Dataframes: Used for tabular data with heterogeneous column types.
# Matrices: A specialized 2D version of arrays, with only numeric or logical elements.
# Lists: Generalized container for diverse objects, including structured data like vectors or matrices.

# Exploring multi-dimensional objects in R: Matrices, Arrays, Dataframes

# (i) 
# Create a 3D array using a sequence from 1 to 24
y <- 1:24
dim(y) <- c(2, 4, 3) # 2 rows, 4 columns, 3 matrices
y  # Output the 3D array

# What does it look like with different dimensions?
dim(y) <- c(3, 2, 4) # Change to 3 rows, 2 columns, 4 matrices
y  # Output the transformed array

# Matrices: Creating 2D Arrays
# (ii) Using the `matrix()` function
X <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)
print(X) 

# Creating a named matrix 
X_named <- matrix(
  c(0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0),
  nrow = 3,
  dimnames = list(c("row1", "row2", "row3"), c("C1", "C2", "C3", "C4"))
)
X_named  # Matrix with named rows and columns

# Using `dim()` to convert vectors into matrices
# Create a vector and convert it into matrices
vector <- c(1, 2, 3, 4, 4, 3, 2, 1)

# Method 1: Using `matrix()` with `byrow = TRUE`
v <- matrix(vector,byrow = TRUE, nrow = 2)  # Fill row-wise
v

# Method 2: Default behavior with `byrow = FALSE`
v2 <- matrix(vector, nrow = 2)  # Fill column-wise (default)
v2

# Method 3: Changing dimensions using `dim()` directly
dim(vector) <- c(4, 2)  # Convert into a 4x2 matrix
vector

# Additional exploration
v3 <- matrix(vector, byrow = FALSE, nrow = 2)  # Filling column-wise explicitly
v3


#Checking if its a matrix
is.matrix(vector)


############### Q6 ###############

# (i) Common functions used with vectors  
vector <- c(1, 2, 3, 4, 5)   # Define vector  
min(vector)        # Finds the minimum value in the vector  
max(vector)        # Finds the maximum value in the vector  
sum(vector)        # Calculates the sum of all elements in the vector  
range(vector)      # Returns the minimum and maximum as a range  
sort(vector)       # Sorts the vector in ascending order  
colMeans(as.matrix(vector)) # Calculates the mean of the column  

# (ii) Outer product  
X 
Y
Z <- X[1:4] %o% Y[1:3]  # Outer product of X[1:4] and Y[1:3]  
Z   # Print the matrix Z  

# Reverse outer product  
YoX <- Y[1:3] %o% X[1:4]  
YoX   # Print the matrix YoX  

# (iii) Matrix multiplication examples  
t(Z)  # Transpose of matrix Z  
t(YoX)  # Transpose of matrix YoX  

# Create a square matrix from vector Y  
Y2 <- Y[1:9]              # Select first 9 elements of Y  
Y2 <- as.matrix(Y2)       # Convert to matrix  
Y2
dim(Y2) <- c(3, 3)        # Reshape into a 3x3 matrix  
Y2                        # Print matrix Y2  

# Dot product  
X %*% Y2                  # Dot product of X and Y2  
sum(X * Y2)               # Alternative way to compute dot product  

# Cross product  
crossprod(X[1:4], Z)      # Cross product of X[1:4] and Z  

# Identity matrix  
diag(4)   # Create a 4x4 identity matrix  

# (iv) Checking the class of variables  
class(X)   # Returns the class/type of X  

