#Lab1 1 Jan 3 2025

#Q1
#Q1.1
print(2.7/2)
#Q1.2
print(2.7%/%2)
#Q1.3
print(10+5i/2)
#Q1.4
print(round(2.5))
#Q1.5
print(round(-2.5))
#Q1.6
print(2%/%4-1)
#Q1.7
print(3*2**2)
#Q1.8
print(3**2*2)
#Q1.9
print(7%/%4)
#Q1.10
print(7%%4)
# Q1.11
print(-7%%4)
# Q1.12
trunc(5.7)
# Q1.13
trunc(-5.7)
#signif() in R
#Rounds a number to n significant digits.
#This method retains meaningful precision across large or small values.

#round() in R
#Rounds a number to n decimal places.
#The precision is controlled by the number of digits after the decimal point.

# Q2: Implementing ceiling 

# Adds 1 instead of 0.5; doesn't replicate the ceiling function correctly.
ceilingNum <- function(x){floor(x + 1)}
ceilingNum(5.7)  # Output: 6 
print(paste('Alternate ceiling function result: ',ceil_alt(5.7)),quote=FALSE)

# Q3: Logical operations, mistakes could arise due to misunderstanding the operators

a <- 1
b <- 2
c <- 4

# Q3.1 Logical AND with incorrect syntax: a && b is technically correct,
# but for numeric values, logical operations consider "non-zero" as TRUE.
a && b  # Output: TRUE (this works, as 1 and 2 are treated as TRUE)

# Q3.2 Logical NOT with OR operation:
!(a < b) | (c > b)  # Output: TRUE (this works)


# Q4: Incorrect vector type conversion

# Q4.1 Vector `x` is a numeric vector (double by default).
x <- c(5, 3, 7, 8)

# Q4.2 Check if `x` is integer. This returns FALSE because it's numeric (double).
is.integer(x)  # Output: FALSE

# Q4.3 Check if `x` is numeric. This should return TRUE.
is.numeric(x)  # Output: TRUE

# Q4.4 Incorrect way of converting to integers:
x <- integer(x)  # Error: `x` is not the length of the vector; numeric values are passed.

# Expected error: "Error in integer(x) : invalid 'length' argument"

# Q4.5 Correct way would be to explicitly convert to integer:
x <- as.integer(x)

# Check if `x` is now an integer type. Should return TRUE after correction.
is.integer(x)  # Output: TRUE
is.numeric(x)  # Output: TRUE


# Q5
# Q5.1 Assign the square root of 2
x <- sqrt(2)

# Q5.2 Test equality with '=='
x * x == 2  # Likely FALSE due to precision errors

# Q5.3 Show the difference
x * x - 2  # A small value close to zero, e.g., 4.44e-16

# Q5.4 Use all.equal for proper comparison
all.equal(x * x, 2)  # Output: TRUE




