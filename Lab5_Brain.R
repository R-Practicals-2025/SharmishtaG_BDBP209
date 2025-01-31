#Lab 5 jan 31 2025


data=read.csv("/home/ibab/R/Lab5/Heart.csv")
print(data)

# 2
# (i) Dimensions of the Data (Number of Rows and Columns)
print(dim(data))  
# This will print a vector with the number of rows and columns in the dataset.

# (ii) Column Names of the Data
print(colnames(data))  
# This will print the names of all columns in the dataset.

# (iii) Row Names of the Data
print(rownames(data))  
# This will print the row names. If not explicitly set, it will show the default row indices.

# (iv) First 30 Rows of the Data
print(head(data, 30))  
# This will display the first 30 rows of the dataset.

# (v) Frequency Table of the Data (for the entire dataset)
# Note: You can't directly use `table()` on a dataframe with mixed data types, so we need to apply it to individual columns
print(table(data))  
# This will give a frequency count for all categorical variables in the dataframe. 
# It’s better to apply this to individual columns for clearer output:
print(table(data$sex))  # Frequency table for 'sex'
print(table(data$diagnosis))  # Frequency table for 'diagnosis'

# (vi) How Many Categorical Variables are Present
# Identify categorical variables by checking column types

#X <- sapply(data, is.numeric)   #...gives numeric
#print(X)  
#y <- data[!X]    #....gives rows containing categorical variables as it removes numeric
#print(colnames(y))   # to print column names...

print(unique(data$ChestPain))
print(unique(data$Thal))
print(unique(data$AHD))
data$ChestPain <- factor(data$ChestPain, levels = c( "typical","asymptomatic", "nonanginal", "nontypical"  ))
data$Thal <- factor(data$Thal, levels = c("fixed","normal","reversable"))
data$AHD <- factor(data$AHD, levels = c("No","Yes"))
categorical_vars <- sapply(data, is.factor)  # Find columns that are factors
print(sum(categorical_vars))  # Sum of TRUE values tells us how many categorical variables are present

# (vii) What Are the Categorical Variables Present
# Print the names of columns that are categorical
print(colnames(data)[categorical_vars])  

# (viii) How Many Levels Does Each Categorical Variable Have and What Are the Levels
# Loop through categorical variables to get the number of levels and the levels themselves
for (var in colnames(data)[categorical_vars]) {
  print(paste("Variable:", var))
  print(paste("Number of levels:", nlevels(data[[var]])))
  print(paste("Levels:", levels(data[[var]])))
}


# 3
print(mean(data$RestECG))
print(mean(data$Chol))
print(median(data$RestECG))
print(mode(data$RestECG))
print(which.max(data$RestECG))
print(sd(data$RestECG))
print(summary(data$RestECG))
hist(data$RestBP)
library(moments)
print(skewness(data$RestECG))
print(kurtosis(data$RestECG))
boxplot(data$RestECG)
boxplot(data$RestECG,xlabel="spread of gtv",ylabel="GTV",horizontal=TRUE,border=c("blue"),col=c
        ("yellow"))

boxplot(data$RestECG,range=0.1,xlabel="spread of gtv",ylabel="GTV",horizontal=FALSE,border=c("blue"),col=c
        ("yellow"))
boxplot(data$RestECG,range=0.2,xlabel="spread of gtv",ylabel="GTV",horizontal=FALSE,border=c("blue"),col=c
        ("yellow"))
boxplot(data$RestECG,range=0.05,xlabel="spread of gtv",ylabel="GTV",horizontal=FALSE,border=c("blue"),col=c
        ("yellow"))
boxplot(data$Chol)
boxplot(data$Fbs)

# 4
filter1=subset(data,data$chol>20)
print(filter1)
print(dim(filter1))
filter2=unique(subset(data,data$chestpain=="typical"))
print(filter2)
filter3=data[c(1,3,8,9,13,14,18,21),]
print(filter3)
filter4_ind=which(data$sex=="0")
print(filter4_ind)