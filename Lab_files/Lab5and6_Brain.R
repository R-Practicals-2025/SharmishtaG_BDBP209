#Lab 5 jan 31 2025


data=read.csv("/Users/sharmishtaganesh/Desktop/Biostat_Lab/Lab_files/BrainCancer.csv", header=TRUE)
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

#vi)categorical 
print(colnames(data)) #"sex" "diagnosis" "loc" "stereo" 
print(head(data,3))

#vii)no. of categorical variables
print(unique(data$sex)) #"Female" "Male"  
print(unique(data$diagnosis)) #"Meningioma" "HG glioma"  "LG glioma"  NA  "Other"
print(unique(data$loc)) #"Infratentorial" "Supratentorial"
print(unique(data$stereo)) #"SRS" "SRT"

#viii)no. of levels in each categorical variable

print(length(unique(data$sex)))
print(length(unique(data$diagnosis)))
print(length(unique(data$loc)))
print(length(unique(data$stereo)))

print(class(data$sex))
data$sex <- factor(data$sex,levels=c("Male","Female"))
print(class(data$sex))
print(is.factor(data$sex))
print(nlevels(data$sex))
print(levels(data$sex))

#Ex3
#statistical parameters of data
#i) mean of GTV column
print(mean(data$gtv)) #8.660795, as mean is more than median.

#ii) mean of time column
print(mean(data$time))

#iii) median of GTV  #6.51
print(median(data$gtv))

#iv) mode of GTV
mode_GTV <- as.numeric(names(which.max(table(RestECG))))
print(mode_GTV)

# Interpretation:
# If mean > median → Right-skewed distribution
# If mean < median → Left-skewed distribution
# If mean ≈ median → Symmetric distribution

#v) standard deviation
print(sd(data$gtv))

#vi) stat summary
print(summary(data$gtv))

#vii) histogram plot
hist(data$gtv) 

#viii) skewness
install.packages("moments")
library(moments)
print(skewness(data$gtv))
# ➡️ Interpretation:
# Skewness > 0 → Right-skewed
# Skewness < 0 → Left-skewed
# Skewness = 0 → Symmetric

#ix) kurtsosis
print(kurtosis(data$gtv))
# Interpretation:
# Kurtosis > 3 → Leptokurtic (sharp peak)
# Kurtosis < 3 → Platykurtic (flatter peak)
# Kurtosis = 3 → Mesokurtic (normal distribution)

#x) 
par(mfrow = c(1, 3))
boxplot(range =0.1,data$gtv,xlabel="spread of GTV",ylabel="GTV",horizontal = FALSE,border=c("blue"),col=c("red"))
boxplot(range =0.2,data$gtv,xlabel="spread of GTV",ylabel="GTV",horizontal = FALSE,border=c("blue"),col=c("red"))
boxplot(range =0.05,data$gtv,xlabel="spread of GTV",ylabel="GTV",horizontal = FALSE,border=c("blue"),col=c("red"))
#range 0.1,0.2 shows fewer outliers
#range 0.05 shows more outliers
par(mfrow = c(1, 1))

#xi)
boxplot(data$time) # most data lie within the IQR
boxplot(data$ki) # 2 Outliers 1st and 2nd quartile almost coincide
boxplot(data$gtv) # more outliers # left skewed

#boxplot(data$time) has the broadest distribution
# Q1 ≈ Median:
# This indicates that 50% of the data is tightly packed in the lower range.
# The dataset likely has low variability in this region.
# Potential Distribution Shape:
# Left-skewed distribution (if the tail extends to the left).
# Heavily clustered data near the lower range.
# Key Observations in Boxplot:
# Q1 ≈ Median → Data is tightly packed in the lower range.
# Whiskers on one side are longer → Suggests skewness.
# Outliers appear as individual points beyond the whiskers.


#Ex4

#i)
filter1 = subset(data, data$gtv > 20)
print(filter1)
print(dim(filter1))
#ii)
filter2 = data[c(1,3,8,9,13,14,18,21),]
print(filter2)
#iii)
filter3_ind = which(data$sex=="Female")
filter3 = data[filter3_ind,]
print(filter3)
#iv)
data$new_col = (data$gtv*data$ki)/234
new_df = data[,c("gtv","ki","new_col")]
print(new_df)

data = data[,c("X",  "sex", "diagnosis" , "loc" , "ki"  , "gtv", "stereo" , "status"  ,"time" )]
print(data)

#v)
write.csv(filter3, file= "Lab5_female_BrainCancer.csv", row.names = FALSE)



### LAB 6 ###

#Ex 5
#(i)
print(class(data$sex))
data$sex<-factor(data$sex,levels=c("Male","Female"))
print(class(data$sex))
print(is.factor(data$sex))
#(ii)
print(nlevels(data$sex))
print(levels(data$sex))
#(iii)
print(unique(data$diagnosis))
data$diagnosis<-factor(data$diagnosis,levels=c("Meningioma","HG glioma","LG glioma","Other"))
print(levels(data$diagnosis))


#Ex 6
# adding new column with same number if riws as dataset
# Print the number of rows in the dataset
print(dim(data)[1])

# Generate the 'Temperature' category with the same number of rows as the dataset
num_rows <- dim(data)[1]  # Get number of rows
temp <- gl(n = 3, k = ceiling(num_rows / 3), length = num_rows, labels = c("Hot", "Cold", "Lukewarm"))

# Add the new 'Temperature' category to the dataframe
new_data <- data.frame(data, Temperature = temp)

# Print the new dataframe
print(new_data)


#Ex 7 

sub=subset(data, ki==70)
sub = sub[order(sub$gtv),]
print(mean(sub$gtv))
# chopping off tails
ind=floor(nrow(sub) * 0.1)
sub2 = sub[(ind + 1):(nrow(sub) - ind), ]
print(sub2)
print(mean(sub2$gtv))
# tapply
print(tapply(data$gtv, data$ki, mean))
print(tapply(data$gtv, data$ki, mean, trim=0.1))

#Ex 8
min=pmin(data$gtv,data$ki,data$time)
print(min)
max=pmax(data$gtv,data$ki,data$time)
print(max)

#Ex 9
# sorting and retrieved data
#(i)
ranks <- rank(data$gtv)
sorted<-sort(data$gtv)
ordered<-order(data$gtv)
view<-data.frame(data$gtv,ranks,sorted,ordered)
print(view)
#(ii)
diagnosis_sorted=data$diagnosis[ordered]
new_data <- data.frame(gtv = sorted, diagnosis = diagnosis_sorted)
write.csv(new_data, "lab6_ordered_data.csv", row.names = FALSE)


#Ex 10
#(i)
fil=data[1:6,3:8]
#(ii)
filter_mat=as.matrix(fil)
print(class(filter_mat))
print(mode(filter_mat))
print(attributes(filter_mat))

# - **`class()`** → Describes the **type of object** (e.g., `matrix`, `data.frame`, `list`).  
# - **`mode()`** → Describes the **data type** stored inside the object (e.g., `numeric`, `character`).  
# - **`attributes()`** → Provides **metadata** like dimensions (`dim`), column names (`names`), etc.  
# 
# **Key Difference:**  
#   - **`class()`** → Focuses on **structure/behavior**.  
# - **`mode()`** → Focuses on **content type**.  
# - **`attributes()`** → Focuses on **additional details** like labels or dimension info.

#(iii)
new_col=data$ki+data$gtv+data$time
#(iv)
new_data2=data.frame(data,new_col)
print(new_data2)
print(colnames(new_data2))
#(v)
new_data3=cbind(data,new_col)
print(colnames(new_data3))
#(vi)
#fil2=data[c(1,7,8,8),]
# Select rows 26 and 35 from the original data
fil2 = data[c(26, 35), ]

# Append the selected rows to the original data
new_data3 = rbind(data, fil2)

# Reset row indices (Fixes incorrect row numbering)
rownames(new_data3) <- NULL

# Print the updated dataframe
print(new_data3)

# Print the new dimensions (number of rows and columns)
print(dim(new_data3))


#Ex 11
#Defining row names and column names
x<-matrix(c(1,0,2,5,3,1,1,3,1,3,3,1,0,2,2,1,0,2,1,0),nrow=4)
print(x)
print(rownames(x))
print(colnames(x))
rownames(x)<-rownames(x,do.NULL = FALSE,prefix='trial.')
drugs<-c("aspirin","paracetamol","nurofen","hedex","placebo")
colnames(x)<-drugs
print(x)
dimnames(x) <- list(rownames(x), paste("drug", 1:5, sep=""))
print(x)
dimnames(x) <- list(NULL, paste("drug", 1:5, sep=""))
print(x)

#Ex 12
# apply for matrix
# Print the matrix
X<-matrix(c(1,0,2,5,3,1,1,3,1,3,3,1,0,2,2,1,0,2,1,0),nrow=4)
print("Matrix X:")
print(X)

# (i) Mean of the 5th column
print("Mean of 5th column:")
print(mean(X[,5]))

# (ii) Variance of the 4th row
print("Variance of 4th row:")
print(var(X[4,]))

# (iii) Row Sums
print("Row Sums:")
print(rowSums(X))

# Alternative method
print("Row Sums (apply function):")
print(apply(X, 1, sum))

# (iv) Column Sums
print("Column Sums:")
print(colSums(X))

# Alternative method
print("Column Sums (apply function):")
print(apply(X, 2, sum))

# (v) Row Means
print("Row Means:")
print(rowMeans(X))

# Alternative method
print("Row Means (apply function):")
print(apply(X, 1, mean))

# (vi) Column Means
print("Column Means:")
print(colMeans(X))

# Alternative method
print("Column Means (apply function):")
print(apply(X, 2, mean))

# (vii) Sum groups of rows within a column
# Define group labels for each row
group <- c("A", "B", "B", "A")

# Using rowsum()
print("Summing rows by group (rowsum function):")
print(rowsum(X, group))

# (b) Using row(X) and col(X)
print("Row indices:")
print(row(X))

print("Column indices:")
print(col(X))

# (c) Alternative method using tapply()
# NOT PREFERRED
print("Summing rows using tapply:")
print(tapply(X, list(group[row(X)], col(X)), sum))

# (d) Alternative method using aggregate()
print("Summing rows using aggregate:")
print(aggregate(X, list(group), sum))


# (viii)  Shuffling of elements
X2<-matrix(c(1,0,2,5,3,1,1,3,1,3,3,1,0,2,2,1,0,2,1,0),nrow=4)
print(X2)
Y<-apply(X2,2,sample) #column wise shuffling
print(Y)
Y1<-apply(X2,1,sample) #row wise shuffling
print(t(Y1))

# (ix)
X3 <- rbind(X2, apply(X2,2,mean))
print(X3)
X4 <- cbind(X3,apply(X3,1,var))
print(X4)
headings <- c(paste("drug.",1:5,sep=""),"var")
dimnames(X4)<- list(rownames(X4),headings)
print(X4)
headings2 <- c(paste("Trial-",1:4,sep=""),"Mean")
rownames(X4)<-headings2
print(X4)


# Ex 13
# (i) and (ii)
# sweep for minusing mean from orginal data
eg_sweep<-data.frame(data$ki, data$gtv, data$time)
cols<-apply(eg_sweep,2,mean)
print(cols)

# (iii)
# sweep from scratch
cols.means <- matrix(rep(cols,rep(dim(eg_sweep)[1],dim(eg_sweep)[2])),
                     nrow=dim(eg_sweep)[1])
print(cols.means)
eg_sweep_alt <- eg_sweep - cols.means
print("Method 1")
print(eg_sweep_alt)


# (iv)
eg_sweep1 <- sweep(eg_sweep,2,cols)
#columns will be departures from relevant col means
print("Method 2")
print(eg_sweep1)


# Ex 14
# Read the data from the file
data <- read.table("/Users/sharmishtaganesh/Desktop/Biostat_Lab/Lab_files/pgfull.txt", header = TRUE)

# Extract columns 1 to 54
species <- data[, 1:54]
print(species)

# Get the column index of the maximum value in each row
max_indices <- max.col(species)
print(max_indices)

# Get the names of the species corresponding to the max indices
species_names <- names(species)[max_indices]
print(species_names)
# Create a frequency table of the species
species_table <- table(species_names)
print(species_table)

# Get the column index of the minimum value in each row
min_indices <- max.col(-species)
print(min_indices)

# Get the names of the species corresponding to the min indices
min_species_names <- names(species)[min_indices]

# Create a frequency table of the species with minimum values
min_species_table <- table(min_species_names)
print(min_species_table)


# EX 15

# (i) Creating a list
apples <- c(4, 4.5, 4.2, 5.1, 3.9)
oranges <- c(TRUE, TRUE, FALSE)
chalk <- c("limestone", "marl", "oolite", "CaCO3")
cheese <- c(3.2 - 4.5i, 12.8 + 2.2i )
items <- list(apples, oranges, chalk, cheese)
print(items)

# Accessing elements
print(items[[3]])        # Access the third element (chalk)
print(items[[3]][3])     # Access the third element of the third list (the third item in chalk)

# Attempting to create a data frame
# df <- data.frame(apples, oranges, chalk, cheese)  # This will throw a warning/error

# Difference between items[3] and items[[3]]
print(items[3])   # Returns a list containing the chalk vector
print(items[[3]])  # Returns the chalk vector itself

# (ii) Creating a named list
items <- list(first = apples, second = oranges, third = chalk, fourth = cheese)
print(items$fourth)  # Access the 'fourth' element (cheese)
print(class(items))   # Check the class of the list

# (iii) Using lapply()
print(lapply(items, length))  # Returns the length of each element in the list
print(lapply(items, class))    # Returns the class of each element in the list
print(lapply(items, mean))     # Will return a list with the mean of apples and NA for others

# (iv) Summary and structure functions
print(summary(items))  # Summary of the list
print(str(items))      # Structure of the list

#str() gives a structural overview, while summary() provides statistical insights.
#Class indicates how R treats an object, mode indicates the data type, and attributes provide extra information about the object.


