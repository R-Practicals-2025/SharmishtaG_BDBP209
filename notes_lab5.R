data=read.csv("/home/ibab/R/Lab5/BrainCancer.csv",header=TRUE)
print(data)
print(dim(data))
print(length(data)) #gives no of columns
print(colnames(data)) #column names of data
print(rownames(data)) #row name if not present it will give serial number
#head gives- first 6 rows + headers
print(head(data,30))
print(table(data$diagnosis)) #gives frequency of each point in data
print(table(data$gtv,data$diagnosis)) #gives frequency of each point in data
print(data$diagnosis[1:5])
print(mean(data$gtv))
print(sd(data$gtv))

#attach data
attach(data)
print(mean(gtv))
print(sd(gtv))

print(summary(data$gtv))
hist(data$gtv) #default histogram
library(moments)
print(skewness)
