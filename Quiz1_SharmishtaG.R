# QUIZ 

# Question 1``

# creating a function to pass input x,y
arithemtic_ops <- function(x,y){
  add <- x+y
  sub <- x-y
  mul <- x*y
  div <- x/y
  mod <- x%%y
  trunc <- x%/%y
  res <- c(add,sub,mul,div,mod,trunc) # created a single object with all the variables
  return(res) # returning the result
}
# passing input to the function 
resV <- arithemtic_ops(4,2) # capturing the result in a variable
paste("Addition: ", resV[1]," Subtraction: ",resV[2], " Multiplication: ",resV[3], " Division: ",resV[4]," Reminder: ",resV[5], " Integer Quotient: ",resV[6]) #concatenation by paste

#Question2

#creating a function to pass the coeffecients
RootsOfQE <- function(a,b,c){
  root1= (-b + sqrt(b^2-(4*a*c))) / 2*a #roots of QE formula
  root2= (-b - sqrt(b^2-(4*a*c))) / 2*a
  res2 <- c(root1,root2) #created a single object with all the variables (root1,root2)
  return(res2)  # returning the result
}
# passing input to the function 
R <- RootsOfQE(1,-7,12)
x1<-R[1] #first root
x2<-R[2] #second root
paste("Roots of the quadratic equation are: x1=",x,", x2=",x2) #concatenation by paste