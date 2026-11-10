#Download and install/update Bioconductor packages
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.20")
#BiocManager::install (c('limma','affy','affyPLM','edgeR','igraph', 'ggplot2', 'reshape2', 'gplots'))

#Download and install/update R packages from CRAN 
install.packages("caret")
#Load an already installed package
library(caret)

#-------------------------------------------------------------------------------------------------------------

#Online help
help("dir.create")		
?"dir.create" #or ?? or F1


#-------------------------------------------------------------------------------------------------------------

#Setting the current working directory
setwd("C:/Users/39339/Desktop/BCG_master_degree/bioinformatics_and_computational_biology/r_bioinfo_esercitazioni")


#-------------------------------------------------------------------------------------------------------------

#Data types

#Vectors 

x <- c('YES', 'YES', 'NO', 'MAYBE', 'NO')
is.character(x)
x
x_f <- factor(x) #or x_f <- as_factor(x)
x_f
is.character(x_f)
is.factor(x_f)
length(x) # unary operator
length(x_f) #always the length of the vector
levels(x_f)
is.character(levels(x_f))

v <- vector(mode = 'character',length = 10)
v[1]<-'hi'
v
v <- v[-1] 
length(v)
v

v1 <- c (1,2,3)
v1[1:2] #Extract first two elements without any change for v1
v1
v[1:3]<-v1 
v #automatic conversion into character values


v2 <- seq(1, 10, 3)
v2

#1. Replicate v1 for 5 times, repeating each element twice (2 times)
?rep
#2. Assign a vector with 10 random numbers, normally distributed with mean=5 and standard deviation=2
?rnorm
#3. Extract from the second to the to last but one element of v4
#?

#Binary operators
vsum <- v1+v2
v1
v2
vsum 

vsum==v2 #example of a logical operator to obtain a predicate
#control functions
any(v1 > 2)
all(v1 > 2)
ifelse(v1 > 2, 2, 0) 
ifelse(!v1 > 2, 2, 0)

v4 <- c(1,6,4,8,-7,2,3,NA) 
#4. Extract v4 values having their square greater than 18 


#5. Extract v4 values having their square greater than 18 excluding NA
?subset


#6. return their indexes on v4
?which


#7. replace such values with 100 in v4



#8. delete from v4 the last 3 elements 



#-------------------------------------------------------------------------------------------------------------

#Matrices

m1 <- cbind(c(1,2,3),1)
m1

m2 <- rbind(c(1,2,3,4),c(5,6,7,8))
m2

m3 <-  rbind(m2, c(2,4,6))
m3 <-  rbind(m2, c(2,3))
m3

m3[1,3]
m3[,3]
m3[1:2,]
m3[1:3,1:3]

m3[,c(1,3)]
#9. replace all the values in the first and third column of m3 with 10


#-------------------------------------------------------------------------------------------------------------

#Lists 

x <- c(2,5,7)
y <- c("a","b","c")
z <- 3
L <- list(x,y,z)
L
names(L)
L1 <- list(A = c(2,4,6), B = c("a","b","c"))
names(L1)
str(L1)

L1$B #vector - chr
L1[["B"]] #vector - chr
L1[[2]] #vector - chr

L1["B"] #sublist made of a vector - chr
is.list(L1["B"])
L1["B"][2] #NULL

#10. extract the second element contained within the B element of the list L1 


L2 <- L1["B"]
str(L2)

x <- c(2,5,7,8)
y <- c("a","b","c")
z <- 3
L3 <- list(x,y,z)
str(L3)
L3
L
L1

#-------------------------------------------------------------------------------------------------------------

#Data frames

D <- data.frame(L3) # unmatching dimensions
D <- data.frame(L) # not named columns
D

D <- data.frame(L1)
D
rownames(D)<- c('first','second','third')
D
summary(D)
D[1,1]
D['first', 'B']
D[,2]
D$B #only for columns
D['first',] # NOT D$first!


#-------------------------------------------------------------------------------------------------------------

#R basics - Data frames working examples
# Load the tab-delimited metadata.txt file and use indexing and predicates to extract required information

?read.table
D <- read.table("metadata.txt", header=T, sep="\t")
typeof(D)
D<-data.frame(D)
str(D)
head(D)

#11. Example 1
# Find all the histone modification samples (HM)
# from the Broad Institute




#12. Example 2
# Find all the filenames of samples coming from
# RAD21 transcription factor in females



#13. Example 3
# Print a table (use table() function) of cell lines (cell column) with their respective 
# number of experiments not involving H3K27me3 histone modification (antibody column)
?table



#14. Example 4
# Order metadata by cell line name (use order() function)

#order returns simply the ordered list of indexes

#-------------------------------------------------------------------------------------------------------------

#R basics - Strings
# Create vectors from strings

strsplit("1;2;3", split = ";") #returns a list

strsplit("1 2 3", " ")[[1]] #returns a character vector

as.numeric(strsplit("1 2 3", " ")[[1]]) #returns a numeric vector

ls()
rm(list=ls())

#15. Strings: a simple exercise
#Change file extension from .tab to .txt using strsplit and paste functions

file_path <- "/home/user1/documents/test.tab"
# \\. Escape sequence to avoid problems with dot pattern in using strsplit function
?strsplit
?paste




#-------------------------------------------------------------------------------------------------------------
#Simple Plots

# Use data() to list sample datasets in R.
# We will use the Puromycin data set.
# This set is stored as a data frame with
# three attributes: conc, rate and state.
# Attribute conc is the concentration of
# substrate, rate is the  speed of the enzimatic 
# reaction, and state is the state of a sample,
# (i.e. untreated or treated with Puromycin)
data()
head(Puromycin)   # have a look at the data frame

summary(Puromycin)    # basic data frame stats

# Analyze the dataset making loess regressions
# and scatter plots for untreated and treated
# samples

# Split the dataset
D1 <- Puromycin[Puromycin$state=="untreated",]
D2 <- Puromycin[Puromycin$state=="treated",]

# example of regression models
?loess
loess1 <- loess(D1$rate ~ D1$conc)
loess2 <- loess(D2$rate ~ D2$conc)
# Model summaries
summary(loess1)
summary(loess2)

# Scatter plots
# open a device for pdf file creation
pdf("plot.pdf")
# divide in 2 rows and 1 column
help(par)
par(mfrow=c(2,1))      

plot(D1$conc, D1$rate, xlab='concentration', ylab='rate', main='scatter plot for untreated state')
plot(D2$conc, D2$rate, xlab='concentration', ylab='rate', main='scatter plot for treated state')

# close the device
dev.off()              


#-------------------------------------------------------------------------------------------------------------

#R basics - Functions

#First example
f2 <-function(arg){
 res<-arg*2
}
f2(5) #does not return anything
a<-f2(5) #a is associated with the value of the last line of the body of the function 
a

f2 <-function(arg){
  res<-arg*2
  return(res)
}
f2(6) #prints the returned value
b<-f2(6) #b is associated with the returned value
b

f2 <-function(arg){
  res<-arg*2
  res # 5 .... be careful
}
f2(7) #prints the value of the last line of the body of the function
c<-f2(7) #c is associated of the last line of the body of the function 
c

#Another example
f <-function(arg1, arg2)
{  
  sum<-arg1+f2(arg2)  
  #k<-10
  k<<-10 #super-assignment to save the variable in the global environment
  return(sum*k)
  #return(8)
}

final<-f(1,2)
final

#Iterations
#apply
M <- rbind(c(3,6,4), c(2,9,8), c(7,5,1))
M
D_M <- data.frame(M)
m_M<-apply(M, 2, mean)
m_D_M<-apply(D_M, 2, mean)
m_M #vector
m_D_M #named vector

#a single-line-body function can be defined directly within apply 
apply(M, 2, function(x) (x*x)^x ) # median(c(3^3,6^6,4^4)),...
apply(M, 2, function(x) all(x >= 3))

#-------------------------------------------------------------------------------------------------------------

#R scripting - EXERCISE_1 - EXPRESSION DATA FINAL NORMALIZATION

# Write an R function to standardize a numeric matrix X.
# The standardized version of X is:
# Z = (X - M) / S
# where M is the mean matrix of X and S the standard deviation matrix of X. 

X <- rbind(c(2,9,8,0,5,2),
           c(3,8,0,1,9,2),
           c(4,3,6,2,8,3))
X
#strategy
apply(X, 2, mean)
U <- rep(1, nrow(X)) # a vector with nrow(X) positions
matrix(rep(1, nrow(X))) # one-column matrix with nrow(X) rows (3,1)
matrix(apply(X, 2, mean)) # one-column matrix with as many rows as the columns of X
U %*% t(apply(X, 2, mean)) #implicit casting to matrices: (3,1)%*%(1,6) = ... 

#function
standardizeData <- function(X) {
  U <- rep(1, nrow(X))
  M <- U %*% t(apply(X, 2, mean)) #(3,1)%*%(1,6) = ...
  S <- U %*% t(apply(X, 2, sd)) #as.matrix(U)
  Z <- (X - M)/S
  return(Z)
}

X_norm=standardizeData(X)
X_norm

#-------------------------------------------------------------------------------------------------------------

#R scripting - EXERCISE_2 - REDUNDANCY REMOVAL

# Write an R function to remove one of the two variables (columns) from an input matrix X 
# with Pearson correlation coefficient greater or equal to a certain threshold value r (e.g. 1).

r<-1
X <- rbind(c(2,9,8,0,5,2),
           c(3,8,0,1,9,2),
           c(4,3,6,2,8,3))
#strategy
R <- cor(X)
R
R[upper.tri(R)] <- 0
diag(R) <- 0
apply(R, 2, function(x) any(x >= r))
!apply(R, 2, function(x) any(x >= r))
X[,-apply(R, 2, function(x) any(x >= r))] #! or - is the same here

#function
removeCorr <- function(X, r=1) {
  R <- cor(X)
  R[upper.tri(R)] <- 0
  diag(R) <- 0
  NEW <- X[, !apply(R, 2, function(x) any(x >= r))] #remove the columns having a value of at least r 
  return(NEW)
}

Y<- removeCorr(X, r=0.9)
Y

#-------------------------------------------------------------------------------------------------------------

#R scripting - FIRST SCRIPT

# Standardize and filter out correlated vars

r <- max(cor(X))
R <- removeCorr(X, r)
Z <- standardizeData(R)

# Exporting to file

Z <- data.frame(Z)
Z
write.table(Z, file = "result_first_script.txt", quote = FALSE,sep = "\t")

