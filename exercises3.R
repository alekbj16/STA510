#File contains weekly R-code examples from lectures as well as weekly exercises

#Remove old variables
rm(list=ls())

#Check that psuedo random number generators are deterministic
runif(10) #10 random uniform between 0 and 1

#Set seed and runif again:
set.seed(98925)  #Lets us "choose" where to start generating random
runif(10) 

#Takeaway: When writing programs that invovles
# random numbers, it's a good idea set a fixed seed
# at the start of the program. In that way one 
# produces reproducable results every time



#Inverse transform examples

#Function that generates Nsim random variables
# with the exponential distribution
# using rate parameter lambda.

genexpdistr <- function(Nsim,lambda) {
  U <- runif(Nsim) #Nsim uniform random numbers
  X <- -log(1-U)/lambda #Inverse transfrom of F(x)
  return(X)
}

Nsim <- 10000
lambda <- 2
x <- genexpdistr(Nsim = Nsim, lambda = lambda)

#Histogram of simulated data with true density on top as red line
hist(x,prob = TRUE, breaks = seq(0,ceiling(max(x)*1.1),length.out = max(10,sqrt(Nsim))), main =  "Histogram of data and true exponential density")
curve(dexp(x, rate = lambda), col = "red", lty = 2, lwd = 2, add = TRUE,xlim=c(0.001,max(x)*1.1))

#Expectation (theoretical mean)
1/lambda 
#Mean of the data 
mean(x)
#Standard deviation of the data
sd(x)


#Inverse transform method, discrete case
# Use sample.function

x <- 0:3 #possible outcomes
n <- 1000 #number of draws
probs <- c(1,3,3,1)/8 #probability of outcome
X <- sample(x=x,size=n,replace=TRUE,prob=probs)


     

#Illustration of Accept-reject algorithm

#a function of the density
f <- function(x){
  return(0.75*(1.0-x^2))
}
C <- 1.5 #Found by calculations by hand

#For plotting purposes
x.grid <- seq(from=-1.0,to=1.0,length.out = 1000)
#Plot the target density f(x)
plot(x.grid,f(x.grid),type="l",lwd=2,xlab="x",ylab="")
#Plot C*g(x)
lines(x.grid,C*dunif(x.grid,min=-1,max=1),col=2,lwd=2)
#We should have f(x) <= C*g(x) for all possible x
#This looks good for now

#Then we actually implement the algorithm
#Single draw
gen.f.AR.draw <- function() {
  proposals <- 0 
  while(TRUE) { #this is scary, should be safeguarded
    x <- runif(1,min=-1,max=1) #proposal from uniform(-1,1)
    a <- f(x)/(C*dunif(x,min=-1,max=1)) #accept prob
    proposals <- proposals + 1 #count the proposals
    if(runif(1) < a) return(c(x,proposals))
  }
}

#Many draws
rf.AR <- function(n = 1000) {
  ret <- vector(mode = "numeric",length = n)
  nprop <- vector(mode = "numeric",length = n)
  for(i in 1:n) {
    tmp <- gen.f.AR.draw()
    ret[i] <- tmp[1]
    nprop[i] <- tmp[2]
  }
  print(paste0(
    "number of proposal per accepted RV: ",
    mean(nprop)),quote = F)
  return(ret)
}

#Draw some random numbers and add to figure
x <- rf.AR(n=10000)
hist(x,probability = T, add = T)



#Remove old variables
rm(list=ls())

#Exercises from the book
#3.3
#Pareto distribution
#cdf: F(x) = 1-(b/x)^a, x>b>0,a>0

#Generating pareto distribution using inverse tranform
genparetodistr <- function(Nsim,a,b){
  U <- runif(Nsim,min=0,max=1) #U = uniform[0,1]
  X <- b/sqrt((1-U)) #F(X)=U solved for X
  return(X)
}

Nsim <- 1
a <- 2
b <- 2
genparetodistr(Nsim = Nsim, a = a, b = b)

#Graph the density historgram with
#the Pareto(2,2) density superimposed

Nsim <- 10000
a <- 2
b <- 2
x <- genparetodistr(Nsim = Nsim, a = a, b = b)

#Histogram of simulated data with true density on top as red line
hist(x,prob = TRUE, breaks = seq(0,ceiling(max(x)*1.1),length.out = max(10,sqrt(Nsim))), main =  "Histogram of data and true pareto density")


#Remove old variables
rm(list=ls())

#3.5
#Discrete random variable 
x <- 0:4 #possible outcomes
probs <- c(0.1,0.2,0.2,0.2,0.3) #pmf
Nsim <- 1000
X <- sample(x=x,size=Nsim,replace=TRUE,prob=probs)
relfreq <- table(X)/Nsim
barplot(relfreq,ylab = "Relative frequency", main = paste(Nsim," simulations"))


# 3.7 
# Write a function to generate 
# a random sample of size n from
# the Beta(a,b) distribution by
# the acceptance-rejection method
# Generate a random sample of size 
# 1000 from the Beta(3,2) dist. 
# Graph the histogram of the sample
# with the theoretical Beta(3,2)
# density superimposed.

#Using relation stated in Rizzo, p. 72
Nsim <- 1000
a <- 3
b <- 2

betagengamma <- function(Nsim, a, b) {
  U <- rgamma(n = Nsim, shape = a, rate = 1)
  V <- rgamma(n = Nsim, shape = b, rate = 1)
  X <- U/(U+V)
  return(X)
}

# Compare with the Beta(3,2) dist using
# a quantile-quantile (QQ) plot
q <- qbeta(ppoints(Nsim),a,b)
qqplot(q, X, cex = 0.25, xlab = "Beta(3, 2)",ylab="Sample")
abline(0,1)




#3.17
#Compare runtime of X from 3.7 with rbeta
Nsim <- 5000
a <- 3
b <- 2
iterations <- 1000

start.time <- Sys.time()
runtimeapprox <- function(Nsim,a,b,iterations) {
  for(i in 1:iterations){
    values <- betagengamma(Nsim = Nsim, a = a, b = b)
  }
  return(values)
}
end.time <- Sys.time()
time.taken <-  end.time - start.time
time.taken


start.time <- Sys.time()
runtimeanalytic <-function(Nsim,a,b,iterations){
  for(i in 1:iterations) {
    values <- rbeta(Nsim, a, b)
  }
  return(values)
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken




#Exercise 1 from paper
Nsim <- 10000
p <- c(0.7,0.95,0.95,0.95,0.99,0.99,0.92,0.92,0.92,0.7)

#a
I1 <- sample(0:1,size = Nsim,replace=TRUE,prob=c(1-p[1],p[1]))
I2 <- sample(0:1,size = Nsim,replace=TRUE,prob=c(1-p[2],p[2]))
sum(I1*I2)/Nsim #The fraction of cases, when both have 1

#b
I1 <- sample(0:1,size=Nsim,replace=TRUE,prob=c(1-p[1],p[1]))
I2 <- sample(0:1,size=Nsim,replace=TRUE,prob=c(1-p[2],p[2]))
sum(pmax(I1,I2))/Nsim

#c
I1 <- sample(0:1,size=Nsim,replace=TRUE,prob=c(1-p[1],p[1]))
I2 <- sample(0:1,size=Nsim,replace=TRUE,prob=c(1-p[2],p[2]))
I3 <- sample(0:1,size=Nsim,replace=TRUE,prob=c(1-p[3],p[3]))
I4 <- sample(0:1,size=Nsim,replace=TRUE,prob=c(1-p[4],p[4]))
sum(I1*I2*pmax(I3,I4))/Nsim

#d
# "The system works if only at least 7 of 10
# components work"

#Create a matrix
Imatrix <- matrix(nrow = Nsim, ncol = 10)
for(i in 1:10) {
  Imatrix[,i] <- sample(0:1,size=Nsim,replace=TRUE,prob=c(1-p[i],p[i]))
}
Nfunctioning <- rowSums(Imatrix) # Calculate number of working for each row
sum(Nfunctioning>=7)/Nsim # Probability for at least 7 working

barplot(table(Nfunctioning)/Nsim)



#Exercise 2

#Remove old variables
rm(list=ls())

#First create a function to simulate values from a 
#general triangle distribution using acceptance-rejection:
rtriangle <- function(Nsim,a,b,c) {
  if(!(a<c & c<b))
    stop("Error, check a, b and c.")
  k <- 0 #counter for accepted
  x <- numeric(Nsim) #vector of accepted
  while(k<Nsim){ #iterate until k have been accepted Nsim times
    u <- runif(1) #generate U
    y <- runif(1,min=a,max=b) #proposed distribution
  }
    
  if(y<c){ #different pdf before and after c
      if(u<=((y-a)/(c-a))) { #then accept
        k <- k + 1 
        x[k] <- y #add proposal to x
      }
    }
  else{
    if(u<=((b-y)/(b-c))) { #then accept
      k <- k + 1 
      x[k] <- y #add proposal to x
    }
  }
  return(x)
}


#a: P(X>8)
#Conditional prob ( condition, we know there is a find)
Nsim <- 1000
probOil <- 0.4
a <- 2
c <- 6
b <- 10
xtriangle <- rtriangle(Nsim = Nsim,a = a, b = b, c = c)
print(sum(as.integer(xtriangle>8.0))/Nsim) #relative frequency

#Unconditional probability that there is a find
finds <- sample(x = c(0,1), size = Nsim, replace = TRUE, prob = c(0.6,0.4)) #0 for no find, 1 for find
unconditionalProb = finds*xtriangle
print(sum(as.integer(unconditionalProb>8.0))/n)

#b
Nsim <- 1000
X1 <- rtriangle(Nsim = Nsim, a = 2, b = 6, c = 4)
X2 <- rtriangle(Nsim = Nsim, a = 3, b = 11, c = 7)
X3 <- rtriangle(Nsim = Nsim, a = 2, b = 6, c = 4)
X4 <- rtriangle(Nsim = Nsim, a = 1, b = 9, c = 5)
X5 <- rtriangle(Nsim = Nsim, a = 8, b = 10, c = 9)
X6 <- rtriangle(Nsim = Nsim, a = 5, b = 9, c = 7)
X7 <- rtriangle(Nsim = Nsim, a = 2, b = 6, c = 4)
X8 <- rtriangle(Nsim = Nsim, a = 3, b = 5, c = 4)
X9 <- rtriangle(Nsim = Nsim, a = 8, b = 12, c = 10)
X10 <- rtriangle(Nsim = Nsim, a = 3, b = 7, c = 5)

#Assuming finding at every reservioir
totalFindings <- X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10
hist(totalFindings, prob = TRUE, density = 10, sqrt(Nsim))
mean(totalFindings)


#Findings with risk
pi <- c(0.8,0.3,0.6,0.6,0.5,0.9,0.5,0.8,0.4,0.4)
nonRiskMatrix <- cbind(X1,X2,X3,X4,X5,X7,X8,X9,X10)
riskMatrix <- nonRiskMatrix
for(i in 1:10) {
  riskMatrix[,i] <- riskMatrix[,i]*sample(x=c(0,1),size=Nsim,prob=c((1-pi[i]),pi[i]),replace = TRUE)
}
risk <- rowSums(riskMatrix)
hist(risk,prob=TRUE,density=10,class = sqrt(Nsim))
mean(risk)




###Problem 3
#Covariance matrix calculated on paper
library(mvtnorm)
mu <- c(1,1) #Expectatons
var <- c(1,3) #Variances
covarMatrix <- matrix(c(1,1.21,1.21,1),nrow = 2) #Create a covariance Matrix
Nsim <- 10000
Xdata <- rmvnorm(n = Nsim, mean = mu, sigma = covarMatrix)

#Probability that the sum is more that 3
probability <- sum((Xdata[,1]+Xdata[,2])>3)/Nsim #divide by Nsim to get relfreq

#b, covariance is 0
covarMatrix <- matrix(c(1,0,0,1),nrow = 2) #Create a covariance Matrix
Xdata <- rmvnorm(n = Nsim, mean = mu, sigma = covarMatrix)

#Probability that the sum is more that 3
probability <- sum((Xdata[,1]+Xdata[,2])>3)/Nsim #divide by Nsim to get relfreq

