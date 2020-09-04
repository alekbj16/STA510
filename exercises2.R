#File contains weekly R-code examples from lectures as well as weekly exercises

#Remove old variables
rm(list=ls())


#Throw dice example from lecture:
#Function to simulate the sum of k throws of a dice
sumdice <- function(Nsim,k){ #Returns a vector. Each value of the vector is the sum of a dice thrown k times. The length of the vector is the number of times the simulation is repeated
  dsum <- vector(length=Nsim) 
  for(i in 1:Nsim)              
    dsum[i] <- sum(sample(1:6,size=k,replace=TRUE)) 
  return(dsum)                
}

sumdice(Nsim=20,k=2)
sumdice(Nsim=20,k=50)
#Average
sumdice(Nsim=20,k=50)/50


#Function to plot histogram of the sum of k throws repeated Nsim times
plotdicesums <- function(Nsim = 10000, k) {
  dicesums <- sumdice(Nsim = Nsim, k = k)
  relfreq <- table(dicesums)/Nsim   # Calculate relative frequency for each outcome
  barplot(relfreq,ylab="Relative frequency", main = paste("Sum of k =",k,"throws"))
}

#Repeat for different values of k
par(mfrow=c(2,3)) #Make a 2x3 display of plots
plotdicesums(Nsim=100000,k=1)
plotdicesums(Nsim=100000,k=2)
plotdicesums(Nsim=100000,k=3)
plotdicesums(Nsim=100000,k=5)
plotdicesums(Nsim=100000,k=10)
plotdicesums(Nsim=100000,k=50)

#Function to plot the histogram of the average of k throws repeated Nsim times
plotdiceaverages <- function(Nsim=10000,k){
  diceaverages <- sumdice(Nsim=Nsim,k=k)/k   # This gives us averages
  hist(diceaverages, prob = TRUE,breaks=seq(0.5,6.5,length.out=6*k+1),
       xlab="Average",main=paste("Average of k =",k,"throws"))
}

#Repeat for different values of k
par(mfrow=c(2,3)) #Make a 2x3 display of plots
plotdiceaverages(Nsim=100000,k=1)
plotdiceaverages(Nsim=100000,k=2)
plotdiceaverages(Nsim=100000,k=3)
plotdiceaverages(Nsim=100000,k=5)
plotdiceaverages(Nsim=100000,k=10)
plotdiceaverages(Nsim=100000,k=5000)
par(mfrow=c(1,1))  # Set back to a 1x1 display



#Simulate the probability of a sum larger than 200 when k=50:
nsim <- 1000000
dsums <- sumdice(Nsim=nsim,k=50)
sum(dsums>200)/nsim
#Calculating the probability based on SGT (see lecture notes):
1-pnorm(200,mean=50*3.5,sd=sqrt(50*2.92))




#Robot assemblyline lifetime
#Simulate times to failure for the system
systemfailuretimes <-function(Nsim) { #R is for robot, A for all other components
  R1 <- rgamma(Nsim, shape = 16, scale = 1)
  R2 <- rgamma(Nsim, shape = 16, scale = 1)
  R3 <- rgamma(Nsim, shape = 12, scale = 2)
  A  <- rlnorm(Nsim, meanlog = 2.3, sdlog = 1.1)
  failuretimes <- pmin(R1,R2,R3,A) #Nsim times to failure
  return(failuretimes)
}

#Fist estimates of mean and sd
ftimes <- systemfailuretimes(Nsim = 1000)
mu1 <- mean(ftimes)
sd1 <- sd(ftimes)
mu1
sd1

#Mean with precision of one month
e1 <- 1 #Timescale = 1 month, corresponds to precision r in lecture  notes
n1 <- 4*sd1^2/e1^2 #how many simulations, using precision
n1
ftimes <- systemfailuretimes(Nsim=n1)
mean(ftimes)

#Mean with precision of one day: 
e2 <- 1/31 
n2 <- 4*sd1^2/e2^2
n2
ftimes <- systemfailuretimes(Nsim=n2)
mean(ftimes) #in months 

#Illustrate various summary measures
Nsim <- 100000
ftimes <- systemfailuretimes(Nsim=Nsim)
#Various summary measures
summary(ftimes)
sd(ftimes) #Standard deviation
var(ftimes) #Variance
quantile(ftimes) #Quantiles
quantile(ftimes, probs = c(0,0.05,0.1,0.25,0.5,0.75,0.9,1)) #Same function, but specify which quantiles we are interested in

#Histogram 
hist(ftimes,nclass = 100, probability =  T)
#Empirical cdf
plot(ecdf(ftimes))
#Density estimate
plot(density(ftimes)) #gives the best estimate of how the pdf of our data looks like


#Various probabilites:
#Failure before 12 months: 
sum(ftimes<12)/Nsim #Divide by nsim to get between 0 and 1

#Failure after 15 months:
sum(ftimes>15)/Nsim

#failure between 10 and 15 months:
sum(ftimes>10 & ftimes<15)/Nsim



#Problem 1
#a)
#Sketch density of exponentially distributed function
prob <- dexp(1,x) #lambda = 1
sequence <- seq(0,1,by=0.01)
prob2 <- dexp(1,sequence) #lambda = 1
plot(prob2)

#Variance
var(prob2)
#Expectation
mean(prob2)

#b)
#Probability of a lightbulb follows exponential with expectation 10000
lifetime <- pexp(rate = 1/10000,x) #Cumulative distribution function

#Lifetime less than 10 000 hours
lifetime1 <- pexp(10000 ,rate = 1/10000,lower.tail =  TRUE)

#Lifetime more than 5000 hours
lifetime2 <- pexp(5000, rate = 1/10000, lower.tail = FALSE)

#Lifetime between 5000 and 1000 hours
lifetime3 <- pexp(10000, rate = 1/10000, lower.tail = TRUE) - pexp(5000, rate = 1/10000, lower.tail =  TRUE)



#Problem 2: Number of accidents on road per year. Poission distributed with lambda = 2

#a) 
#Exactly 2 accidents during one year P(X=2)
lambda = 2
a1 <- dpois(2,lambda = lambda) #density function

#At least 3 accidents P(X>=3)
a2 <- 1 - ppois(3,lambda = lambda) #Same as ppois(3,lambda=lambda,lower.tail=FALSE)

#b)
#Two accients in half a year p(2;0.5)
lambda = 1/2
b1 <- dpois(2,lambda = lambda)

#c) At least 20 accidents in 10 years 
lambda = 2*10 
c1 <- ppois(20,lambda=lambda, lower.tail = FALSE)


#d) More than two accidents in a year when we know there has at least been one 
lambda = 2
d1 <- (ppois(2,lambda = lambda, lower.tail = FALSE) / ppois(1,lambda = lambda, lower.tail = FALSE))

#Problem 3 
#Implement beta distributions with i) alpha=beta=1, ii)alpha=2,beta=1 and iii)alpha=1,beta=2
#Draw samples of different sizes from each of the beta distributions and calculate the average of these samples
#Compare the average to the expectation which for the beta-distribution is given by E(X)=aplha/(aplha+beta)

alpha1 = beta1 = 1
alpha2 = 2
beta2 = 1
alpha3 = 1
beta3 = 2

expectation <- function(alpha,beta) {
  return(alpha/(alpha+beta))
}

expectation(alpha = alpha1, beta = beta1)
expectation(alpha = alpha2, beta = beta2)
expectation(alpha = alpha3, beta = beta3)


#Run for Nsamples = 10,100,1000,10000,100000 and se what happens
drawSampleSizes <-function(Nsamples) { 
  B1 <- rbeta(n = Nsamples, shape1 = alpha1, shape2 = beta1)
  B2 <- rbeta(n = Nsamples, shape1 = alpha2, shape2 = beta2)
  B3 <- rbeta(n = Nsamples, shape1 = alpha3, shape2 = beta3)
  avgB1 <- mean(B1)
  avgB2 <- mean(B2)
  avgB3 <- mean(B3)
  averages <- c(avgB1,avgB2,avgB3)
  return(averages)
}


#Problem 4
#b) Let n be 5,20 and 100 and use as many replications in the simulation
#as you deem necessary to obtain a histogram that gives a good approximation
#to the distribution for X. 

Un  <- function(n) {
  return(runif(n = n, min = 0, max = 2))
}


findSumsUn <-function(Nsim,n){ #Input argument is number of simulations 
  sumVector <- vector(length = Nsim)
  for(i in 1:Nsim) {
    sumVector[i] <- sum(Un(n))
  }
  return(sumVector)
} 

plotResults <-function(Nsim,n){
  sums <- vector(length = Nsim) 
  for(i in 1:Nsim) {
    sums[i] <- findSumsUn(Nsim=Nsim,n=n)
  }
  hist(sums)
  return(sums)
}



