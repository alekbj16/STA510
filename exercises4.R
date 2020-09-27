###Examples from lectures week 21.09.20-25.09.20
# as well as solution to exrcises 
# Disclaimer: This weeks solutions are highly inspired  from the solutions given by the course instructor


##Remove old variables
rm(list=ls())

#Homogenous poisson process (HPP)

#First implementation of a poisson process
#Slow but clear
plotHPP <- function(lambda, stoptime){
  timesto <- 0.0
  t <- 0.0
  while(TRUE){ # run while the process has not finished
    timebetween <- rexp(1,lambda) # time between events
    t <- t + timebetween # process time at last event
    if(t < stoptime) {
      timesto <- c(timesto,t) # resizing vectors should be avoided !!
    } else {
      break # done
    }
  } #end while
  Nevents <- length(timesto)
  plot(timesto,1:Nevents,type="s",xlab = "arrival time", 
       ylab = "Event number",lwd=1.5,ylim=c(0,Nevents))
  points(timesto,rep(0,Nevents),pch=21,bg="red")
}

plotHPP(lambda = 0.05,stoptime=365)


# vectorized implementation based on simulating times between events
plotHPPv <- function(lambda,stoptime){
  expectednumber <- stoptime*lambda  # Expected number of event until stoptime
  Nsim <- 3*expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
  timesbetween <- rexp(Nsim,lambda) # Simulate interarrival times
  timesto <- cumsum(timesbetween)   # Calculate arrival times
  timesto <- timesto[timesto<stoptime] # Dischard the times larger than stoptime
  Nevents <- length(timesto) # Count the number of events
  plot(timesto,1:Nevents,type="s",xlab = "arrival time", 
       ylab = "Event number",lwd=1.5,ylim=c(0,Nevents))
  points(timesto,rep(0,Nevents),pch=21,bg="red")
}

plotHPPv(lambda=0.05,stoptime=365)

# implementation based on the relation to the uniform distribution
plotHPPu <- function(lambda,stoptime){
  Nevents <- rpois(1,lambda*stoptime) # total number of events
  timesto <- sort(runif(Nevents,min=0,max=stoptime)) # the events are uniformly distributed conditional on the number of events
  plot(timesto,1:Nevents,type="s",xlab = "arrival time", 
       ylab = "Event number",lwd=1.5,ylim=c(0,Nevents))
  points(timesto,rep(0,Nevents),pch=21,bg="red")
}
plotHPPu(lambda = 0.05, stoptime = 365)



### Renewal process (RP)

# Use a gamma distribution for times between events 
# (often called gamma renewal process)
plotgammaRP <- function(alpha,beta,stoptime){
  expectedtimebetween <- alpha*beta # Expected interarrival time
  print(paste("Expected time between events: ",expectedtimebetween, " with variance:",expectedtimebetween*beta))
  expectednumber <- stoptime/expectedtimebetween # Expected number of event until stoptime
  Nsim <- 3*expectednumber # Simulate more than the expected number to be certain to exceed stoptime
  timesbetween <- rgamma(Nsim,shape = alpha, scale = beta) # Simulate interarrival times
  timesto <- cumsum(timesbetween)     # Calculate arrival times
  timesto <- timesto[timesto<stoptime]  # Dischard the times larger than stoptime
  Nevents <- length(timesto)  # Count the number of events
  plot(timesto,1:Nevents,type="s",xlab = "arrival time", 
       ylab = "Event number",lwd=1.5,ylim=c(0,Nevents),
       main=paste("expectation=",expectedtimebetween, ", variance=",expectedtimebetween*beta))
  points(timesto,rep(0,Nevents),pch=21,bg="red")
}
plotgammaRP(alpha=0.2,beta=30,stoptime=365)
plotgammaRP(alpha=2,beta=3,stoptime=365)
plotgammaRP(alpha=20,beta=0.3,stoptime=365)

# Generat data from four different RPs with the same
# expected interarrival time but different variances
par(mfrow=c(2,2))
plotgammaRP(alpha=0.02,beta=100,stoptime=100)
plotgammaRP(alpha=0.2,beta=10,stoptime=100)
plotgammaRP(alpha=2,beta=1,stoptime=100) # HPP
plotgammaRP(alpha=20,beta=0.1,stoptime=100)
par(mfrow=c(1,1))
# Repeat the lines above several times


# NHPP describing arrivals of cars 


# Function for simulating arrival times for a NHPP between a and b using thinning
simtNHPP <- function(a,b,lambdamax,lambdafunc){
  # Simple check that a not too small lambdamax is set
  if(max(lambdafunc(seq(a,b,length.out = 100)))>lambdamax)
    stop("lambdamax is smaller than max of the lambdafunction")
  # First simulate HPP with intensity lambdamax on a to b
  expectednumber <- (b-a)*lambdamax  
  Nsim <- 3*expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
  timesbetween <- rexp(Nsim,lambdamax) # Simulate interarrival times
  timesto <- a+cumsum(timesbetween)   # Calculate arrival times starting at a
  timesto <- timesto[timesto<b] # Dischard the times larger than b
  Nevents <- length(timesto) # Count the number of events
  # Next do the thinning. Only keep the times where u<lambda(s)/lambdamax
  U <- runif(Nevents)
  timesto <- timesto[U<lambdafunc(timesto)/lambdamax]  
  timesto  # Return the remaining times
}

# Specify the intensity function for the traffic example
lambdatraffic <- function(t)
  200+190*sin(2*t-1)
# Plot the intensity function
tvec <- seq(0,2,by=0.01)
plot(tvec,lambdatraffic(tvec),type="l",ylim=c(0,400))

# Generate data with the traffic intensity and plot them
NHPPtimes <- simtNHPP(a=0,b=2,lambdamax=390,lambdafunc=lambdatraffic)
plot(NHPPtimes,1:length(NHPPtimes),type="s",xlab = "time", 
     ylab = "Event number",lwd=1.5)
points(NHPPtimes,rep(0,length(NHPPtimes)),pch=21,bg="red")
# Rerun the lines above several times



# Specify a different intensity function
lambdacyclic <- function(t)
  200*cos(t*2)^2
# Plot the intensity function
tvec <- seq(0,2,by=0.01)
plot(tvec,lambdacyclic(tvec),type="l",ylim=c(0,200))

# Generate and plot them
NHPPtimes <- simtNHPP(a=0,b=2,lambdamax=200,lambdafunc=lambdacyclic)
plot(NHPPtimes,1:length(NHPPtimes),type="s",xlab = "arrival time", 
     ylab = "Event number",lwd=1.5)
points(NHPPtimes,rep(0,length(NHPPtimes)),pch=21,bg="red")
# Rerun the lines above several times



# The web-page example
lambdaexp <- function(t)
  800*exp(-0.1*t)
# Plot the intensity function
tvec <- seq(0,10,by=0.01)
plot(tvec,lambdaexp(tvec),type="l",ylim=c(0,800))

# One way to verify the calculations for the number of vistors during 10 hours
# is to simulate data from the NHPP on [0,10] many times and count the number
# of arrivals each time:
Nsim <- 200
NHPPnumbers <- vector(length=Nsim)
for(i in 1:Nsim)
  NHPPnumbers[i] <- length(simtNHPP(a=0,b=10,lambdamax=800,lambdafunc=lambdaexp))
# Average
mean(NHPPnumbers)
# Estimated probability
mean(NHPPnumbers>5000)


#################################
##Remove old variables
rm(list=ls())


#Problem 3.1 in Rizzo:
#Quantile
quants <- function(q,lambda,eta){
  -(log(1-q))/lambda + eta
}

#Random number
rquant <- function(n, lambda, eta) {
  quants(runif(n),lambda,eta)
}

n <- 100
lambda <- 2.0
eta <- 0.5
x <- rquant(n = n, lambda = lambda, eta = eta)

# The quantiles to compute
alphas <- seq(from=0.01,to=0.99,by=0.01)

#Theoretical quantiles
theoretical <- quants(q = alphas,lambda = lambda, eta = eta)

#The corresponding empirical quantiles using built-in function quantile
empquants <- quantile(x=x,probs=alphas)

# Compare the numbers
round(rbind(theoretical,empquants),3)

# Plot the empirical aganist theoretical quantiles, should be approximately
# on the x=y-line (try googling qq-plot for explanation)
plot(x=empquants,y=theoretical)
# add x=y-line
abline(0,1 ,col="red")

# Plot of probabilities versus quantiles:
plot(empquants,alphas,type="l",xlab="quantile",ylab="probability",lwd=2.5)
lines(theoretical,alphas,col="red",lwd=2.5)


#################################
##Remove old variables
rm(list=ls())


#Problem 3.11 in Rizzo:
#Mixture generation
Nsim <- 1000
p <- 0.25
meanvec <- c(0,3)
sdvec <- c(1,1)
probsvec <- c(p,(1-p))


#Create function
genmixofnorm <- function(Nsim, meanvec, sdvec, probsvec) {
  m <- length(probsvec) # Number of distributions
  # Determine which distribution to sample from in each simulation
  whichdist <- sample(1:m,size=Nsim,replace = TRUE,prob = probsvec)
  # Simulate the mixture
  Y <- rnorm(Nsim,mean = meanvec[whichdist],sd = sdvec[whichdist])
  Y
}

Y <- genmixofnorm(Nsim = Nsim, meanvec = meanvec, sdvec = sdvec, probsvec = probsvec)

hist(Y,prob = TRUE,nclass=max(10,sqrt(Nsim)),main="Histogram of normal mixture")
curve(probsvec[1]*dnorm(x, mean=0, sd=1)+probsvec[2]*dnorm(x, mean=3, sd=1),from=-4,to=7,add=TRUE, col='red',lwd=1.6)


#################################
##Remove old variables
rm(list=ls())


#Problem 4.1 in Rizzo:

score <- 10+sample(c(-1,1),size=1) # Step 1
Sn <- score  # Step 1
# Next steps - stop when score reaches 0 or 20
while(score>0 & score<20){
  score <- score+sample(c(-1,1),size=1)
  Sn <- c(Sn,score)
}
tstop <- length(Sn)
plot(1:tstop,Sn[1:tstop],type="b",ylim=c(0,20),
     xlab="n",ylab="Sn")


#################################
##Remove old variables
rm(list=ls())

#Problem from sheet
### Problem 1

# a)

# Function to simulate Homogenous possion proecss
simHPP <- function(lambda,stoptime){
  narrivals <- rpois(1,lambda*stoptime)
  timesto <- sort(runif(narrivals,min=0,max=stoptime))
  return(timesto)
}

# Test
timesto <- simHPP(lambda=0.2,stoptime=365)
Nevents <- length(timesto)
plot(timesto,1:Nevents,type="s",xlab = "arrival time", 
     ylab = "Event number",lwd=1.5,ylim=c(0,Nevents))
points(timesto,rep(0,Nevents),pch=21,bg="red")

# b) 
Nsim <- 10000 
N10vector <- vector(length=Nsim)
for(i in 1:Nsim)
  N10vector[i] <- length(simHPP(lambda=2,stoptime=10))
mean(N10vector>25)

# Or more directly whithout generating all the times to events:
N10vector <- rpois(n=Nsim,20)
mean(N10vector>25)

# Exact
1-ppois(25,2*10)
