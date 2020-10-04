#Weekly examples from lecture as well as exercises

#Remove old variables
rm(list=ls())

#Define some constants
gam <- 0.1
phi <- 0.9
sigma <- 0.1

#The stationary mean and standard dev as calculated analytically
stationaryMean <- gam/(1-phi)
stationarySD <- sqrt((sigma^2)/(1-phi^2))
print(paste0("Stationary mean: ",stationaryMean,". Stationary SD: ",stationarySD))

T <- 1000  #Length of simulations
N <- 10000 #Number of simulations 

a <- -1.0 #Initial configuration, same in ann relaizations

simulations <- matrix(0.0,nrow = N, ncol = T)

#Each simulations starts at the value a:
simulations[,1] <- a

#The itereate and calculate next step
for( n in 2:T){
  simulations[,n] <- gam + phi*simulations[,n-1] + rnorm(N,mean=0.0, sd = sigma)
}

#Visualizations: 
if(F){ #change to T
  par(mfrow=c(2,1))
  plot(1:T,simulations[1,],type="l",xlab = "time n", ylim = c(min(simulations),max(simulations))) # first trajectory
  for(i in 2:min(N,50)){ # take only first 50 realizations to save computing time
    lines(1:T,simulations[i,],col=i)
  }
  
  hist(simulations[,T],probability = TRUE,breaks=21)
  if(abs(phi)<1.0){ # only stationary case
    xgrid <- seq(from=stationaryMean-4*stationarySD,
                 to=stationaryMean+4*stationarySD,
                 length.out = 1000)
    lines(xgrid,dnorm(xgrid,mean=stationaryMean,stationarySD))
  }
}
if(T){ #change to F
  meanEstAR <- rowMeans(simulations)
  meanEstIId <- rowMeans(matrix(rnorm(N*T,mean=stationaryMean,sd=stationarySD),N,T))
  msd <- max(sd(meanEstAR),sd(meanEstIId))
  hist(meanEstAR,xlim=c(mean(meanEstAR)-3*msd,mean(meanEstAR)+3*msd),breaks=21)
  hist(meanEstIId,xlim=c(mean(meanEstAR)-3*msd,mean(meanEstAR)+3*msd),breaks=21)
  
  print(paste0("Standard deviation based on AR(1) samples: ",sd(meanEstAR)))
  print(paste0("Standard deviation based on iid samples: ",sd(meanEstIId)))
}



##################
#Monte Carlo integration

#clear variables
rm(list=ls())

#Integrate cos(x) from 0 to 1

x <- runif(10000) #Draw a large number of uniform random numbers
mean(cos(x)) #Calculate the average of the cos function applied to the random uniform numbers

#Exact:
sin(1) #Since the integral is so simple, we know the exact value



###Integrate cos(x) from a to b 

cosintegral <- function(Nsim = 10000, a, b) {
  #Simulate uniform variables between a and b
  x <- runif(Nsim,a,b)
  integralresult <- (b-a)*mean(cos(x))
  print(paste("Integrating cos(x) from ",a,"to ", b))
  print(integralresult)
}
  
cosintegral(a=0,b=1)
cosintegral(a = 0, b = 2)
cosintegral(a = pi, b = 2*pi)



###Find the cumuluative distribution of gamma(4,0.5)-dist
#Version 1 implementing unif numbers
approxgamma <- function(xvec,Nsim=10000){
  cdf <- numeric(length(xvec))
  for(i in 1:length(xvec)) {
    u <- runif(Nsim, 0, xvec[i])
    cdf[i] <- xvec[i]*mean((8/3)*u^3*exp(-2*u))
  }
  cdf
  
}

#Version 2, using indicator function I(u<x)
approxgamma2 <- function(xvec,Nsim =10000){
  cdf <- numeric(length(xvec))
  for(i in 1:length(xvec)) {
    x <- rgamma(Nsim, shape = 4, scale = 0.5)
    cdf[i] <- mean(x<xvec[i])
  }
  cdf
}

xvec <- seq(0.25,5,by=0.25)
cdf1 <- approxgamma(Nsim = 10000, xvec = xvec)
cdf2 <- approxgamma2(Nsim = 10000, xvec = xvec)
exact <- pgamma(xvec,shape=4, scale = 0.5)

plot(xvec,cdf1,type="l",col="blue")
lines(xvec,cdf2,col="red")
lines(xvec,exact)
legend("bottomright",c("Version 1","Version 2", "Exact"), col = c("blue","red","black"),lty=1)



#### integral of x^5e^-3x from 0 to infinity
x <- rexp(1000000, rate = 3)
mean((x^5)/3)


### integrate function of three RVs

n <- 100000
x <- runif(n,-1,1)
y <- runif(n,2,3)
z <- runif(n,0,1)
int <- 2*1*1*mean(x*y+z*y)
int

###Weekly exercises in rizzo
rm(list=ls())
#4.3)


nonhomogenouspois <- function(low = 4, up = 5, maxval = 12) {
  # First simulate HPP with intensity maxval on low to up
  expectednumber <- (up-low)*maxval 
  Nsim <- 3*expectednumber  # Simulate more than the expected number to be certain to exceed stoptime
  Tn <- rexp(Nsim, maxval) # Simulate interarrival times
  Sn <- low + cumsum(Tn) #only interested in times over low
  Sn <- Sn[Sn<up] #Drop values larges than max
  Nevents <- length(Sn) #counter
  Un <- runif(Nsim)
  keep <- (Un <= (2*Sn + 2)/maxval)
  Sn[keep]
  return(round(Sn[keep], 4))
}

nonhmgpois <- nonhomogenouspois()

#######
#6.1)
x <- runif(10000,min = 0, max = pi/3)
mean(sin(x))*(pi/3)

#Real solution
cos(pi/3)

#6.2)
approxpnorm <- function(xvec,Nsim=10000){
  cdf <- numeric(length(xvec)) 
  for(i in 1:length(xvec)){
    u <- runif(Nsim,0,xvec[i])
    cdf[i] <- 0.5+xvec[i]*mean(exp(-u^2/2))/sqrt(2*pi)
  }
  cdf  
}

# Run and compare to pnorm
xvec <- seq(0.25,2.5,by=0.25)
cdfest <- approxpnorm(xvec)
Phi <- pnorm(xvec)
round(rbind(xvec,cdfest,Phi),4)

#6.3)
#Using unif
x <- runif(1000,min = 0, max = 0.5)
mean1 <- mean(exp(-x))*0.5
var1 <- var(mean1)
mean1
var1

#Using exponential dist
x <- rexp(1000,rate=1)
(theta.hat <- mean(x<=0.5) ) # I.e. proportion of x's <= 0,5
# Standard error:
sqrt(theta.hat*(1-theta.hat)/1000 )


#####
#Problems on excercise sheet
#1

#Simulate two nonhomogenous poisson processes

#Start by defining intesity functions

intensity1 <- function(x) {
  lambda <- numeric(length(x))
  lambda[ceiling(x) %% 2 == 0] <- 2 #Using ceiling as hinted of
  lambda[ceiling(x) %% 2 == 1] <- 1
  lambda
}


intensity2 <- function(x) {
  lambda <- numeric(length(x))
  lambda[ceiling(x) %% 2 == 0] <- 1.5 #Using ceiling as hinted of
  lambda[ceiling(x) %% 2 == 1] <- 1
  lambda
}

NHHsimulator <- function(low,up,maxval,intensity) {
  #First simulate HPP with intensity lambdamax on a to b
  expect <- (up-low)*maxval
  Nsim <- 3*expect #generate enough simulations
  Tn <- rexp(Nsim,maxval) # Simulate interarrival times
  timesto <- low+cumsum(Tn)   # Calculate arrival times starting at low
  timesto <- timesto[timesto<up] # drop the times larger than b
  counter <- length(timesto) #counter
  U <- runif(counter) #thinning
  timesto <- timesto[U<intensity(timesto)/maxval]  
  timesto  # Return the remaining times
}



NHHPtimes1 <- NHHsimulator(low=0,up=20,maxval = 2,intensity = intensity1)
NHPPtimes1
plot(NHHPtimes1)

NHHPtimes2 <- NHHsimulator(low=0,up=20,maxval = 1.5,intensity = intensity2)
NHHPtimes2
plot(NHHPtimes2)

######
#Problem 2
#Gaussian AR process
ARsim <- function(T=100,mu=0.0,sigma=1.0,phi=0.0) {
  x <- numeric(length(t))
  x[1] <-  rnorm(n=1,mean=mu,sd=sigma)
  innov.sd <- sigma*sqrt(1.0-phi^2)
  for(t in 2:T) x[t] <- mu + phi*(x[t-1]-mu) + rnorm(n=1,mean=0.0,sd=innov.sd)
  return(x)
}

phi <- -0.9 # change to get remaining cases
Nsim <- 1000
resM1 <- vector(mode="numeric",length=Nsim)
resM2 <- vector(mode="numeric",length=Nsim)

for(i in 1:Nsim){
  x <- ARsim(phi=phi)
  resM1[i] <- mean(x)
  resM2[i] <- mean(x^2)
}
print("part b")
print(paste0("phi = ",phi),quote=FALSE)
print(paste0("M_1: mean: ",mean(resM1),"  variance: ",var(resM1)),quote=FALSE)
print(paste0("M_2: mean: ",mean(resM2),"  variance: ",var(resM2)),quote=FALSE)


