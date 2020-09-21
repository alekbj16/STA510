#Remove old variables
rm(list=ls())

#Problem 1
#c)

# Using the inverse transform method
genrandnumb <- function(n) {
  #Generate a random u from Uniform(0,1)
  u <- runif(n = n, min = 0, max = 1)
  #Deliver x = F-1(u)
  x <- 2-sqrt(-8*u+9)
  return(x)
}

#Historgram to verify
Nsim <- 1000000
stimeINV <- Sys.time() #start timer
X <- genrandnumb(n = Nsim) 
etimeINV <- Sys.time() #end timer
print(etimeINV - stimeINV) #print runtime
hist(X,prob = TRUE)



#e)
#Using the acceptance-rejection method
aceptreject <- function(n) {
  k <- 0 #counter for accept
  y <- vector(mode = "numeric",length = n) #double-precision vector
  while(k < n) {
    u <- runif(1,min = 0, max = 1)
    x <- runif(1,min = -1, max = 1) #random variate from g
    if (((2-x)/3)>u){ 
      #accept x
      k <- k + 1
      y[k] <- x
      }
  }
  return(y)
} 


Nsim <- 1000000

stimeAR <- Sys.time() #start timer
Xs <- aceptreject(n = Nsim)
etimeAR <- Sys.time() #end timer
print(etimeAR-stimeAR) #print runtime
hist(Xs,prob = TRUE)

#Which method would I prefer?

#By using Sys.time() before each function,
#I get that the runtime for the inverse method is 
#0.05243492 secs, while for the acceptance-
#rejection the runtime is 4.730012 secs. 
#Both cases when n = 1000000. The AR method is slower
#because it requires more iterations. I would prefer
#the inverse transform method

#f) 
#Using inverse method
mean(X)           #E(X), gives -0.1679223
var(X)            #Var(X), gives 0.3048559
mean(0.5*(X+1))   #E(0.5(X+1)), gives 0.4160388
mean(X^3)         #E(X^3), gives -0.1007686
#Using AR method
mean(Xs)          #E(X), gives -0.1665917
var(Xs)           #Var(X), gives 0.3056142
mean(0.5*(Xs+1))  #E(0.5(X+1)), gives 0.4167041
mean(Xs^3)        #E(X^3), gives -0.1000548
#For both the inverse method and the AR method
#the results match my calculations

#95% confidence interval
#using inverse method
confidence <- function(N){
  X <- aceptreject(n = N)
  upper <- mean(X) + 1.96*(sd(X)/sqrt(N))
  lower <- mean(X) - 1.96*(sd(X)/sqrt(N))
  
  interval <- c(lower,upper)
  return(interval)
}

confidence(100000) #Try for rising values, we see that the confidence interval gets reduced ("higher certainty")



#Problem 2
#b)

#Define distribution A as a function
distA <- function(n,mu,sigma) {
  Ai <- sum(rnorm(n = n, mean = mu, sd = sigma))
  return(Ai)
}


#Constants
n = 100
mu = 1
sigma = 2
Nsim = 100000
As <- vector(mode = "numeric", length = Nsim)
Bs <- vector(mode = "numeric", length = Nsim)
#Simulate the Distribution of A Nsim times. Each sim creates 100 values
for(i in 1:Nsim)
  As[i] <- distA(n,mu,sigma)
  
print("Mean and variance for A: ")
mean(As) #Yields 100.0089, which is close to the analytical n*mu = 100*1 = 100.
var(As) #Yields 400.5617, which is close to the analytical n*sigma^2 = 100*2^2 = 400.
sd(As) #sqrt of var

#Distribution B
Bs <- As/n
print("Mean and variance for B: ")
mean(Bs) #Yields 1.000089, which is close to the analytical mu = 1.
var(Bs) #Yields 0.04005617, which is close to the analytical sigma^2/n = 0.04.
sd(Bs) #sqrt of var

#Distribution C
print("Mean and variance for C: ")
Cs <- (sqrt(n)*(Bs-mu))/sigma
mean(Cs) #Yields 0.0004468529, which is close to the analytical mu = 0.
var(Cs) #Yields 1.001404, which is close to the analytical mu = 1.
sd(Cs) #sqrt of var
 
#d)
#Simulate many outcomes of of X1, X2, ...,Xn for the n obtained. 

n <- 6400
mu <- 1
sigma <- 2
Nsim <-1000 #exercise said "many" outcomes


CIfunc <- function(n, mu, sigma,Nsim){
  c <- 0 #Create a counter for how many in the interval
  for(i in 1:Nsim){
    Xs <- rnorm(n,mu,sigma)
    B <- mean(Xs) 
    CI <- -1.96*sd(Xs)/sqrt(n)
    if(B-CI>mu & mu>B+CI) { #if within interval
      print(c(B+CI,B-CI)) 
      c <- c + 1 #counter for percent within interval
    }
  }
  print(c("Percent of times that interval contains true mean", c / Nsim * 100))
} 

hello <-CIfunc(n,mu,sigma,Nsim)

#Redo with n = 10, what do you observe? 
#When redoing with n = 10, I observe that the confidence intervals become larger
#and that the percent also is reduced a little bit

#Problem 3

#Resolution
res = 10^6

#Testing loading the data
source("seir_model.R")
f <- seir_sim()
names(f)
head(f$ode)
head(f$dates)

#Suceptible - 15.03.2021
t <- which(f$dates=="2021-03-15") # the row index of 15th march of 2021
numb <- f$ode[t,"S"]
print(c("Persons in state S (suceptible) on 15th of March 2021: ",numb*res))

#Max CC in entire period
days <- length(f$dates) #Number of days in f
ccmax <- 0
daynumber <- 0
for(i in 1:days) {
 ccdaily <- f$ode[i,"CC"]
 if(ccdaily > ccmax){
   ccmax <- ccdaily
   daynumber <- i
 }
}
ccmaxPerMillion <- ccmax * res
date <- f$dates[daynumber]
print("Date at which CC is highest (per million): ")
date
print("At that date, CC is (per million): ")
ccmaxPerMillion


#Avg number of persons in CC in 2021
startyear <- which(f$dates=="2021-01-01")
endyear <- which(f$dates=="2021-12-31") 
days2021 <-endyear - startyear +1 #+1 to get entire year (365 days)

inCC <- 0 
for(d in 1:days2021){
  inCC <- inCC + f$ode[d+366,"CC"] #+366 to start in 2021
}

avgNumber <- (inCC*res)/days2021

print("The average number of people in CC in 2021 is: ")
print(avgNumber)


#b)
f2 <- seir_sim(upper_thresh = 10000)

#Max CC in entire period when upper_thresh=10000
days <- length(f2$dates) #Number of days in f2
ccmax2 <- 0 #variable to store the max
daynumber <- 0 #daily value
for(i in 1:days) {
  ccdaily <- f2$ode[i,"CC"] #get the daily value
  if(ccdaily > ccmax2){  
    ccmax2 <- ccdaily #store max in the ccmax2 variable
    daynumber <- i #store the daynumber at which the max value occurs
  }
}
ccmaxPerMillion2 <- ccmax2 * res
date <- f2$dates[daynumber]
print("Date at which CC is highest (per million, when upper_tresh=10000): ")
date
print("At that date, CC is (per million): ")
ccmaxPerMillion2
print("We see that CC is notecibly higher than when social distancing is turned off")

#An epidemic ends when fewer than 1 per 1  million are in the E state
#Virus is introduced on 13th of march 2020
#Model runs from 1st January 2020 until 1st January 2025

days <- length(f2$dates)
virus <- which(f2$dates=="2020-03-13")
endepidemic <- 0
for(d in virus+1:days){ #+1 to get the date after which the virus is introduced
  if(1 > (f2$ode[d,"E"])*res){
    endepidemic <- d
    break
  }
}

freedom <- f2$dates[endepidemic] # the row index of new years eve 2021
print("The epidemic ends (with upper_thres = 10000):")
print(freedom) 


# Portion in RR, RH or RC when the epidemic ends
infected <- f2$ode[endepidemic,"RR"]*res + f2$ode[endepidemic,"RH"]*res + f2$ode[endepidemic,"RH"]*res
print("Portion of the population that has been infected when the pandmic ends (per million): ")
infected


#c)
#Uncertain about reproduction 
n <- 1000 #Number of simluations, at least 100
#mean and standard deviations of the maximum numbers
#of critical care beds required in 2020. 

startDate <-which(f$dates=="2020-01-01")
endDate <-which(f$dates=="2020-12-31")


days <- endDate - startDate + 1 #+1 to account for ending day
results <- matrix(nrow = n, ncol = days)
listOfR0Max <- vector(mode = "numeric",length = n)
listOfSeason_amp <- vector(mode = "numeric",length = n)
for(x in 1:n) {#Double for loop, should be avoided but is hard to avoid when working with matrix
  R0max <- runif(1, min = 1.8, max = 3.0)
  listOfR0Max[x] <- R0max
  season_amplitude <- runif(1,min = 0.0, max = 0.4)
  listOfSeason_amp[x] <- season_amplitude
  ff <- seir_sim(R0max = R0max, season_amplitude = season_amplitude)
  for(j in 1:days) { #for each day in 2020, find the CC
    results[x,j] <- res*ff$ode[j,"CC"]
    }
  }
results

#CCmax <- function(inputMatrix,Nsim){ #Should be same Nsim as in CCdays
maximumPerSim <- vector(mode="numeric",length=n)

for(y in 1:n){
  maximumPerSim[y] <- max(results[y,])
  }
maximumPerSim

#Mean max CC
mean(maximumPerSim)
#Standard deviation
sd(maximumPerSim)
#Quantiles
quantiles <- quantile(maximumPerSim, probs = c(0.90,0.95,0.99)) #Source number 1 (see end of pdf)
quantiles

#P(more than 100 CC beds)
iterations <- length(maximumPerSim)
k <- 0
for(i in 1:iterations) {
  if(maximumPerSim[i]>100) {
    k <- k + 1
  }
}
probability <- k/iterations
print("More than 100 beds per million: ")
probability
print("Or in (%): ")
probability*100

#Plot for density (relative frequency)
hist(maximumPerSim,main = "Maximum required\n CC beds in 2021 ", probability  = TRUE)
#"Regular" plot
hist(maximumPerSim,main = "Maximum required\n CC beds in 2021 ",sub = "n = 1000 simulations")
#A normal distribution should be symmetrically distributed around the mean.
#Here, we do not see this behavior. 
#Therefore, this is not well approximated by the mean. 


plot(listOfR0Max,maximumPerSim,main="Max CC 2021\n vs\n R0max~U(1.8,3.0)",xlab="R0Max ~Unif(1.8,3.0)",ylab="Maximum required CC in 2021",sub="n = 1000 simulations")

plot(listOfSeason_amp,maximumPerSim,main="Max CC 2021\n vs\n seanson_amplitude~U(0.0,0.4)",xlab="season_amplitude~Unif(0.0,0.4)",ylab="Maximum required CC in 2021",sub="n = 1000 simulations")

#d)
#Between max CC beds and R0max, a linear relationsship
#can be interpreted in the range [2.2,2.8]. However,
#when examining the plot as a whole, there seems to be 
#a resemblence to a exponential distribution, or at least
#the start of it. What I mean by this is that is seems like 
#maxCC starts to increase quite a lot when R0max increases,
#and one may expect maxCC to continue increasing if the R0Max 
#range would have been larger, by the look of the current plot. 

#In the plot between max CC beds vs seasonal amplitude, we
#observe that MaxCC is reduced as seasonal_amplitude increases
#This is due to the effect that the seasonal amplitude has on
#R0max, which is R0max*(1-seasonal_amplitude) in July. 
#One might expect these to be linear. However, the seasonal_amplitude
#has a different effect in january, where R0max is unaffected by seasonal
#amplitude. The seasonal_amplitude effect increases towards the summer
#and lowers towards fall and the end of the year. Thus, we do not 
#see any strong linear relationship. 


#In total, both of these numbers affect the number of critical care beds,
#as seasonal_amplitude affects R0max. However, R0max has more effect on
#the number of CC beds, while seasonal_amplitude works as a "damper". 



#Using pearson method to check for linear relationship
#Source number 2 and 3 (see end of pdf)
cor(listOfR0Max,maximumPerSim,method="pearson")
#A value of 0 implies that there is no linear correlation between 
#the variables, while 1 implies that a linear equation describes the 
#relationship between listOfR0Max and maximumPerSim. A value of 0.77
#therefore implies a somewhat noticeable linear correlation, but the relationship
#cannot be expressed by a linear equation. 



#Using pearson method to check for linear relationship
cor(listOfSeason_amp,maximumPerSim,method="pearson") 
#A value of -0.343 implies that the relationsship is negative.
#It cannot be described by a linear equation. 
#The correlation was more noticeable in the correlation
#study of listOfR0Max and maximumPerSim, meaning that
#R0max has greater effect on CC beds than season_amplitude.


