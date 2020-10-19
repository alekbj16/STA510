#Mandatory 2
#Student number 236227
#Aleksander B. Jakobsen

#Remove variables
rm(list=ls())

#1c)
nrep <- 1000
estimates <- vector(mode = "numeric", length = nrep) 
for(est in 1:nrep){
  n = 1000
  x <- runif(n, min = -1, max = 1)
  y <- runif(n, min = -1, max = 1)
  Ids <- vector(mode="integer",length = n)
  #value <- replicate(n,0.0)  
  for(i in 1:n) {
    value <- x[i]^2 + y[i]^2
    if(value>1){
      Ids[i] <- 0
    }
    else{
      Ids[i] <- 1
    }
  }
  estimates[est] <-(4)*mean(Ids)
}


#Check that the results are approximately normally distributed:

print(paste0("Simulated mean: ",mean(estimates)))
print(paste0("Analytical mean: ",  pi))
print(paste0("Simulated variance: ",var(estimates)))
print(paste0("Analytical variance: ",(pi/n)*(4-pi)))  
hist(estimates,probability = TRUE)
curve(dnorm(x, mean=mean(estimates), sd=sd(estimates)), 
      col="red", lwd=2, add=TRUE, yaxt="n")
#Comment: When running this with nrep = 10000, the resulting histogram shows a approximatley 
#normal distribution. The simulated and analytical mean are equal unit the fourth decimal point. 
#The simulated and analytical variance are equal until the fifth decimal point. 
#However, this is a very large number of repetitions and is computationally slow since a double for loop is used.
#When reducing nrep to for example 100, one might not exactly recognize the distribution as normal.
#Further, when reducing nrep, the simulated values are more off compared to the analytical values. 

#Finally, by adding a dnorm curve to the plot, we actually see that the
#resulting historgram is well approximated by a normal distribution.


#1d)
counter <- 0
for(est in 1:length(estimates)){
  current <- trunc(estimates[est] * 100)
  if(current == 314){
    counter <- counter + 1 
  }
}
print(paste0("Probability that estimate is correct to the second decimal: ", counter/length(estimates)))


#1f) #Verify answer from 1e)
nrep = 1000 
antithetic_verify <- vector(mode = "numeric", length = nrep)
for(ant in 1:nrep) {
  a = -1 
  b = 1
  x <- runif(n = n/2, min = a, max = b)
  v <- a + b - x 
  x_v <- c(x,v) #antithetic variables consist of original and shifted
  y <- runif(n = n/2 ,min = a, max = b)
  w <- -y
  y_w <- c(y,w) #same here
  Ids_ant_verify <- vector(mode="integer",length = n) #calculate the function using antithetic
  for(i in 1:n) {
    val <- x_v[i]^2 + y_w[i]^2
    if(val>1){ #check conditions
      Ids_ant_verify[i] <- 0 
    }
    else{
      Ids_ant_verify[i] <- 1
    }
  }
  antithetic_verify[ant] <- (4)*mean(Ids_ant_verify) 
}

mean(antithetic_verify)
var(antithetic_verify) #Note that variance is higher than in 1c), and thus my conclusion from point e) is correct.


#1g) 
shift <- function(u) {
  return(((u+2.0)%%2.0)-1.0)
}

#Antithetic
nrep <- 1000
antithetic <- vector(mode = "numeric", length = nrep)
for(ant in 1:nrep) {
  n = 1000
  x <- runif(n/2, min = -1, max = 1)
  y <- runif(n/2, min = -1, max = 1)
  v <- shift(x)
  w <- shift(y)
  x_v <- c(x,v)
  y_w <- c(y,w)
  Ids_ant <- vector(mode="integer",length = n)
  for(i in 1:n) {
    val <- x_v[i]^2 + y_w[i]^2
    if(val>1){
      Ids_ant[i] <- 0
    }
    else{
      Ids_ant[i] <- 1
    }
  }
  antithetic[ant] <- (4)*mean(Ids_ant)
}

mean(antithetic)
var(antithetic) #better 

#We see that the variance is reduced compared to crude Monte Carlo


#1h)
#Importance MC
n <- 1000
sigma_1 <- 0.3
sigma_2 <- 0.624127
sigma_3 <-1.0
sigmas <- c(sigma_1,sigma_2,sigma_3)

fxy <- function(x,y,sd) {
  (exp(-(x^2)/(2*sd^2)))*(exp(-(y^2)/(2*sd^2)))/(2*pi*sd^2)
}


imp_samp <- function() { #Importance sampling
  tests <- vector(mode = "numeric", length = length(sigmas))
  for(s in 1:length(sigmas)) {
    x <- rnorm(n,mean = 0, sd = sigmas[s])
    y <- rnorm(n,mean = 0, sd = sigmas[s])
    proposal <- fxy(x,y,sigmas[s])
    val <- (x^2 + y^2 <= 1)*(x>= -1 & x<=1)*(y>=-1 & y<=1)
    test <- val/proposal
    tests[s] <- mean(test)
  }
  return(tests)
}
tests <- imp_samp()
tests


#Now lets repeat this several times:
Nsim <- 1000
testsNsim <- matrix(nrow = Nsim, ncol = 3)
for(i in 1:Nsim){
  testsNsim[i,] <- imp_samp()
}

mean(testsNsim[,1])
var(testsNsim[,1])
mean(testsNsim[,2])
var(testsNsim[,2]) #We see that this sd has the smalles variance, and therefore gives the best possible imporance sampling
mean(testsNsim[,3])
var(testsNsim[,3])


#In the crude Monte Carlo estimator for 1000 repetitions, I got 
#a simulated variance of 0.00255819179179179. In simulating 1000 times
#with the importance sampling routine, I get a variance of 0.005390187. 
#Thus, the variance of the crude Monte Carlo is smaller, and therefore performs better
#in these terms. 

#########################


#Problem 2:
#2a) 
# X has a Poisson distribution, ref lecture note on chap4, page 10.

intensity <- function(t) {
  return((297*(1+cos(2*pi*(t+(1/10))))*(1-(1/2)*exp(-t/10))/10 + 3/5))
}

#Expected number in 2025
expected_storms25 <- integrate(intensity,5,6)
print(paste0("Expected number of storms in 2025: ",expected_storms25$value))

#Expected number in 2020 and 2021 combined
expected_storms20_21 <-integrate(intensity,0,2)
print(paste0("Expected total number of storms in 2020 and 2021 combined: ",expected_storms20_21$value))
sd <- sqrt(expected_storms20_21$value)
print(paste0("Standard deviation: ", sd))


#2c) 
#NHPP describing the times storms occur in some time interval [a,b]
intensitymax <- 60 #Found analytically in 2b)

# Function for simulating arrival times for a NHPP between a and b using thinning
simtNHPP <- function(a,b,intensitymax,intensityfunc){ 
  if(a<0 )
    stop("a is smaller than 0")
  if(b<a)
    stop("b is not larger than a")
  # First simulate HPP with intensitymax on a to b
  expn <- (b-a)*intensitymax  
  Nsim <- 3*expn  # Simulate more than the expected number to be certain to exceed stoptime
  timesbetween <- rexp(Nsim,intensitymax) # Simulate interarrival times
  time_until <- a+cumsum(timesbetween)   # Calculate arrival times starting at a
  time_until <- time_until[time_until<b] # Dischard the times larger than b
  Nevents <- length(time_until) # Count the number of events
  # Next do the thinning. Only keep the times where u<intensity/intensitymax
  U <- runif(Nevents)
  time_until <- time_until[U<intensityfunc(time_until)/intensitymax]  
  time_until  # Return the remaining times
}

#One way to verify the calculations in point a), and thus verifying the algorithm,
#is to simulate data from the NHPP on [5,6] many times and count the number of storms 
#each time. 

#Storms in 2025:
Nsim <- 1000
NHPPnumbers25 <- vector(length=Nsim)
for(i in 1:Nsim)
  NHPPnumbers25[i] <- length(simtNHPP(a=5,b=6,intensitymax = 60,intensityfunc = intensity))
# Average
avg25 <- mean(NHPPnumbers25)
print(paste0("The simulated number of storms in 2025 is ",avg25,". This corresponds nicely to the value obtained in 2a), which was ",expected_storms25$value))


#Storms in 2020 and 2021
NHPPnumbers20_21 <- vector(length=Nsim)
for(i in 1:Nsim)
  NHPPnumbers20_21[i] <- length(simtNHPP(a=0,b=2,intensitymax = 60,intensityfunc = intensity))
# Average
avg20_21 <- mean(NHPPnumbers20_21)
print(paste0("The simulated number of storms in 2020 and 2021 is ",avg20_21,". This corresponds nicely to the value obtained in 2a), which was ",expected_storms20_21$value))

#Standard deviation
sd20_21 <- sqrt(var(NHPPnumbers20_21))
print(paste0("The simulated standard deviation of the total number of storms in 2020 and 2021 combined is ",sd20_21,". This corresponds nicely to the value obtained in 2a), where the standard deviation was found to be ",sd))



#2d)
claim_func <- function(t){ #the claim amount as defined in equation 5 in the sheet. 
  return(10*exp(5*t/100)) 
}

claims <- function(a,b,intensmax, intens){ #A function that simulates total claim size.
  times_storms_occur <- simtNHPP(a = a, b = b, intensitymax = intensmax, intensityfunc = intens)
  amount_claim <- vector(mode="numeric", length = length(times_storms_occur))
  for(t in 1:length(times_storms_occur)){
    amount_claim[t] <- claim_func(times_storms_occur[t])
  }
  total_claims <- sum(amount_claim)
  return(total_claims)
}

claimsmoney <- claims(a = 0, b = 1,intensmax = intensitymax, intens = intensity)

Nsim <- 10000
claims_Nsim <- vector(mode="numeric",length(Nsim))
for(i in 1:Nsim) {
  claims_Nsim[i] <- claims(a = 0, b = 1, intensmax = intensitymax, intens = intensity)
}

#Mean
mean_claims <- mean(claims_Nsim)
#SD
sd_claims <- sd(claims_Nsim)

#95% confidence interval
lower_side <- mean_claims - qnorm(0.975)*sd_claims/sqrt(length(claims_Nsim))
upper_side <- mean_claims + qnorm(0.975)*sd_claims/sqrt(length(claims_Nsim))
print(paste0("95% confidence interval for the mean u: [",lower_side,", ",upper_side,"]"))


#Be 97.5% certain to cover 
certanty <- mean_claims + qnorm(0.9875)*sd_claims/sqrt(length(claims_Nsim))
print(paste0("Must set aside at least: ",certanty, " Million NOK"))

#2e)
#Number of storms
y_func <- function(N,t) { #The function Y (equation 7), from the sheet. 
  ifelse(N==0,0,sum(claim_func(t)))
}

N <- 10000
costs <- vector(mode="numeric",length=N)
for(n in 1:N) {
  storms_time_occurence <- simtNHPP(a=0,b=1,intensitymax = 60,intensityfunc = intensity)
  l <- length(storms_time_occurence)
  costs[n] <- y_func(l,storms_time_occurence)
}


#Values from 2d:
mean_claims
sd_claims
#Values using MC
mean(costs) 
sd(costs)

#Compared to the values obtained in d), the mean is slightly lower and the standard deviation is slightly higher. 
#This may be due to the fact that these data are not analytical. Therefore, one can conclude that they are similar. 
#This is also reflected by the confidence interval by found in 2d).

#2f)
M <- 100
variances <- vector(mode="numeric",length=M)
means <- vector(mode="numeric",length=M)
for(m in 1:M) {
  N <- 100
  costs <- vector(mode="numeric",length=N)
  for(n in 1:N) {
    storms_time_occurence <- simtNHPP(a=0,b=1,intensitymax = 60,intensityfunc = intensity)
    l <- length(storms_time_occurence)
    costs[n] <- y_func(l,storms_time_occurence)
  }
  variances[m] <- var(costs)
  means[m] <- mean(costs)
}



finalvar <- mean(variances) + var(means)









##################
#Problem 3
#Bootstrapping for a linear regression model

#Load data into environement
load("prob23.dat")

#Run a single regression on these data
lm.obj <- lm(y ~ x1 + x2 + x3 + x4 + x5, data = df)

#Extract R^2 
Rsquared <- summary(lm.obj)$r.squared # original parameter estimate

#How to  resample complete rows of data sets
df.star <- df[sample.int(n,replace = TRUE),]

#3a)
#Function that computes B bootsrap samples of R^2 for the given dataset
n <- length(df[,1]) #I checked that n is the same for df[,1],df[,2],df[,3],df[,4],df[,5],df[,6], meaning all colums have same amount of rows, so safe to use length(df[,1]) to find n

rsquared_bs <- function(B) {
  rsquares <- vector(mode="numeric",length = B)
  for(i in 1:B) {
    #Resample data
    df.star <- df[sample.int(n=n, replace = TRUE),]
    #Run a single regression on these data
    lm.obj <- lm(y ~ x1 + x2 + x3 + x4 + x5, data = df.star)
    #Extract R square
    Rsquared <- summary(lm.obj)$r.squared
    rsquares[i] <- Rsquared
  }
  return(rsquares)
}

B <- 5000
bootstrapped <- rsquared_bs(B)

#standard deviaton
stddev <- sd(bootstrapped)
#Check the standard deviation to confirm the hint
print(paste0("Checking if standard deviation is approx 0.017. Simulated value: ",stddev))

#Mean of the data
meanboot <- mean(bootstrapped)

#Create histogram
hist(bootstrapped,prob=TRUE,main=paste0("Historgram of B=",B,"\nbootstrap samples of R^2"))
curve(dnorm(x, mean=meanboot, sd=stddev), 
      col="darkblue", lwd=2, add=TRUE, yaxt="n")

#Plotting a normally distributed curve on top of the histogram with mean and standard deviation as calculated by the simulations
#We see that the bootstrap samples are well approximated by a normal distribution. 


#3b)
#Bootstrap estimate of the estimated Bias
bias <- meanboot - Rsquared 

#Bootstrap estimate of standard deviation of the estimated R^2
#Implementing equation (8.1) from Rizzo, p. 215. [Source 2, p. 215]

thesum <- 0
for(B in 1:length(bootstrapped)) {
  val <- (bootstrapped[B] - meanboot)^2
  thesum <- thesum + val
}
boot_sd <- sqrt(thesum/(B-1))
print(paste0("The bootstrap estimate of the standard deviation of the estimated r^2 is: ",boot_sd))

#Standard normal bootstrap interval 
error <- qnorm(0.99)*boot_sd
lowerN <- Rsquared - error
upperN <- Rsquared + error
print(paste0("99% confidence interval for R^2 using standard normal bootstrapping: [",lowerN,",",upperN,"]"))


#Percentile bootstrap interval
#First find alpha
#(1-alpha)*100%=99%
#alpha = 1-99%/100% = 0.01
#so percentile interval is [values of alpha/2, values of 1-alpha/2]
alpha <- 0.01
p_boot_sd_low <- quantile(bootstrapped,probs=alpha/2)
p_boot_sd_high <- quantile(bootstrapped,probs=1-alpha/2)
print(paste0("99% confidence interval for R^2 using percentile bootstrap interval: [",p_boot_sd_low[[1]],",",p_boot_sd_high[[1]],"]"))

#Comment
#The standard normal bootstrap interval uses the assumption that the estimator is unbiased.
#Also, the boot_sd used is a estimated value, but is treated as a known parameter.
#The percentile interval uses the empirical distribution as the reference distribution. The quantiles
#of the empirical distribution are estimators of the sampling dist. These random quantiles may match
#the true distribution when the dist of the estimator is not normal. Further, percentile interval has some
#theoretical advantages as proven by Efron and Tibishirani in [source 1, 13.3].




##############
#Problem 4

#4b)
invtransform <- function(n,l,mu,sigma){ #Inverse transfrom for drawing random numbers with density equation (8) from the sheet
  U <- runif(n, min = 0, max = 1)
  R <- qnorm(pnorm(l,mean = mu, sd = sigma) + U*(1-pnorm(l, mean = mu, sd = sigma)),mean = mu, sd = sigma)
  return(R)
}

l <- 1.8
mu <- 2.2
sigma <- 1.5
n <- 1000
r<- invtransform(n = n, l = l, mu = mu, sigma = sigma)
r
hist(r,probability = TRUE) 

#mean
mean_r <- mean(r)
sd_r <- sd(r)
mean_r
sd_r

#Find mean and sd using integrate
Rfunction <- function(r,l,mu,sigma) {
  if(r < l) {
    return(0)
  }
  else{
    top <- exp(-((r-mu)^2)/(2*sigma^2))
    bottom <- sqrt(2*pi)*sigma*(1-pnorm(l, mean = mu, sd = sigma))
    return(top/bottom)
  }
}


fr_mean <- function(r) {
  top <- r*exp(-((r-2.2)^2)/(2*1.5^2))
  bottom <- sqrt(2*pi)*1.5*(1-pnorm(1.8,mean=2.2,sd = 1.5))
  return((top/bottom)*(r>=1.8))
} 

analytic_mean <- integrate(fr_mean,-Inf,Inf)
print(paste0("The mean using integrate is: ",analytic_mean))

fr_var <- function(r) {
  top <- r*r*exp(-((r-2.2)^2)/(2*1.5^2))
  bottom <- sqrt(2*pi)*1.5*(1-pnorm(1.8,mean=2.2,sd = 1.5))
  return((top/bottom)*(r>=1.8))
} 

analytic_var <- integrate(fr_var,-Inf,Inf)$value - (analytic_mean$value)^2 
print(paste0("The standard deviation using integramte is: ",sqrt(analytic_var)))

#4c)
fScondR <- function(s,n){ #f S conditional R. 
  R <- invtransform(n,l = 1.8, mu = 2.2, sigma= 1.5)
  fScondR_s <- (runif(n,min=0.5/R,max=0.8/R))*((0.5/R) <= s & s <= (0.8/R))
  return(fScondR_s)
}

#Write a function that simulates n random vectors with density 9
#Simulate R first, then from S cond R.

#Inv transform to find S
n = 1000
R <- invtransform(n,l = 1.8, mu = 2.2, sigma= 1.5)
U <- runif(n)
S_R <- U*(3/(10*R))+(1/(2*R)) #This steps to reach this expression is shown by hand in the accompanying pdf handed in. 


RSmatrix <- matrix(nrow = n, ncol = 2)
RSmatrix[,1] <- R
RSmatrix[,2] <- S_R
colMeans(RSmatrix)
cor(RSmatrix[,1],RSmatrix[,2])
cov(RSmatrix[,1],RSmatrix[,2])

plot(R,S_R,ylab="S",main = )



#4d)

S_R <- sort(S_R)

marg_S  <- vector(mode="numeric",length = length(S_R))
for(s in 1:length(S_R)) {
  current_s <- S_R[s]
  f_RS <- function(r) {
      top <- exp(-((r-2.2)^2)/(2*1.5^2))
      bottom <- sqrt(2*pi)*1.5*(1-pnorm(1.8, mean = 2.2, sd = 1.5))
      other_term <- (10*r)/3
      return((top*other_term/bottom)*(0.5/r <= current_s & 0.8/r >= current_s))
  }
  marg_S[s] <- integrate(f_RS,1.8,Inf)$value
}


hist(S_R,probability=TRUE) #S conditional of R is uniform with varying R. 
lines(S_R,marg_S,col="red") #When plotting S marginal, we see that it is not uniform.


#4e) 
source("seir_model.R")

#Antithetic
N <- 100
U <- runif(N)
V <- 1 - U 
Y <- runif(N)
W <- 1-Y
R <- qnorm(pnorm(1.8, mean = 2.2, sd = 1.5) + U*(1-pnorm(1.8, mean = 2.2, sd = 1.5)),mean = 2.2, sd = 1.5)
R_ant <- qnorm(pnorm(1.8, mean = 2.2, sd = 1.5) + V*(1-pnorm(1.8, mean = 2.2, sd = 1.5)),mean = 2.2, sd = 1.5)
allR <- c(R,R_ant) 
ScondR_reg <- Y*(3/(10*R))+(1/(2*R))  
ScondR_ant <- (W)*(3/(10*R_ant))+(1/(2*R_ant))
allS <- c(ScondR_reg,ScondR_ant)

#Using antithetic
res = 10^6 #resolution
maxBeds <- vector(mode="numeric",length = length(allR))
startTime <- Sys.time()
for(t in 1:length(allR)) {
  f <- seir_sim(R0max = allR[t], soc_dist_Rfac = allS[t])
  days <- which(grepl("2021",f$dates))
  days2021 <- which(f$dates=="2021-12-31") - which(f$dates=="2021-01-01") +1 #+1 to get entire year
  numberOfBeds <- vector(mode = "numeric",length = days2021)
  for(d in days) {
    numberOfBeds[d-366] <- res*f$ode[d,"HC"] + res*f$ode[d,"CC"] + res*f$ode[d,"HH"]
  }
  maxBeds[t] <- max(numberOfBeds)
}
endTime <- Sys.time()
runtime <- endTime - startTime
print(paste0("Runtime antithetic: ",runtime))
mean(maxBeds)
var(maxBeds)
sd(maxBeds)
mean(maxBeds^2)

#################
#Using regular variables
N <- 200
U <- runif(N)
Y <- runif(N)
R <- qnorm(pnorm(1.8, mean = 2.2, sd = 1.5) + U*(1-pnorm(1.8, mean = 2.2, sd = 1.5)),mean = 2.2, sd = 1.5)
S <- Y*(3/(10*R))+(1/(2*R))
maxBeds <- vector(mode="numeric",length = length(R))
startTime <- Sys.time()
for(t in 1:length(R)) {
  f <- seir_sim(R0max = R[t], soc_dist_Rfac = S[t])
  days <- which(grepl("2021",f$dates))
  days2021 <- which(f$dates=="2021-12-31") - which(f$dates=="2021-01-01") +1 #+1 to get entire year
  numberOfBeds <- vector(mode = "numeric",length = days2021)
  for(d in days) {
    numberOfBeds[d-366] <- res*f$ode[d,"HC"] + res*f$ode[d,"CC"] + res*f$ode[d,"HH"]
  }
  maxBeds[t] <- max(numberOfBeds)
}
endTime <- Sys.time()
runtime <- endTime - startTime
print(paste0("Runtime regular: ",runtime))
mean(maxBeds)
var(maxBeds)
sd(maxBeds)
mean(maxBeds^2)

#Even though runtime and mean is fairly similar, we see that the variance in higher with antithetic.
#The function is not monotonic, thus antithetic does not reduce variance. 


###########
#4f)

#Using regular variables since antithetic had higher variance
N <- 200
U <- runif(N)
Y <- runif(N)
R <- qnorm(pnorm(1.8, mean = 2.2, sd = 1.5) + U*(1-pnorm(1.8, mean = 2.2, sd = 1.5)),mean = 2.2, sd = 1.5)
ScondR_reg <- Y*(3/(10*R))+(1/(2*R))  
allR <- R
allS <- ScondR_reg


days2021 <- which(f$dates=="2021-12-31") - which(f$dates=="2021-01-01") +1 #+1 to get 365 days

daily2021BedsMatrix <- matrix(nrow=length(allR),ncol = days2021)
for(t in 1:length(allR)) {
  f <- seir_sim(R0max = allR[t], soc_dist_Rfac = allS[t])
  days <- which(grepl("2021",f$dates))
  days2021 <- which(f$dates=="2021-12-31") - which(f$dates=="2021-01-01") +1 #+1 to get entire year
  numberOfBeds <- vector(mode = "numeric",length = days2021)
  for(d in days) {
    daily2021BedsMatrix[t,d-366] <- res*f$ode[d,"HC"] + res*f$ode[d,"CC"] + res*f$ode[d,"HH"]
  }
}

dailyMeans <- vector(mode="numeric",length = length(daily2021BedsMatrix[1,]))
dailyUpper <- vector(mode="numeric",length = length(daily2021BedsMatrix[1,]))
dailyMedians <- vector(mode="numeric",length = length(daily2021BedsMatrix[1,]))
dailyLower <- vector(mode="numeric",length = length(daily2021BedsMatrix[1,]))
for(day in 1:length(daily2021BedsMatrix[1,])) {
  dailyMeans[day] <- mean(daily2021BedsMatrix[,day])
  dailyUpper[day] <- quantile(daily2021BedsMatrix[,day], probs=0.95)
  dailyLower[day] <- quantile(daily2021BedsMatrix[,day],probs=0.05)
  dailyMedians[day] <- median(daily2021BedsMatrix[,day])
}

plot(dailyMeans,col="red",ylim=c(min(dailyLower)-10,max(dailyUpper)+30),type="l",main="Uncertanty of daily beds",xlab = "Days",ylab="Beds")
lines(dailyMedians,col="blue")
lines(dailyUpper,col="green")
lines(dailyLower,col="black")
legend("top",lty=1,col=c("red","blue"),legend = c("Mean","Median"),bg="transparent",bty="n")
legend("topright",lty=1,col=c("green","black"),legend = c("95% quantile","5% quantile"),bg="transparent",bty="n")

#The plot shows that there are large differences in the quantiles and the mean and median are also different. 
#This thus illustrates the uncertainty. 
  