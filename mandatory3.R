#Mandatory 3
#Student number 236227
#Aleksander B. Jakobsen


rm(list=ls())

#Implement libraries
library(mvtnorm)
library(coda)
library(boot)


#Note: Check that datasets are in same directory as .R-file,
#and that working directory is set to source file location. 

#################################
#Problem 1
#1a)

#Distribution given in the problem            
p_x_a <- function(x) {
  return(((2^(1/4)*gamma(3/4))/pi)*exp((-x^4)/2))
}

#Implementing independent MH sampler because 
#the exercise asks to generate independent sample.
#Not using RWMH to avoid auto-correlation.
#Note: Sigma was identified as 0.8 (yielding highest accept rate),
#this sigma was identified by trying various values. 

sigma <- 0.8 
n <- 10000
out <- numeric(n) #Allocate memory
out[1] <- 0.0 
#Old importance weight
weight.old <- p_x_a(out[1])/dnorm(out[1],sd=sigma)
#Accept counter
Nacc <- 0
for(i in 2:n) {
  #proposal, independent of past
  thetaStar <- rnorm(1,sd=sigma)
  # new importance weight
  weight.star <- p_x_a(thetaStar)/dnorm(thetaStar,sd=sigma)
  #accept prob
  alpha <- min(1,weight.star/weight.old)
  #Accept/reject
  if(runif(1)<alpha){
    out[i] <- thetaStar
    weight.old <- weight.star
    Nacc <- Nacc+1
  } else {
    out[i] <- out[i-1]
  }
}

print(paste0("Accept rate : ",Nacc/(n-1))) #We see that sigma gives high accept prob, note: n-1 becaus of samples
hist(out,probability = TRUE,breaks=20,main="p(x)",xlab="x") 
xg <- seq(from=-2,to=2,length.out = 10000)
lines(xg,p_x_a(xg),col="red")



#1b)
#Implementing inverse transform method, so that for/while loops are avoided: 
p_x_b_function <- function(x) {
  return(2*x*exp(-x^2))
}

n <- 10000
U <- runif(n)
p_x_b <- numeric(length=n)
p_x_b <- sqrt(-log(1-U))
hist(p_x_b,probability = TRUE,main="p(x)",xlab = "x") 
xg <- seq(from=0,to=3,length.out = 1000)
lines(xg,p_x_b_function(xg),col="blue")
#We see that produced histogram fits well with the probability density.



#1c)
#Evaluate integral from 0 to inf using importance sampling:

#The function in the integral
p_x_c <- function(x){
  return(exp(sqrt(x))*exp(-20*(x-4)^2))
}
#See pdf note to find normal dist, using this as importance function when evaluating.
n <- 10000
x <- rnorm(n,mean=4,sd=sqrt(1/40))

for(i in 1:length(x)) { #Only interested in values above 0
  if(x[i] < 0) {
    x[i] <- 0 
  }
}


#Evaluate the integral for every value
evaluated <- p_x_c(x)/dnorm(x,4,sqrt(1/40))
#Mean
evaluated <-  mean(evaluated)

print(paste0("The integral evaluated is: ",evaluated," using importance sampling."))




#######################################  
#Problem 2
#Implementing the function given in the exercise. Takes in two dimensional theta
log_g <- function(theta){
  return(((-(theta[1])^2/2 ))-((theta[2]-theta[1]^2)^2/2)# log-target (i.e. log g)
  )
}


#Two dimensional implementation of RWMH that I have written myself. 
#Using this since it is the most used algorithm to sample multi-dimensional distributions,
#as mentioned in lecture note on Markov Chain Monte Carlo part 1.
#Since no special requirements or information is given,
#I do not implement any of the special cases of the random walk. 

#Note, the rmvnorm package contains cholesky decomposition for efficiency
RWMH_2D <- function(lprob, #logprob density
                    Sigma, 
                    thetaInit=c(0.0,0.0), #Initial values
                    n=10000){ #iterations
  
  #Allocate output memory
  result <- matrix(0.0,n,2)
  #Set initial condition thetaInit
  result[1,] <- thetaInit
  
  #previous log probability
  old <- lprob(thetaInit)
  
  #covariance equal to identity matrix
  idcovar <- diag(Sigma,2,2)
  
  #accept counter
  c <- 0
  # main iteration
  for(i in 2:n){
    
    #new proposal
    thetaStar <- result[(i-1),] + rmvnorm(1,thetaInit,idcovar)
    #log of new proposal
    now <- lprob(thetaStar)
    
    #acceptance probability
    alpha <- exp(min(0.0,now-old))
    
    #acceptance-rejection step
    if(runif(1)<alpha && is.finite(alpha)){
      #accept this value
      result[i,] <- thetaStar
      #current log value becomes old log value in next iteration.
      old <- now #This saves times. Dont have to re-calculate old value in next iteration.
      #update counter
      c <- c+1
    } 
    else {
      #rejection
      result[i,] <- result[(i-1),]      
    }
  }
  print(paste0("RWMH done, accept rate : ",c/(n-1))) #note n-1 since sample size
  return(result) 
  
}
output<-RWMH_2D(log_g,Sigma=2,thetaInit = c(0.0,0.0),n=100000) #Sigma 2 to get accpetance rate 
output <- output[10000:length(output[,1]),] #Burn in removal. This comes from geweke plot analysis done below.
#Plot to see smile-shaped target
plot(output,xlab="Theta_1",ylab="Theta_2",main = "Theta",col="lemonchiffon4")
hist(output)
#Historgrams to visualize marginal distributions
hist(output[,1],probability = TRUE,xlab="Theta_1",main="Theta_1 marginal distribution",col="lightcoral") #theta1
curve(dnorm(x, mean=0, sd=1), col="darkblue", lwd=2, add=TRUE, yaxt="n")
#We see from the histogram and curve that the marginal theta_1 is N(0,1),
#and can therefore assume that the implementation is correct. 


#Evaluate effective sample size, ESS
#Using coda, source [1]. 
#Evaluating for both columns
ESS1 <- coda::effectiveSize(output[,1])
ESS2 <- coda::effectiveSize(output[,2])

print(paste0("We see that for theta1, effective sampling size is ",round(ESS1,digits = 2)," > 1000."))
print(paste0("We see that for theta2, effective sampling size is ",round(ESS2,digits = 2)," > 1000."))

mcmc_output <- mcmc(output)
geweke.plot(mcmc_output) #From the geweke plot, we see that most
#iterations have Z-scores within accetable limits, meaning that it 
#identifies same distribution among various values of the set.
#Burn in is therefore OK. 


###########################
#Problem 3 
#Make sure that you have 
#Start by implementing code given in problem description


#load the data set
df <- data.frame(read.table("logistic_regression_data.txt"))
x <- df$x
y <- df$y

#function returning a log-posterior kernel for theta=(alpha,beta)
logistic.lp <- function(theta){
  alpha <- theta[1]
  beta <- theta[2]
  #log-likelihood
  Eeta <- exp(alpha+beta*x)
  p <- Eeta/(1.0+Eeta)
  log.like <- sum(dbinom(y,size=1,prob=p,log=TRUE))
  
  #priors
  log.prior <- dnorm(alpha,sd=10,log=TRUE) + dnorm(beta,sd=10,log=TRUE)
  
  #log-posterior kernel
  return(log.like+log.prior)
}


#3a)

#Implementing independent metropolis hastings (IMH) algorithm
IMH <- function(logprob,
                thetaInit,
                sigma=diag(2),
                n=10000,
                scale=1) {
  #Allocate memory
  result <- matrix(0,n,2)
  #Initial conditions
  result[1,] <- thetaInit
  #previous log probability 
  old <- logprob(thetaInit) - dmvnorm(thetaInit,mean=thetaInit,sigma=sigma,log=TRUE)
  #counter
  c <- 0
  for(i in 2:n) {
    #new proposal
    thetaStar <- rmvnorm(1,thetaInit,scale*sigma)
    #weight of new proposal as log
    new <- logprob(thetaStar) - dmvnorm(thetaStar,mean=thetaInit,sigma=sigma,log=TRUE)
    alpha <- exp(min(1.0,new-old))
    #acceptance rejection step
    if(runif(1)<alpha && is.finite(alpha)) {
      #accept
      result[i,] <- thetaStar
      #new becomes old in next iteration
      old <- new 
      #update counter 
      c <- c+1
    } 
    else {
      result[i,] <- result[i-1,]
    }
  }
  print(paste0("Accept rate: ",c/(n-1))) #n-1 because samples 
  return(result)
}

sigma <- matrix(c(0.00653,-0.00058,-0.00058,0.01689),2,2)

values <- IMH(logistic.lp,thetaInit = c(-0.102,1.993),sigma = sigma, n = 10000, scale = 0.6) #scale is tuning parameter. Is set to 0.6 because this gave a high acceptance rate as well as good effective sample sizes. 
values <- values[2:length(values[,1]),] #Note, only removing first values
#as burn in because this is independent Metropolis Hasting.
#The independence between values means that it is only necessary to burn the first value which is set as input theta. 
plot(values,xlab="Alpha",ylab = "Beta",col = "darksalmon") 

hist(values[,2],probability=TRUE,main="Marginal beta",xlab="Beta",col="cyan4")
hist(values[,1],probability=TRUE,main="Marginal alpha",xlab = "Alpha",col="gold")

print(paste0("Alpha mean posterior:", mean(values[,1])))
print(paste0("Beta mean posterior: ",mean(values[,2])))

#Evaluate effective sample size
ESSalpha <- coda::effectiveSize(values[,1])
ESSbeta <- coda::effectiveSize(values[,2])

print(paste0("Effective sample size alpha: ",ESSalpha))
print(paste0("Effective sample size beta: ",ESSbeta))

#3b)
#Predictions of a random variable m(x*) = exp(alpha+betax*)/1+exp(alpha+betax*)
#for some new observation y* with associated covariate x*. (alpha,beta) are distributed
#according to posterior dist from a). 

#Note that alphas are values[,1] and betas are values[,2] from 3a).
#alphas <- values[,1]
#betas <- values[,1]

#I interpret this problem as if I am asked to calculate m(x*) for 
#x* on the given interval [-5,5] for every alpha and beta value. 
#I therefore make a function which return the  desired quantiles
#and median. 


getinfo <- function(x_star) {
  m_x_star <- exp(values[,1]+values[,2]*x_star)/(1+exp(values[,1]+values[,2]*x_star))
  quantile0.05 <- quantile(m_x_star, 0.05)
  quantile0.95 <- quantile(m_x_star, 0.95)
  median <- median(m_x_star)
  return(c(median,quantile0.05,quantile0.95))
}

#create x* values on the interval [-5,5] as given in the exercise
x_star <- seq(from = -5, to = 5, by = 0.01)

#Now I have a function and the x_stars, and need to iterate over every x_star
#I therefore create a matrix to store information, where rows are information
#about media and the quantiles, and columns is for each x value
infomatrix <- matrix(0.0,nrow=3,ncol=length(x_star))
#iterate and find info for each x
for(i in 1:length(x_star)) {
  infomatrix[,i] <- getinfo(x_star[i])
}

plot(x_star,infomatrix[1,],type="l",col="cadetblue",main="m(x*)",xlab = "x*")
#Note, this plot will vary with various delta (scale) values.
lines(x_star,infomatrix[2,],type="l",col="darkolivegreen")
lines(x_star,infomatrix[3,],type="l",col="firebrick")
legend("topleft",lty=1,col=c("cadetblue","darkolivegreen","firebrick"),legend = c("Median","0.05 Quantile","0.95 Quantile"),bg="transparent",bty="n")

#3c)
get0.99quantile2 <- function(x_star) {
  m_x_star <- exp(values[,1]+values[,2]*x_star)/(1+exp(values[,1]+values[,2]*x_star))
  quantile0.99 <- quantile(m_x_star,0.01)
  return(c(x_star,quantile0.99))
}


matrix0.99 <- matrix(0.0,nrow=2,ncol=length(x_star))
for(i in 1:length(x_star)) {
  matrix0.99[,i] <- get0.99quantile2(x_star[i])
}

v <- matrix0.99[2,]

c <- matrix0.99[1,Position(function(x) x > 0.8, v)] #source [2]

print(paste0("For x* values from ",c, " and above, we are at least 99% certain that m(x*)>0.8."))




##################
#Problem 4

#4b) 
#Load the data set we are considering using the
#comands given in problem statement.
df <- data.frame(read.table(file="linear_regression_data.txt"))
x <- df$x
y <- df$y

#Gibbs sampler (3 blocks) targeting the joint posterior 
#p(theta|y) associated with the data set, using updates
#according to alpha|beta,tau,y , beta|alpha,tau,y and tau|alpha,beta,y.
#defining function that takes data input and theta, as well as iteration number. 
gibbs_3_block <- function(x_data,y_data,theta=c(1,-1.2,1),n=10000) {
  #Note: theta[1] = alpha, theta[2]=beta, theta[3] = tau
  #Allocate memory for output matrix
  result <- matrix(0.0,3*n+1,3)
  #Initial state
  result[1,] <-theta
  #counter
  c <- 2
  for(i in 1:n) {
    #block 1: normal distribution alpha
    b_1 <- -theta[3]*sum(theta[2]*x_data - y_data)
    c_1 <- -(theta[3]*length(x_data)/2 + 1/200)
    theta[1] <- rnorm(1,-b_1/(2*c_1),-1/(2*c_1))
    result[c,] <- theta
    c <- c+1
    #block 2: normal distribution beta
    b_2 <- -theta[3]*sum(x_data*(theta[1]-y_data))
    c_2 <- -((theta[3]/2)*(sum(x_data^2) + 1/200))
    theta[2] <- rnorm(1,-b_2/(2*c_2), (-1/(2*c_2)))
    result[c,] <- theta
    c <- c+1
    #block 3: gamma dist tau (note: shape and scale as given in problem formulation)
    g_shape <- length(x_data)/2 + 1
    g_scale <- 1/((1/2)*sum((y_data-theta[1]-x_data*theta[2])^2) + 1)
    theta[3] <- rgamma(1,shape = g_shape, scale = g_scale)
    result[c,] <- theta
    c <- c+1
    
  }
  return(result)
}

#4c) 
#Run sampler function from b) for at least 10000 iters.

gibbs_results <- gibbs_3_block(x_data = x, y_data = y)
#Visualisations 
plot(gibbs_results[,1],gibbs_results[,2],type = "l",xlab="Alpha",ylab="Beta",col = "darkgoldenrod1") 
plot(gibbs_results[,2],gibbs_results[,3],type = "l",xlab="Beta",ylab="Tau", col = "darkgoldenrod2") 
plot(gibbs_results[,1],gibbs_results[,3],type = "l",xlab="Alpha",ylab="Tau",col="darkgoldenrod3") 
#We see from the plots that burn-in needs to be romve

#Remove burn in as needed:
burn <- gibbs_results[500:length(gibbs_results[,1]),] #Burn in of 500 removed after seeing initial traceplots and Geweke plot below.  
plot(burn[,1],burn[,2],type="l",xlab="Alpha",ylab="Beta",col="coral") 
plot(burn[,2],burn[,3],type="l",xlab = "Beta",ylab="Tau",col="coral1")
plot(burn[,1],burn[,3],type="l",xlab="Alpha",ylab="Tau",col="coral2") 

#Compare with classical linear regression
colMeans(burn)
lm.out <- lm(y~x,data=df)
lm.out
#We see that the posterior means are approximately the same as classical regresison

burn_ESS <- coda::effectiveSize(burn)
burn_ESS #We see that we have at least 1000 effective sample size for alpha, beta and tau

burn_mcmc <- mcmc(burn)

#Traceplot
burn_mcmc_alpha <- burn_mcmc[,1]
burn_mcmc_beta <- burn_mcmc[,2]
burn_mcmc_tau <- burn_mcmc[,3]
traceplot(burn_mcmc_alpha,col="deepskyblue",ylab = "Alpha",main="Traceplot of alpha")
traceplot(burn_mcmc_beta,col="firebrick",ylab="Beta",main="Traceplot of beta")
traceplot(burn_mcmc_tau,col="darkolivegreen1",ylab="Tau",main="Traceplot of tau")
#We see from all the traceplots that the values are converged within
#some finite range, and can therefore be satisfied with this. 

#Geweke plots 
varnames(burn_mcmc) <- c("Alpha", "Beta", "Tau") 
geweke.plot(burn_mcmc)
#We see that in most runcases there is good values,
#in some cases some outliers, but nothing too dramatic.


#4d)

gibbs_2_block <- function(x_data, y_data, theta=c(1.0,-1.2,1.0),n=10000) {
  #Allocate memory for output matrix
  result <- matrix(0.0,2*n+1,3)
  #Initial state
  result[1,] <- theta
  counter<-2
  #Same as before, but now block1 is [alpha,beta]
  for(i in 1:n) {
    #C matrix
    c_matrix <- matrix(c(-length(x_data)*theta[3]-0.01, -theta[3]*sum(x_data), -theta[3]*sum(x_data), -theta[3]*sum(x_data^2)-0.01), nrow=2, ncol = 2)
    #B vector
    b_vector <- c(theta[3]*sum(y_data), theta[3]*sum(y_data*x_data))
    #Alpha, beta
    theta[1:2] <- rmvnorm(1,-solve(c_matrix)%*%b_vector, -solve(c_matrix))
    result[counter,] <-theta
    counter <- counter + 1
    #tau
    theta[3] <- rgamma(1,shape = length(x_data)/2+1, scale = 1/((1/2)*sum((y_data-theta[1]-x_data*theta[2])^2) + 1))
    result[counter,] <- theta
    counter <- counter + 1
  }
  return(result) 
  }



gibbs_results2 <- gibbs_2_block(x_data = x, y_data = y)
gibbs_results2 <- gibbs_results2[501:length(gibbs_results2[,1]),] #Burn in of 500 to keep similar as 3block. Geweke plot below also showes that this is acceptable. 
hist(gibbs_results2)

#Visualization
plot(gibbs_results2[,1],gibbs_results2[,2],type="l",col="cornflowerblue",xlab="Alpha",ylab = "Beta",main="Gibbs sampler, 2 block") 
#Check colMeans to see that output is somewhat expected.
colMeans(gibbs_results2)
gibbs2_mcmc <- mcmc(gibbs_results2)
#Geweke plot
varnames(gibbs2_mcmc) <- c("Alpha", "Beta", "Tau") 
geweke.plot(gibbs2_mcmc) 

#effective sample sizes
ess_2_block <- coda::effectiveSize(gibbs_results2)
ess_2_block 
#Observe that effective sample sizes has increased. 
#This is due to the multinormal distribution for alpha and beta
#using the covariance matrix, whereas the three block gibbs function
#sampled values with dependency. 
#Note that the result is reflecting lecture note "Markoc chain Monte
#Carlo, part 3", where it is stated that more blocks -> more conditional
#distribution. Since we have reduced block, we have less dependency. 

#4e) 
#Does the analysis of data suggest that alpha+beta is not equal to 0?
# 3 block gibbs data 
aplusb_3block <- burn[,1] + burn[,2] 
hist(aplusb_3block)
abline(v = c(quantile(aplusb_3block, 1-0.995), quantile(aplusb_3block, 0.995)), col="red")
#We see that the three block alpha+beta does not encapture 0 when considering
#a 99% confidence interval. Therefore, the analysis here suggest that 
#we may discard alpha+beta=0, and indeed confirm that alpha+beta not equal to 0.

#2 block:
aplusb_2block <- gibbs_results2[,1] + gibbs_results2[,2]
hist(aplusb_2block)
abline(v = c(quantile(aplusb_2block, 1-0.995), quantile(aplusb_2block, 0.995)), col="red")
#However, here we see that a 99% confidence interval encaptures 0. 
#Therefore, we cannot reject that alpha+beta=0.

#In light of both these tests, it is inconclusive. 



#################################
#Problem 5: Bootstrapping
#Same linear regression situation as in previous problem
#including the same data set, but bootstrap approach.

linreg <- function(d,i){ #Function for bootstrap
  lm(d[i,1]~d[i,2],data = df)$coefficients
}
#Bootstrap
bootstrap <- boot(data = df, statistic = linreg, R=5000)
plot(bootstrap$t[,1],bootstrap$t[,2],type="l",xlab="Alpha",ylab="Beta",main="Bootstrap",col="darkcyan") 
#We see that this plot is similar to the plot in exercise 4. 

#Claim 1
sigmas <- 1/sqrt(rgamma(10000, shape = length(x)/2 + 1, scale = 1/((1/2)*sum((y-mean(bootstrap$t[,1]) - mean(bootstrap$t[,2])*x)^2)+1)))
hist(sigmas,breaks=40,main="Sigmas",xlab="Sigma",col="darkslategray3") 
#As we can see from the histogram, sigma does in fact,
#go over 1.0, and is therefore not strictly smaller than 1.0.

#Claim 2
alpha <- bootstrap$t[,1]
hist(alpha,main="Alpha distribution",col = "dodgerblue3") #We see that alpha is not equal to 0. 
#This may be illustrated with the 99% confidence interval:
#Note: 99% CI gives alpha = 0.005, which gives lower value 0.005 and upper value 1-0.005 = 0.995.
abline(v=c(quantile(alpha,0.995),quantile(alpha,0.005)) ,col="firebrick")
#We see that 99% confidence interval does not include alpha = 0.
#And we therefore reject the claim that alpha is equal to 0.
 
#Claim 3
#Same argumentation as claim 2:
aplusb <- alpha + bootstrap$t[,2]
hist(aplusb,main="Distribution of alpha+beta",xlab="alpha + beta",col="forestgreen") 
#Confidene intervals 
abline(v = c(quantile(aplusb,0.995),quantile(aplusb,0.005)),col="firebrick")
#The confidence intervals include the value 0.
#As a consequence, we can not discard claim 3. 




######################
#Problem 6

source("seir_mcmc_funs.R")
#Implementing almost identical function as in problem 2, with minor tweaks such as cholesky
RWMH2D_covid  <- function(logprob, #Notice log prob density kernel
                          Sigma=diag(2),
                          thetaInit = c(3.18,0.20),
                          n=1000) {
  #Allocate memory
  result <- matrix(0.0,n,2)
  #Append initial values
  result[1,] <- thetaInit
  #old logprob
  oldlog <- logprob(thetaInit)
  #counter
  c <- 0
  #cholesky for sampling efficiency
  cholesky = t(chol(Sigma)) 
  
  for(i in 2:n) {
    thetaStar <- result[i-1,] + cholesky%*%rnorm(2)
    nowlog <- logprob(thetaStar) 
    alpha <- exp(min(0.0,nowlog-oldlog))
    #acceptance/rejection step
    if(runif(1)<alpha && is.finite(alpha)) {
      #the value is accepted
      result[i,] <- thetaStar
      #current value becomes previous value in next step
      oldlog <- nowlog
      c <- c+1
    }
    else{
      result[i,] <- result[i-1,]
    }
  }
  print(paste0("Acceptance rate: ",c/(n-1))) #Note n-1 due to sample size.
  return(result) #Notice: no burn in removed. Tuning process will be described in steps. 
}

output <- RWMH2D_covid(seir.lp,Sigma=matrix(c(3e-4,0,0,5e-7),nrow = 2, ncol = 2),thetaInit = c(3.18,0.20),n=1000) #Note, matrix values as given in problem dscription

#Lets tune this beast 
#Start by plotting initial output to get a visualization
plot(output[,1],output[,2],type="l",main="RWMH 2D",xlab="R",ylab="S",col="goldenrod4")                 

#Now study R and S traceplots more closely
 
plot(output[,1],type="l",main="Traceplot of R",ylab="R",col="olivedrab4")
output[1,1] #The first value is as initialized.
#Then, from this initial value, we see that the 
#values converges after approximatley 110 iterations
#to a R value of about 3.0 < R < 3.07. 

plot(output[,2],type="l",main="Traceplot of S",ylab="S",col="midnightblue")
output[1,1] #The first value is as initialized. 
#Them from this value, we see in the plot that the S
#values converges to a value of approximately 0.17
#after nearly 200 iterations. 

#I therefore re-initialize the function with these new values,
#choosing R of 3.035 and sigma of 0.17
output_2 <- RWMH2D_covid(seir.lp,Sigma=matrix(c(3e-4,0,0,5e-7),nrow=2,ncol=2),thetaInit = c(3.035,0.17),n=1000)

#Study plots once again:
#Overall plot
plot(output_2[,1],output_2[,2],type="l",main="RWMH 2D",xlab="R",ylab="S",col="goldenrod4")  
plot(output_2[,1],type="l",main="Traceplot of R",ylab="R",col="olivedrab4")
plot(output_2[,2],type="l",main="Traceplot of S",ylab="S",col="midnightblue")
#The chosen value for R seems to be spot on. 
#The value chosen for S was too high, as can be seen in the traceplot of S.
#Therefore, I reduce S to 0.167:

output_3 <- RWMH2D_covid(seir.lp,Sigma=matrix(c(3e-4,0,0,5e-7),nrow=2,ncol=2),thetaInit = c(3.035,0.167),n=10000)
plot(output_3[,1],output_3[,2],type="l",main="RWMH 2D",xlab="R",ylab="S",col="goldenrod4")  
plot(output_3[,1],type="l",main="Traceplot of R",ylab="R",col="olivedrab4")
plot(output_3[,2],type="l",main="Traceplot of S",ylab="S",col="midnightblue")


#Check the effective sample sizes: #Note that in output_3, I increased n to 10000 to get enough ESS. 
ESS_R <- coda::effectiveSize(output_3[,1])
ESS_S <- coda::effectiveSize(output_3[,2])
print(paste0("Effective sample size of R is: ",ESS_R))
print(paste0("EFfective sample size of S is: ",ESS_S))

#These values do meet the criteria for effective sample sizes
#of 200 for both S and R as described in the problem formulation.

#Note: I have not removed any burn in in this case. The reason for this is 
#because I have tuned to start within the convergence interval, and this process
#would have been unecessary if I would have moved burn in. 

#6b)
#Compare posterior dist p(theta|y) with prior p(theta)
#Sample from prior dist first
#Using some functions from the solution r-code to problem 4, mandatory 2:

rltrunc_normal <- function(n,l,mu,sigma){
  Fl <- pnorm(q=l,mean=mu,sd=sigma)
  return(qnorm(p=(Fl + runif(n)*(1.0-Fl)),mean=mu,sd=sigma))
}

rRS <- function(n){
  R <- rltrunc_normal(n,1.8,2.2,1.5)
  S <- runif(n,min=0.5/R,0.8/R)
  return(cbind(R,S))
}

#Using this, prior samples S,R can be obtained, as in:
prior_samples <- rRS(10000) #Same amount of values as simulations in a).

plot(prior_samples[,1],prior_samples[,2],xlab="R",ylab="S",col="darkolivegreen3")
points(output_3[,1],output_3[,2],col="orange3")
legend("topright",lty=0,col=c("darkolivegreen3","orange3"),legend = c("Prior","Posterior"),bg="transparent",bty="n",pch=c(1,1))

#From the plot, we see that observing data that the hospital 
#provides us with much information. 


#6c) 
#Redo 4f from mandatory 2.
#Compare plot with and without conditional data.

#Note: Some code is from the solution to mandatory 2, but tweaked to fit.


#Tweak code to choose values from prior samples.
nsim <- 500 #Taken from mandatory 2 solutions
M <- numeric(nsim)
H <- matrix(0.0,nsim,365) # prepare for next point
for(i in 1:nsim){ 
  ff <- seir_sim(R0max = prior_samples[i,1],soc_dist_Rfac = prior_samples[i,2]) #This line of code has been tweaked. Was originally: ff <- seir_sim(R0max = R[i],soc_dist_Rfac = S[i])
  t1 <- which(ff$dates=="2021-01-01")
  t2 <- which(ff$dates=="2021-12-31")
  H[i,] <- rowSums(1000000*ff$ode[t1:t2,c("HH","HC","CC")])
  M[i] <- max(H[i,])
}

# make a data frame to facilitate plotting with dates
df <- data.frame(t(H),dates=ff$dates[t1:t2]) #also from mandatory 2

# in the upper plot, simply show some of the trajectories (also from mandatory 2):
plot(df$dates,df$X1,type="l", 
     xlab="dates in 2021",
     ylab="required hospital beds",
     main="example trajetories, prior") 
for(i in 2:50) lines(df$dates,df[,i],col=i)

# compute daily statistics for plotting, note: also from mandatory 2
daily.stats <- matrix(0.0,365,6)
for(t in 1:365){
  daily.stats[t,1] <- mean(H[,t])
  daily.stats[t,2] <- median(H[,t])
  daily.stats[t,3:6] <- quantile(H[,t],probs=c(0.005,0.1,0.9,0.995))
}

df2 <- data.frame(daily.stats,ff$dates[t1:t2])
colnames(df2) <- c("mean","median","q05","q10","q90","q995","dates")

plot(df2$dates,df2$mean,ylim=c(0,300),col=0,
     xlab="dates in 2021",
     ylab="required hospital beds",
     main="distribution representation, prior") # simply sets up plot window
polygon(c(df2$dates,rev(df2$dates)),c(df2$q05,rev(df2$q995)),
        col="red",density=100)
polygon(c(df2$dates,rev(df2$dates)),c(df2$q10,rev(df2$q90)),
        col="blue",density=100)
lines(df2$dates,df2$mean,lwd=2)
lines(df2$dates,df2$median,col="green")
legend("bottomright",bty="n",bg="transparent",lty=c(1,1,1,1),lwd=c(10,10,2,1),
       legend=c("99%","90%","mean","median"),
       col=c("red","blue","black","green"),)



#Now do all this again for posterior:

#Tweak code to choose values from posterior sample.
nsim <- 500 #Taken from mandatory 2 solutions
M <- numeric(nsim)
H <- matrix(0.0,nsim,365) # prepare for next point
for(i in 1:nsim){ 
  ff <- seir_sim(R0max = output_3[i,1],soc_dist_Rfac = output_3[i,2]) #This line of code has been tweaked. Was originally: ff <- seir_sim(R0max = R[i],soc_dist_Rfac = S[i])
  t1 <- which(ff$dates=="2021-01-01")
  t2 <- which(ff$dates=="2021-12-31")
  H[i,] <- rowSums(1000000*ff$ode[t1:t2,c("HH","HC","CC")])
  M[i] <- max(H[i,])
}

# make a data frame to facilitate plotting with dates
df <- data.frame(t(H),dates=ff$dates[t1:t2]) #also from mandatory 2

# in the upper plot, simply show some of the trajectories (also from mandatory 2):
plot(df$dates,df$X1,type="l", 
     xlab="dates in 2021",
     ylab="required hospital beds",
     main="example trajetories, posterior") 
for(i in 2:50) lines(df$dates,df[,i],col=i)

# compute daily statistics for plotting, note: also from mandatory 2
daily.stats <- matrix(0.0,365,6)
for(t in 1:365){
  daily.stats[t,1] <- mean(H[,t])
  daily.stats[t,2] <- median(H[,t])
  daily.stats[t,3:6] <- quantile(H[,t],probs=c(0.005,0.1,0.9,0.995))
}

df2 <- data.frame(daily.stats,ff$dates[t1:t2])
colnames(df2) <- c("mean","median","q05","q10","q90","q995","dates")

plot(df2$dates,df2$mean,ylim=c(0,300),col=0,
     xlab="dates in 2021",
     ylab="required hospital beds",
     main="distribution representation, posterior") # simply sets up plot window
polygon(c(df2$dates,rev(df2$dates)),c(df2$q05,rev(df2$q995)),
        col="red",density=100)
polygon(c(df2$dates,rev(df2$dates)),c(df2$q10,rev(df2$q90)),
        col="blue",density=100)
lines(df2$dates,df2$mean,lwd=2)
lines(df2$dates,df2$median,col="green")
legend("bottomright",bty="n",bg="transparent",lty=c(1,1,1,1),lwd=c(10,10,2,1),
       legend=c("99%","90%","mean","median"),
       col=c("red","blue","black","green"),)

#Comparing all of the 4 plots produced in this exercise, 
#it becomes clear that the uncertainty is significantly reduced in 
#the posterior case. 

#6d)
#Use additional randomness characterized by (1) in the problem formulation

#Allocate memory for the simulations
n <- 1000
res <- 10^6
des20 <- matrix(0.0,n,31) #31 days in desember
jan21 <- matrix(0.0,n,31) #In january as well
feb21 <- matrix(0.0,n,28) #but not in february

for(i in 1:n) { #Inspired by code lines 505-510 in mandatory2 solution
  ff <- seir_sim(R0max = prior_samples[i,1], soc_dist_Rfac = prior_samples[i,2])
  t1 <- which(ff$dates=="2020-12-01")
  t2 <- which(ff$dates=="2020-12-31")
  des20[i,] <- rowSums(res*ff$ode[t1:t2,c("HH","HC","CC")])
  
  t1 <- which(ff$dates=="2021-01-01")
  t2 <- which(ff$dates=="2021-01-31")
  jan21[i,] <- rowSums(res*ff$ode[t1:t2,c("HH","HC","CC")])
  
  t1 <- which(ff$dates=="2021-02-01")
  t2 <- which(ff$dates=="2021-02-28")
  feb21[i,] <- rowSums(res*ff$ode[t1:t2,c("HH","HC","CC")])
  
}

des20_prior <- rowSums(des20)
jan21_prior <- rowSums(jan21)
feb21_prior <- rowSums(feb21)

data_prior <- data.frame(des20_prior, jan21_prior, feb21_prior)
#dataframe created to make a boxplot. Inspired by source [3]
boxplot(data_prior,ylab="HH+HC+CC",main="Uncertanty of required bed \n using prior data",col="red")

#Now repeat for posterior 
des20 <- matrix(0.0,n,31) #31 days in desember
jan21 <- matrix(0.0,n,31) #In january as well
feb21 <- matrix(0.0,n,28) #but not in february
for(i in 1:n) { #Inspired by code lines 505-510 in mandatory2 solution
  ff <- seir_sim(R0max = output_3[i,1], soc_dist_Rfac = output_3[i,2])
  t1 <- which(ff$dates=="2020-12-01")
  t2 <- which(ff$dates=="2020-12-31")
  des20[i,] <- rowSums(res*ff$ode[t1:t2,c("HH","HC","CC")])
  
  t1 <- which(ff$dates=="2021-01-01")
  t2 <- which(ff$dates=="2021-01-31")
  jan21[i,] <- rowSums(res*ff$ode[t1:t2,c("HH","HC","CC")])
  
  t1 <- which(ff$dates=="2021-02-01")
  t2 <- which(ff$dates=="2021-02-28")
  feb21[i,] <- rowSums(res*ff$ode[t1:t2,c("HH","HC","CC")])
  
}

des20_posterior <- rowSums(des20)
jan21_posterior <- rowSums(jan21)
feb21_posterior <- rowSums(feb21)

data_posterior <- data.frame(des20_posterior,jan21_posterior,feb21_posterior)
boxplot(data_posterior,ylab="HH+HC+CC",main="Uncertanty of required bed \n using posterior data",col="green")

#From the two boxplots, it becomes clear that the need for beds
#can be determined with less uncertainty using the posterior data. 


#For a person who do not have knowlegde about statistical analysis or the SEIR
#model, I would first start by explaining the meaning of S and R. 
#Then I would like to inform that the current boxplots shows the predicted 
#need for hospital beds. I would then stress the meaning of a boxplot: 
# - The fact that the black line in the box is the expected value
# - The fact that worst case scenarios are visalized by the lines
# - The fact that everything in between are other scenarios.
#I would then ask the hospital staff to note the differences betweent the boxplots
#for prior and posterios cases, and further inform that posterior is taking 
#data into account. 
#Lastly, I would stress the fact that the posterior illustrates a much lower difference
#between the expected value (black line) and worst case values, and with this
#let the hospital personell know that they can make a better budget descision 
#using this data rather than the posterior case, where the gap between worst case and 
#mean is bigger. 


#Comment about the boxplot in generally: 
#The box plot graphically depicts the data. 
#The box plots are useful when comparing distributions 
#of various data or groups. 
#The box plot is used to simpler visualize the 
#data to hospital personnel. 
#The boxplot may be compared to a normal distribution
#as in source [4], which also inspired the use of boxplot.
#From the illustration in [4] we clearly see the statistical meaning of
#such plots, and how it describes statistical distributions. 

