#### Power anaylsis simulation practice
#### Source: http://egap.org/methods-guides/10-things-you-need-know-about-statistical-power

## REMOVE ALL ##
rm(list=ls())

###############################
###### Calculating power ######
###############################

mu_t # the average outcome in treatment group
mu_c # the average outcome in control group
sigma # standard deviation of outcomes (σ)
N # total number of subjects

power_calculator <- function(mu_t, mu_c, sigma, alpha=0.05, N){
  lowertail <- (abs(mu_t - mu_c)*sqrt(N))/(2*sigma)
  uppertail <- -1*lowertail
  beta <- pnorm(lowertail-qnorm(1-alpha/2), lower.tail=TRUE) + 1 - pnorm(uppertail-qnorm(1-alpha/2), lower.tail=FALSE)
  return(beta)
}

power_calculator(mu_t = 65, mu_c = 60, sigma = 20, N = 500)
# OR
power_calculator(mu_t = 65, mu_c = 60, sigma = 20, alpha = 0.05, N = 500) # You can change the alpha level if you want

##########################################
###### Simulations to estimate power######
##########################################

####### The stanard desgin : randomly assigns subjects to either treatment or control with probability 0.5

possible.ns <- seq(from = 100, to = 2000, by=50) # The sample sizes we'll be considering
powers <- rep(NA, length(possible.ns)) # Empty object to collect simulation estimates
alpha <- 0.05
sims <- 500 # Number of simulations to conduct for each N

#### Outer loop to vary the number of subjects #### 
for (j in 1:length(possible.ns)){ N <- possible.ns[j] # Pick the jth value for N 

Y0 <- rnorm(n=N, mean=60, sd=20) # control potential outcome 
tau <- 5 # Hypothesize treatment effect 
Y1 <- Y0 + tau # treatment potential outcome
significant.experiments <- rep(NA, sims) # Empty object to count significant experiments

#### Inner loop to conduct experiments "sims" times over for each N #### 
for (i in 1:sims){
  Z.sim <- rbinom(n=N, size=1, prob=.5) # Do a random assignment 
  Y.sim <- Y1*Z.sim + Y0*(1-Z.sim)  # Reveal outcomes according to assignment 
  fit.sim <- lm(Y.sim ~ Z.sim) # Do analysis (Simple regression)
  p.value <- summary(fit.sim)$coefficients[2,4] # Extract p-values
  significant.experiments[i] <- (p.value <= alpha) # Determine significance according to p <= 0.05
}
powers[j] <- mean(significant.experiments) # store average success rate (power) for each N 
}
powers

plot(possible.ns, powers, ylim=c(0,1))
abline(h=0.8, col="red")

###### For covariate control

rm(list=ls())
possible.ns <- seq(from=100, to=2000, by=50)
powers <- rep(NA, length(possible.ns))
powers.cov <- rep(NA, length(possible.ns)) # Need a second empty vector
alpha <- 0.05
sims <- 500

for (j in 1:length(possible.ns)){
  N <- possible.ns[j]
  
  significant.experiments <- rep(NA, sims)
  significant.experiments.cov <- rep(NA, sims) # Need a second empty vector, too
  
  for (i in 1:sims){
    gender <- c(rep("F", N/2), rep("M", N/2)) # Generate "gender" covariate
    age <- sample(x=18:65, size=N, replace=TRUE) # Generate "age" covariate
    EffectOfGender <- 10 # Hypothesize the "effect" of gender on income
    EffectOfAge <- 2 # Hypothesize the "effect" of age on income
    
    ## Hypothesize Control Outcome as a function of gender, age, and error
    Y0 <- EffectOfGender*(gender=="M") + EffectOfAge*age + rnorm(n=N, mean=100, sd=20)
    
    ## This is all the smae ##
    tau <- 5
    Y1 <- Y0 + tau
    Z.sim <- rbinom(n=N, size=1, prob=.5)
    Y.sim <- Y1*Z.sim + Y0*(1-Z.sim)
    fit.sim <- lm(Y.sim ~ Z.sim)
    
    ## This is the novel analysis -- including two covariates to increase precision ##
    fit.sim.cov <- lm(Y.sim ~ Z.sim + (gender=="M") + age)
    
    ## extract p-values and calculate signifiance ##
    p.value <- summary(fit.sim)$coefficients[2,4]
    p.value.cov <- summary(fit.sim.cov)$coefficients[2,4]
    significant.experiments[i] <- (p.value <= alpha)
    significant.experiments.cov[i] <- (p.value.cov <= alpha)
  }
  
  powers[j] <- mean(significant.experiments)
  powers.cov[j] <- mean(significant.experiments.cov)
}

plot(possible.ns, powers, ylim=c(0,1))
points(possible.ns, powers.cov, col="blue")
abline(h=0.8, col="red")

###### For multiple treatments

rm(list=ls())
#install.packages("randomizr")
library(randomizr) # ramndomizr package for complete random assignment

possible.ns <- seq(from = 100, to = 5000, by = 100)
power.atleastone <- rep(NA, length(possible.ns)) # The probability that at least one of the treatments turns up significant.
power.bothtreatments <- rep(NA, length(possible.ns)) # The probability that all of the treatments turn up significant.
power.fullranking <- rep(NA, length(possible.ns))# The probability that the treatment effects will be in the hypothesized ranking, and that all the differences are significant.
alpha <- 0.1 # (one-tailed test at. .05 level)
sims <- 100

#### Outer loop to vary the number of subjects ####
for (j in 1:length(possible.ns)){
  N <- possible.ns[j]
  p.T1vsC <- rep(NA, sims)
  p.T2vsC <- rep(NA, sims)
  p.T2vsT1 <- rep(NA, sims)
  c.T1vsC <- rep(NA, sims)
  c.T2vsC <- rep(NA, sims)
  c.T2vsT1 <- rep(NA, sims)
  #### Inner loop to conduct experiments "sims" times over for each N ####
  for (i in 1:sims){
    Y0 <- rnorm(n=N, mean=60, sd=20)
    tau_1 <- 2.5
    tau_2 <- 5
    Y1 <- Y0 + tau_1
    Y2 <- Y0 + tau_2
    Z.sim <- complete_ra(N=N, num_arms=3) #num_arms = the number of treatment arms
    Y.sim <- Y0*(Z.sim=="T3") + Y1*(Z.sim=="T1") + Y2*(Z.sim=="T2")
    frame.sim <- data.frame(Y.sim, Z.sim)
    fit.T1vsC.sim <- lm(Y.sim ~ Z.sim=="T1", data=subset(frame.sim, Z.sim!="T2"))
    fit.T2vsC.sim <- lm(Y.sim ~ Z.sim=="T2", data=subset(frame.sim, Z.sim!="T1"))
    fit.T2vsT1.sim <- lm(Y.sim ~ Z.sim=="T2", data=subset(frame.sim, Z.sim!="T3"))
    
    ### Need to capture coefficients and pvalues (one-tailed tests, so signs are important)
    c.T1vsC[i] <- summary(fit.T1vsC.sim)$coefficients[2,1]
    c.T2vsC[i] <- summary(fit.T2vsC.sim)$coefficients[2,1]
    c.T2vsT1[i] <- summary(fit.T2vsT1.sim)$coefficients[2,1]
    p.T1vsC[i] <- summary(fit.T1vsC.sim)$coefficients[2,4]
    p.T2vsC[i] <- summary(fit.T2vsC.sim)$coefficients[2,4]
    p.T2vsT1[i] <- summary(fit.T2vsT1.sim)$coefficients[2,4]
  }
  power.atleastone[j] <- mean(c.T1vsC>0 & c.T2vsC>0 & (p.T1vsC < alpha/2 | p.T2vsC < alpha/2))
  power.bothtreatments[j] <- mean(c.T1vsC>0 & c.T2vsC>0 & p.T1vsC < alpha/2 & p.T2vsC < alpha/2)
  power.fullranking[j] <- mean(c.T1vsC>0 & c.T2vsC>0 & c.T2vsT1 > 0 & p.T1vsC < alpha/2 & p.T2vsT1 < alpha/2)
  print(j)
}

plot(possible.ns, power.atleastone, ylim=c(0,1))
points(possible.ns, power.bothtreatments, col="purple")
points(possible.ns, power.fullranking, col="blue")
abline(h=0.8, col="red")