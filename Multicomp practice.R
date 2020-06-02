#### Multiple comparisons practice
#### Source: http://sia.webpopix.org/statisticalTests2.html


#####################
# Set directory
#####################

rm(list=ls())

setwd("C:/R")
options(scipen=999)

data <- read.csv("data/ratSurvival.csv")
str(data)
summary(data)

levels(data$regimen)

#####################
# Single comparison between 2 groups
#####################

data.control <- subset(data, regimen=="control" & gender=="female") # Extract the control group of females
data.control$time

# Only 2 rats of this group died before the end of the experiment. On the other hand, 7 females of the group fed with 22% of maize NK693 died during the experiment.
data.test <- subset(data, regimen=="NK603-22%" & gender=="female")
data.test$time

# Because of the (right) censoring process, we cannot just compare the mean survival times using a t-test. On the other hand, we can use the Wilcoxon-Mann-Whitney test which precisely aims to compare the ranks of the survival times in both groups.
wilcox.test(data.test$time, data.control$time, alternative="less") #  tie가 있어 정확한 p값을 계산할 수 없습니다 = 절대값 때문

#####################
# Single comparison...among others
#####################

library(ggplot2) ; theme_set(theme_bw())

pval <- stat <- gender <- regimen <- NULL # Zero vectors

for (g in levels(data$gender)) {
  data.control <- subset(data, regimen=="control" & gender==g)
  for (r in levels(data$regimen)) {
    if (r != "control") {
      data.test <- subset(data, gender==g & regimen==r)
      wt <- wilcox.test(data.test$time, data.control$time, alternative = "less")
      pval <- c(pval, wt$p.value)
      stat <- c(stat, wt$statistic)
      gender <- c(gender, g)
      regimen <- c(regimen, r)
    }
  }
}

R <- data.frame(gender=gender, regimen=regimen, stat=stat, p.value=pval)
R
R <- R[order(R$p.value),]
R

# Plotting...
ggplot(data=R) + geom_point(aes(x=1:18, color=regimen, y=p.value, shape=gender), size=4) +
  scale_y_log10() + xlab("p-value") + scale_x_continuous(breaks=NULL)
# How can we conclue that these diferences are statistically significant?

#####################
# Permutation test
#####################

# imagine that the null hypothesis is true. We can then randomly exchange the labels (i.e. the regimen) and perform the 18 comparisons between the experimental groups and the control groups.
set.seed(100)
dperm.m <- subset(data, gender=="male")
dperm.m
n.m <- dim(dperm.m)[1] # the dimension to be permutated
dperm.m$regimen <- dperm.m$regimen[sample(n.m)] # permute
dperm.f <- subset(data, gender=="female")
n.f <- dim(dperm.f)[1]
dperm.f$regimen <- dperm.f$regimen[sample(n.f)] # permute
dperm <- rbind(dperm.m, dperm.f)

pval <- gender <- regimen <- NULL

for (g in levels(dperm$gender)) {
  dperm.control <- subset(dperm, regimen=="control" & gender==g)
  for (r in levels(dperm$regimen)) {
    if (r != "control") {
      dperm.test <- subset(dperm, gender==g & regimen==r)
      wt <- wilcox.test(dperm.test$time, dperm.control$time, alternative="less")
      pval <- c(pval,wt$p.value)
      gender <- c(gender, g)
      regimen <- c(regimen, r)
    }
  }
}
R.p <- data.frame(gender=gender, regimen=regimen, p.value=pval)
Ro.p <-  R.p[order(R.p$p.value),]

# The test statistics and the p-values now really behave how they are supposed to behave under the null hypothesis
ggplot(data=Ro.p) + geom_point(aes(x=1:18,color=regimen,y=p.value,shape=gender), size=4) + 
  scale_y_log10() + xlab("regimen") + ylab("p-value") + scale_x_continuous(breaks=NULL)

# If we now repeat the same experiment using many different permutations, we will be able to estimate the m distributions of the m test statistics as well as the m distributions of the m p-values under the null hypothesis.
L <- 1000 # the number of repeated experiments
dperm.m <- subset(data,  gender=="male")
n.m <- dim(dperm.m)[1]
dperm.f <- subset(data,  gender=="female")
n.f <- dim(dperm.f)[1]
PV <- ST <- NULL

for (l in (1:L)) {
  dperm.m$regimen <- dperm.m$regimen[sample(n.m)]
  dperm.f$regimen <- dperm.f$regimen[sample(n.f)]
  dperm <- rbind(dperm.m, dperm.f)
  pval <- stat <- NULL
  for (g in levels(dperm$gender)) {
    dperm.control <- subset(dperm, regimen=="control" & gender==g)
    for (r in levels(dperm$regimen)) {
      if (r != "control") {
        dperm.test <- subset(dperm, gender==g & regimen==r)
        wt <- wilcox.test(dperm.test$time, dperm.control$time, alternative="less")
        pval <- c(pval, wt$p.value)
        stat <- c(stat, wt$statistic)
      }
    }
  }
  PV <- rbind(PV, sort(pval))
  ST <- rbind(ST, sort(stat))
}

# We can estimate, for instance, prediction intervals of level 90% for the m=18 ordered p-values
q <- apply(PV, MARGIN = 2, quantile, probs = c(0.05, 0.5, 0.95))
q <- as.data.frame(t(q))
names(q) <- c("low", "median", "up")
q$rank <- 1:dim(q)[1]

# and plot them...
pl <- ggplot(data=q) +
  geom_errorbar(aes(x=rank, ymin=low, ymax=up), width=0.2,size=1.5,colour="grey50") +
  scale_y_log10() + xlab("regimen") + ylab("p-value") + scale_x_continuous(breaks=NULL) 
pl

# ..with original p-values
pl +  geom_point(data=R, aes(x=1:18,color=regimen,y=p.value,shape=gender), size=4) 

# Here, all the p-values, including the smallest ones, belong to the 90% prediction intervals: all the observed p-values behave individually how they are expected to behave under the null hypothesis.
# In particular, when 18 comparisons are performed, it’s not unlikely under the null hypothesis to obtain a smallest p-value less than or equal to the observed one (0.011).
# The probability of such event can easily be estimated by Monte Carlo simulation. Let p(1),ℓ be the smallest p-value obtained from the ℓ-th replicate of the Monte Carlo.
mean(ST[,1]<R$stat[1])

#####################
# Controlling the Family Wise Error Rate
# - The Bonferroni correction
#####################

m <- nrow(R) # the number of null hypothesis
R$pv.bonf <- pmin(1, R$p.value*m)
R
# The Bonferroni correction is appropriate when a single false positive in a set of tests would be a problem. It is mainly useful when there are a fairly small number of multiple comparisons and very few of them might be significant. The main drawback of the Bonferroni correction is its lack of power: it may lead to a very high rate of false negatives.

#####################
# Controlling the False Discovery Rate
# - Detecting assoications
#####################

data <- read.csv("data/dietary.csv")
data
# We can see that five of the variables show a significant p-value (<0.05). However, because García-Arenzana et al. (2014) tested 25 dietary variables, we would expect one or two variables to show a significant result purely by chance, even if diet had no real effect on mammographic density.
# Applying the Bonferroni correction, we divide α=0.05 by the number of tests (m=25) to get the Bonferroni critical value, so a test would have to have p<0.002 to be significant. Under that criterion, only the test for total calories is significant.
# An alternative approach is to control the false discovery rate, i.e the expected proportion of ``discoveries" (significant results) that are actually false positives. FDR control offers a way to increase power while maintaining some principled bound on error.

#####################
# Controlling the False Discovery Rate
# - The Benjamini-Hochberg procedure
#####################

p.bh <- data$p.value
m <- length(p.bh)

for (i in ((m-1):1)) 
  p.bh[i] <- min(data$p.value[i]*m/i , p.bh[i+1])

data$p.bh <- p.bh

# the largest p-value that has P(i)<(i/m)α is significant, and all of the P-values smaller than it are also significant, even the ones that aren’t less than their Benjamini-Hochberg critical value α×i/m

alpha <- 0.25
m <- dim(data)[1]
data$critical.value <- (1:m)/m*alpha
data

library(gridExtra)
pl1 <- ggplot(data) + geom_line(aes(x=1:m,y=p.value), colour="blue") +
  geom_line(aes(x=1:m,y=critical.value), colour="red") +xlab("(i)")
grid.arrange(pl1,pl1 + xlim(c(1,8)) + ylim(c(0,0.21)) + geom_vline(xintercept=5.5))

# The FDR is not bounded by α, but by (m0⋅/m)α. We could increase the global power of the tests and get a FDR equal to the desired level α, either by defining the critical values as (i/m0⋅)α, or by multiplying the adjusted p-values by m0⋅/m.
# Unfortunately, m0⋅ is unknown… but it can be estimated, as the number of non significant tests for instance.

m0.est <- sum(data$p.bh>alpha)
data$crit.valc <- round(data$critical.value*m/m0.est,4)
data$p.bhc <- round(data$p.bh*m0.est/m,4)
head(data,10)

#####################
# Controlling the False Discovery Rate
# - The Benjamini-Hochberg procedure
#####################

nx <- 50
ny <- 50
m0 <- 120
m1 <- 20
mu.x <- 0
mu.y <- c(rep(0, m0), seq(0.3, 0.6, length=m1))

# For each of the L=1000 simulated replicate of the same experiment, we will randomly sample observation x and y from the model and perform a t-test for each of the m=140 variables. We therefore get m=140 p-values for each of these L=1000 replicates.
L <- 1000
m <- m0+m1

set.seed(1234)

pval <- matrix(ncol=L, nrow=m)

for (l in (1:L)) {
  x.sim <- matrix(rnorm(nx*m, mu.x), ncol=nx)
  y.sim <- matrix(rnorm(ny*m, mu.y), ncol=ny)
  dat <- cbind(x.sim, y.sim)
  pval[,l] <- apply(dat, 1, function(dat) {
    t.test(x = dat[1:nx], y = dat[(nx +  1):(nx + ny)])$p.value})
}


alpha=0.2
m01 <- colSums(pval[1:m0,] < alpha)
mean(m01/m0) # the proportion of wrongly rejected null hypotheses

m11 <- colSums(pval[(m0+1):m,] < alpha)
mean(m11/m1) # The proportion of correctly rejected null hypotheses is an estimate of the power of the test

mean(m01/(m01+m11)) #  about 60% of them are false discoveries.

# Let us now apply the Bonferroni correction, and compute the proportion of wrongly and correctly rejected null hypotheses and the proportion of false discoveries
pval.b <- apply(pval,2, function(pval) {p.adjust(pval, method="bonferroni")})
m01.b <- colSums(pval.b[1:m0,] < alpha)
m11.b <- colSums(pval.b[(m0+1):m,] < alpha)
c(mean(m01.b/m0), mean(m11.b/m1), mean(m01.b/(m01.b+m11.b), na.rm=TRUE))

# Very few of the true null hypotheses are rejected, which may be a good point, but the price to pay is a very low power (less that 20%).
# The familywise error rate (FWER) is the probability P(m01≥1) to reject at least one of the true null hypothesis. This probability remains quite small when the Bonferroni correction is used.
mean(m01.b>=1)

# The Benjamini-Hochberg correction increases the power and controls the FDR as expected since the proportion of false discoveries remains below the level α.
pval.bh <- apply(pval,2, function(pval) {p.adjust(pval,method="BH")})
m01.bh <- colSums(pval.bh[1:m0,] < alpha)
m11.bh <- colSums(pval.bh[(m0+1):m,] < alpha)
c(mean(m01.bh/m0), mean(m11.bh/m1),mean(m01.bh/(m01.bh+m11.bh),na.rm=TRUE))

# Lastly, we can slightly improve the BH procedure by multiplying the p-values by the ratio m^0⋅/m
pval.bhc <- pval.bh

for (l in (1:L)) {
  m0e <- sum(pval.bh[,l]>alpha)
  pval.bhc[,l] <- pval.bh[,l]*m0e/m
}

m01.bhc <- colSums(pval.bhc[1:m0,] < alpha)
m11.bhc <- colSums(pval.bhc[(m0+1):m,] < alpha)
c(mean(m01.bhc/m0), mean(m11.bhc/m1),mean(m01.bhc/(m01.bhc+m11.bhc),na.rm=TRUE))

#############################
###### addressing multiple comparisons
###### Source: https://sas-and-r.blogspot.com/2012/05/example-930-addressing-multiple.html
#############################

pvals = c(.001, .001, .001, .02, .22, .59, .87)
BONF = p.adjust(pvals, "bonferroni")
BH = p.adjust(pvals, "BH")
res = cbind(pvals, BH=round(BH, 3), BONF=round(BONF, 3))

matplot(res, ylab="p-values", xlab="sorted outcomes")
abline(h=0.05, lty=2)
matlines(res)
legend(1, .9, legend=c("Bonferroni", "Benjamini-Hochberg", "Unadjusted"), 
       col=c(3, 2, 1), lty=c(3, 2, 1), cex=0.7)