#Demo Squeezy with hyperparameter shrinkage
#Squeezy: written by Mirrelijn van Nee; squeezy_hyper: Mark van de Wiel
library(squeezy)
library(multiridge)
library(glmnet)

#Set working directory containing the squeezy_hyper.R script
setwd("C:/Synchr/Rscripts/CodataRegrReview/")
source("squeezy_hyper.R")

#####################
# Simulate toy data #
#####################

p<-200 #number of covariates
n<-100 #sample size training data set
n2<-100 #sample size test data set
G<- 5 #number of groups

taugrp <- rep(c(0.05,0.1,0.2,0.5,1),each=p/G) #ridge prior variance
groupIndex <- rep(1:G,each=p/G) #groups for co-data
groupset <- lapply(1:G,function(x){which(groupIndex==x)}) #group set with each element one group
sigmasq <- 2 #linear regression noise
lambda1 <- sqrt(taugrp/2) #corresponding lasso penalty
#A Laplace(0,b) variate can also be generated as the difference of two i.i.d.
#Exponential(1/b) random variables
betas <-   rexp(p, 1/lambda1) -  rexp(p, 1/lambda1) #regression coefficients
X <- matrix(rnorm(n*p),n,p) #simulate training data
Y <- rnorm(n,X%*%betas,sd=sqrt(sigmasq))

######################
# Fit squeezy, lasso #
######################

res.squeezy <- squeezy(Y,X,groupset=groupset, model="linear",alpha=1)

######################
# Fit squeezy, lasso #
# hyper              #
######################

#fit squeezy with one feature group to obtain target lambda
groupset1 <- list(1:p)
res.squeezy1 <- squeezy(Y,X,groupset=groupset1, alpha=1,compareMR=FALSE)
lambdas1 <- res.squeezy1$lambdaMR
print("Target lambda:")
print(lambdas1)
res.squeezyshr <- squeezy_hyper(Y,X,groupset=groupset, alpha=1,compareMR=FALSE, penhyper=TRUE, targlam =lambdas1)
print("Shrunken lambdas:")
print(res.squeezyshr$lambdaMR)

#Compare with unshrunken ones:
print(res.squeezy$lambdaMR)
