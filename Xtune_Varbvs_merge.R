library(varbvs)
library(xtune)
library(pROC)

#Simulating genetics data as in Carbonetto, P. and Stephens, M. (2012). 

p <- 10000
n <- 500
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) + (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)

# Generate the ground-truth regression coefficients for the variables
# (X) and additional 2 covariates (Z).
u    <- c(-1,1)
beta <- as.matrix(c(0.5*rnorm(150),rep(0,9850)),ncol=1)
Zgr <- c(sample(c(rep(0,50),rep(1,100))),sample(c(rep(0,9450),rep(1,400))))
pbeta <- c(rbeta(75,0.1,10),rbeta(75,1,5),runif(9850))
resp <- as.numeric(beta != 0)
y <- as.double(X %*% beta + rnorm(n))

# Fit the basic spike-and-slab model using varbvs.
q <- 1/100
fit0 <- varbvs(X=X,Z=NULL,y=y,family="gaussian",logodds = log10(q/(1-q)))
print(summary(fit0))

#Apply xtune to estimate co-data adapted ridge penalties
#First, both co-data sources; may take a few minutes (have a coffee)
mm <- cbind(log(pbeta),Zgr)
system.time(res.xtune <- xtune(X,y,Z=mm,family="linear",c=0))

#Repeat, but use only pbeta as co-data
mm1 <- cbind(log(pbeta)) 
system.time(res.xtune1 <- xtune(X,y,Z=mm1,family="linear",c=0))

#Repeat, but use only Zgr as co-data
mm2 <- cbind(Zgr)
system.time(res.xtune2 <- xtune(X,y,Z=mm2,family="linear",c=0))

#Estimate feature specific inclusion probabilities for the spike-and-slab
#from the penalties. Normalize w.r.t. overall q. 
  q1 <- 1/res.xtune$penalty.vector
  q1m <- median(q1)
  q1[q1>=100*q1m] <- 100*q1m
  q2 <- q*q1/mean(q1)
  
#summaries for the relevant and irrelevant features
summary(q2[1:150])  
summary(q2[-(1:150)])  

#fit spike-and-slab with feature-specific prior inclusion probabilities
lodds <- matrix(log10(q2*(1-q2)),ncol=1)
fit <- varbvs(X=X,Z=NULL,y=y,family="gaussian",logodds = lodds)

#Repeat, but use only pbeta as co-data
q1 <- 1/res.xtune1$penalty.vector
q1m <- median(q1)
q1[q1>=100*q1m] <- 100*q1m
q2 <- q*q1/mean(q1)
lodds <- matrix(log10(q2*(1-q2)),ncol=1)
fit1 <- varbvs(X=X,Z=NULL,y=y,family="gaussian",logodds = lodds)

#Repeat, but use only Zgr as co-data
q1 <- 1/res.xtune2$penalty.vector
q1m <- median(q1)
q1[q1>=100*q1m] <- 100*q1m
q2 <- q*q1/mean(q1)
lodds <- matrix(log10(q2*(1-q2)),ncol=1)
fit2 <- varbvs(X=X,Z=NULL,y=y,family="gaussian",logodds = lodds)

#Save results
save(res.xtune,res.xtune1,res.xtune2, fit0,fit1,fit2,fit, file="resXtunevarbvs.Rdata")

#Plot roc-curves for variable selection using all four spike-and-slab fits 
par(mfrow=c(2,2))
r0 <- roc(resp,fit0$pip)
r1 <- roc(resp,fit1$pip)
r2 <- roc(resp,fit2$pip)
r <- roc(resp,fit$pip)
par(mfrow=c(1,4), mar=c(4,2,2,1))
plot.roc(r0,print.auc=TRUE, ylim=c(0,1),main="No Co-data",auc.polygon=TRUE,auc.polygon.col="lightblue")
plot.roc(r1,print.auc=TRUE, main = "log-p only",auc.polygon=TRUE,auc.polygon.col="lightblue")
plot.roc(r2,print.auc=TRUE, main = "Groups only",auc.polygon=TRUE,auc.polygon.col="lightblue")
plot.roc(r,print.auc=TRUE, main = "Both",auc.polygon=TRUE,auc.polygon.col="lightblue")


