library(ecpc)
#set the directory where you stored the data
#wd <- "C:/Synchr/Rscripts/CodataRegrReview/"
setwd(wd)
show(load("XYcodata.Rdata"))
load("XYcodata.Rdata")

abund <- codata$abund; minlogFDR1 <- codata$minlogFDR1; minlogFDR2 <- codata$minlogFDR2

#proprocessing: create spline basis for co-data
Z1 <- createZforSplines(values=abund)
S1.Z1 <- createS(orderPen=2, G=dim(Z1)[2]) #create difference penalty matrix
Con1 <- createCon(G=dim(Z1)[2], shape="positive+monotone.i") #create constraints

Z2 <- createZforSplines(values=minlogFDR1)
S1.Z2 <- createS(orderPen=2, G=dim(Z2)[2]) #create difference penalty matrix
Con2 <- createCon(G=dim(Z2)[2], shape="positive+monotone.i") #create constraints

Z3 <- createZforSplines(values=minlogFDR2)
S1.Z3 <- createS(orderPen=2, G=dim(Z3)[2]) #create difference penalty matrix
Con3 <- createCon(G=dim(Z3)[2], shape="positive+monotone.i") #create constraints

Z.all <- list(Z1=Z1,Z2=Z2,Z3=Z3)
paraPen.all <- list(Z1=list(S1=S1.Z1),Z2=list(S1=S1.Z2), Z3=list(S1=S1.Z3))
paraCon <- list(Z1=Con1,Z2=Con2, Z3=Con3)

#fit ecpc: ridge regression with co-data adaptive penalties and posthoc selection
fit <- ecpc(Y,X,Z = Z.all, paraPen = paraPen.all, paraCon = paraCon, model="logistic", maxsel=50)

wh <- which(fit$betaPost != 0)

#plotting
par(mfrow=c(2,3),mar=c(4,4,3,1), cex=1.2)
i <-1
groupsetNO <- c(unlist(sapply(1:length(Z.all),function(i) rep(i,dim(Z.all[[i]])[2]))))
rvk1 <- as.vector(Z.all[[i]]%*%fit$gamma[groupsetNO==i])*fit$tauglobal
plot(abund,rvk1, ylab="Inverse penalty factor", xlab="abundance",xlim=c(-0.2,0.2), type="p")
par(mar=c(4,3,3,1))
i <-2
groupsetNO <- c(unlist(sapply(1:length(Z.all),function(i) rep(i,dim(Z.all[[i]])[2]))))
rvk2 <- as.vector(Z.all[[i]]%*%fit$gamma[groupsetNO==i])*fit$tauglobal
plot(minlogFDR1,rvk2, xlab="- log FDR1",xlim=c(0,10))

i <-3
groupsetNO <- c(unlist(sapply(1:length(Z.all),function(i) rep(i,dim(Z.all[[i]])[2]))))
rvk3 <- as.vector(Z.all[[i]]%*%fit$gamma[groupsetNO==i])*fit$tauglobal
plot(minlogFDR2,rvk3, xlab="- log FDR2", xlim=c(0,10))

par(mar=c(2,4,1,1))
d1 <- density(abund)
d2 <- density(abund[wh])
ylim <- c(0, max(d1$y,d2$y))
plot(density(abund),lwd=2,ylim=ylim, main="")
points(d2, type="l", col="grey",lwd=2)

par(mar=c(2,3,1,1))
d1 <- density(minlogFDR1)
d2 <- density(minlogFDR1[wh])
ylim <- c(0, max(d1$y,d2$y))
plot(d1,lwd=2,ylim=ylim, main="")
points(d2, type="l", col="grey",lwd=2)
legend(2.2,0.38, legend=c("All features", "Selected features"), lwd=c(2,2),
       col=c("black","grey"))
d1 <- density(minlogFDR2)
d2 <- density(minlogFDR2[wh])
ylim <- c(0, max(d1$y,d2$y))
plot(d1,lwd=2,ylim=ylim, main="")
points(d2, type="l", col="grey",lwd=2)
#dev.off()

#evaluate predictive performance and compare with basic ridge model; may take a while

#library that includes function to create balanced folds (alternative: caret)
library("multiridge") 

folds <- CVfolds(Y,model="logistic", kfold=10,nrepeat=3)
MSEs <- c()
for(k in 1:30){
  #k <-1
  outs <- folds[[k]]
  ins <- (1:n)[-outs]
  fit <- ecpc(Y=Y[ins],X=X[ins,],Z = Z.all, paraPen = paraPen.all, paraCon = paraCon,
              model="logistic", Y2 = Y[outs], X2=X[outs,],postselection=FALSE)
  MSEs <- rbind(MSEs,c(fit$MSEecpc,fit$MSEridge))
  print(MSEs)
}
colMeans(MSEs)

