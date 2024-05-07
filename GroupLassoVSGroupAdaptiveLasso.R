#Comparison of sparse group lasso and group-adaptive lasso in terms of 
#feature selection.NOTE: Results my deviate somewhat from those in the manuscript, 
#as the training data and the true betas are regenerated. 
#They should, however, qualitatively agree.
#Code: Mark van de Wiel

library(squeezy)
library(SGL)
library(glmnet)

####################### 
# Simulation 1
# 1/3rd of groups
# with signal
#######################

n <- 200         ## Nb of observations  
p <- 2000
Nblock <- 1   
CorX <- 0.0     ## Correlation level between variables

#generate equicorrelated data matrix X
pblock <- p/Nblock
X <- Reduce(cbind,lapply(1:Nblock, function(z) matrix(rep(rnorm(n,sd=sqrt(CorX/(1-CorX))),times=pblock),n,pblock))) + 
  matrix(rnorm(n*p),n,p)
X <- t((t(X) - apply(t(X),1,mean))/apply(t(X),1,sd))


# Number of feature groups
Gs <- c(3,6,9,15,24,39,60,99)

#Number of repeats
nrep <- 25
allres <- list()
for(G in Gs){
  print(G)
  resG <- list()
  for(k in 1:nrep){
    print(k)
    Gnon0 <- G/3
    betaprop <- rbeta(Gnon0,2,6)
    Gsize <- ceiling(p/G)
    nbeta <- ceiling(Gsize*betaprop)
    nbeta
    
    #generate betas
    betainit <- sapply(1:Gnon0,function(i) {
      ni <- nbeta[i];
      betas <- c(rt(ni,3),rep(0,Gsize-ni)) 
      return(betas)
      })
    betainitall <- c(unlist(betainit),rep(0,p-Gsize*Gnon0))
    betasc <- betainitall/(sqrt(p)*sd(betainitall))
    nnonzero <- sum(betasc!=0)
    nnonzero
    
    #generate response
    Y <- as.numeric(X %*% betasc + rnorm(n))
    G0 <- G*2/3
    Gsize0 <- floor(p/G) 
    Gsizelast <- p-Gnon0*Gsize - (G0-1)*Gsize0
    
    #group indices (required by sgl)
    groupind <- c(unlist(lapply(1:Gnon0,function(i) rep(i,Gsize))),
                          unlist(lapply(Gnon0+(1:(G0-1)),function(i) rep(i,Gsize0))),
                          rep(G,Gsizelast))
    
    #convert to group sets (required by squeezy)
    groupset <- lapply(1:G,function(x){which(groupind==x)})
    
    #apply squeezy
    res.squeezy <- squeezy(Y,X,groupset=groupset, alpha=1,compareMR=FALSE)
    datXY <- list(x=X,y=Y)
    
    #apply sparse group lasso
    res.sgl <- SGL(data=datXY,index=groupind,nlam=100,standardize=FALSE)
    res <- list(G=G,k=k,nnonzero=nnonzero,betatrue=betasc,res.sq = res.squeezy,res.sgl=res.sgl)
    resG <- c(resG,list(res))
  }
  allres <- c(allres,list(resG))
  save(allres,file="allres.Rdata")
}



#### Analyse results: Compare feature selection 
load("allres.Rdata")
Gs <- c(3,6,9,15,24,39,60,99)
nrep <- 25
ntosels <- c(25,50)

par(mfrow=c(2,1),mar=c(3,2,1,1))
for(ntosel in ntosels){
allgF1sq <- allgF1sgl <- c()
allgR2sq <- allgR2sgl <- c()
for(g in 1:length(Gs)){
  print(Gs[g])
  allF1sq <- allF1sgl <- c()
  allR2sq <- allR2sgl <- c()
  for(k in 1:nrep){
    res.squeezy <- allres[[g]][[k]]$res.sq
    res.sgl <-  allres[[g]][[k]]$res.sgl
    
    #true betas
    betasc <- allres[[g]][[k]]$betatrue
    p <- length(betasc)
     
    #selection by squeezy
    sqbetas <- res.squeezy$glmnet.fit$beta
    nselsq <- apply(sqbetas,2,function(vec) sum(vec!=0))
    nselsq
    whsq <- which.max(nselsq-ntosel > 0) - 1
    betasqsel <- sqbetas[,whsq]
    
    #selection by sgl
    betassgl <- res.sgl$beta
    nsel <- apply(betassgl,2,function(vec) sum(vec!=0))
    whsgl <- which.max(nsel-ntosel > 0) - 1
    betasglsel <- betassgl[,whsq]
     
    betarandnon0 <- sample(1:p,ntosel)
    betarandsel <- rep(0,p)
    betarandsel[betarandnon0] <- 1
    
    #indices true non-zero betas
    whpos <- which(betasc!=0)
    
    #compute F1 score squeezy
    TPsq <- sum(betasqsel[whpos]!=0)
    FPsq <- sum(betasqsel[-whpos]!=0)
    FNsq <- sum(betasqsel[whpos]==0)
    recallsq <- TPsq/(FNsq + TPsq)
    precsq <- TPsq/(FPsq + TPsq )
    F1sq <- 2*precsq*recallsq/(precsq + recallsq)
    
    #compute F1 score SGL
    TPsgl <- sum(betasglsel[whpos]!=0)
    FPsgl <- sum(betasglsel[-whpos]!=0)
    FNsgl <- sum(betasglsel[whpos]==0)
    recallsgl <- TPsgl/(FNsgl + TPsgl)
    precsgl <- TPsgl/(FPsgl + TPsgl )
    F1sgl <- 2*precsgl*recallsgl/(precsgl + recallsgl)

    
    allF1sq <- c(allF1sq,F1sq)
    allF1sgl <- c(allF1sgl,F1sgl)
  
  }
  allgF1sq <- cbind(allgF1sq,allF1sq)
  allgF1sgl <- cbind(allgF1sgl,allF1sgl)
}

#plotting
colnames(allgF1sq) <- colnames(allgF1sgl) <- paste("G",Gs,sep="")
ng <- length(Gs)
allgF1s <- matrix(data=NA,nrow=nrep,ncol=2*ng)
allgF1s[,seq(1,2*ng,by=2)] <- allgF1sq
allgF1s[,seq(2,2*ng,by=2)] <- allgF1sgl
colnames(allgF1s) <- paste("G",sort(rep(Gs,2)),sep="")
maint <- paste("F1 score, p_sel=", ntosel,sep="")
boxplot(allgF1s,col=rep(c("green","blue"),8), main=maint)
abline(v=seq(2,2*(ng-1),by=2)+0.5,lty=3)
}


####################### 
# Simulation 2: 
# group-sparse setting
#######################

#As above, but now the group sparse setting. Larger p, five non-zero feature groups.
n <- 200         ## Nb of observations  
ntest <- 1000
nboth <- n+ntest
p <- 10000
Nblock <- 1   
CorX <- 0.0     ## Correlation level between variables 

#generate equicorrelated data matrix X
pblock <- p/Nblock
X <- Reduce(cbind,lapply(1:Nblock, function(z) matrix(rep(rnorm(n,sd=sqrt(CorX/(1-CorX))),times=pblock),n,pblock))) + 
  matrix(rnorm(n*p),n,p)
X <- t((t(X) - apply(t(X),1,mean))/apply(t(X),1,sd))

Gs <- c(60,99)
nrep <- 25
allres <- list()
for(G in Gs){
  print(G)
  resG <- list()
  for(k in 1:nrep){
    print(k)
    Gnon0 <- 5
    betaprop <- rbeta(Gnon0,10,90)
    Gsize <- ceiling(p/G)
    nbeta <- ceiling(Gsize*betaprop)
    nbeta
    betainit <- sapply(1:Gnon0,function(i) {
      ni <- nbeta[i];
      betas <- c(rt(ni,3),rep(0,Gsize-ni)) 
      return(betas)
    })
    betainitall <- c(unlist(betainit),rep(0,p-Gsize*Gnon0))
    betasc <- betainitall/(sqrt(p)*sd(betainitall))
    nnonzero <- sum(betasc!=0)
    nnonzero
    Y <- as.numeric(X %*% betasc + rnorm(n))
    G0 <- G - Gnon0
    Gsize0 <- floor(p/G) 
    Gsizelast <- p-Gnon0*Gsize - (G0-1)*Gsize0
    groupind <- c(unlist(lapply(1:Gnon0,function(i) rep(i,Gsize))),
                  unlist(lapply(Gnon0+(1:(G0-1)),function(i) rep(i,Gsize0))),
                  rep(G,Gsizelast))
    groupset <- lapply(1:G,function(x){which(groupind==x)}) #group set with each element one group
    res.squeezy <- squeezy(Y,X,groupset=groupset, alpha=1,compareMR=FALSE)
    datXY <- list(x=X,y=Y)
    res.sgl <- SGL(data=datXY,index=groupind,nlam=100,standardize=FALSE)
    res <- list(G=G,k=k,nnonzero=nnonzero,betatrue=betasc,res.sq = res.squeezy,res.sgl=res.sgl)
    resG <- c(resG,list(res))
  }
  allres <- c(allres,list(resG))
  save(allres,file="allresmoregroupsparse.Rdata")
}


load("allresmoregroupsparse.Rdata")
Gs <- c(60,99)
nrep <- 25
ntosels <- c(25,50)
par(mfrow=c(2,1),mar=c(3,2,1,1))
for(ntosel in ntosels){
  allgF1sq <- allgF1sgl <- c()
  allgR2sq <- allgR2sgl <- c()
  for(g in 1:length(Gs)){
    print(Gs[g])
    allF1sq <- allF1sgl <- c()
    allR2sq <- allR2sgl <- c()
    for(k in 1:nrep){
      res.squeezy <- allres[[g]][[k]]$res.sq
      res.sgl <-  allres[[g]][[k]]$res.sgl
      
      #true betas
      betasc <- allres[[g]][[k]]$betatrue
      p <- length(betasc)
      
      #selection by squeezy
      sqbetas <- res.squeezy$glmnet.fit$beta
      nselsq <- apply(sqbetas,2,function(vec) sum(vec!=0))
      nselsq
      whsq <- which.max(nselsq-ntosel > 0) - 1
      betasqsel <- sqbetas[,whsq]
      
      #selection by sgl
      betassgl <- res.sgl$beta
      nsel <- apply(betassgl,2,function(vec) sum(vec!=0))
      whsgl <- which.max(nsel-ntosel > 0) - 1
      betasglsel <- betassgl[,whsq]
      
      betarandnon0 <- sample(1:p,ntosel)
      betarandsel <- rep(0,p)
      betarandsel[betarandnon0] <- 1
      
      #indices true non-zero betas
      whpos <- which(betasc!=0)
      
      #compute F1 score squeezy
      TPsq <- sum(betasqsel[whpos]!=0)
      FPsq <- sum(betasqsel[-whpos]!=0)
      FNsq <- sum(betasqsel[whpos]==0)
      recallsq <- TPsq/(FNsq + TPsq)
      precsq <- TPsq/(FPsq + TPsq )
      F1sq <- 2*precsq*recallsq/(precsq + recallsq)
      
      #compute F1 score SGL
      TPsgl <- sum(betasglsel[whpos]!=0)
      FPsgl <- sum(betasglsel[-whpos]!=0)
      FNsgl <- sum(betasglsel[whpos]==0)
      recallsgl <- TPsgl/(FNsgl + TPsgl)
      precsgl <- TPsgl/(FPsgl + TPsgl )
      F1sgl <- 2*precsgl*recallsgl/(precsgl + recallsgl)
      
      
      allF1sq <- c(allF1sq,F1sq)
      allF1sgl <- c(allF1sgl,F1sgl)
      
    }
    allgF1sq <- cbind(allgF1sq,allF1sq)
    allgF1sgl <- cbind(allgF1sgl,allF1sgl)
  }
  
  #plotting
  colnames(allgF1sq) <- colnames(allgF1sgl) <- paste("G",Gs,sep="")
  ng <- length(Gs)
  allgF1s <- matrix(data=NA,nrow=nrep,ncol=2*ng)
  allgF1s[,seq(1,2*ng,by=2)] <- allgF1sq
  allgF1s[,seq(2,2*ng,by=2)] <- allgF1sgl
  colnames(allgF1s) <- paste("G",sort(rep(Gs,2)),sep="")
  maint <- paste("F1 score, p_sel=", ntosel,sep="")
  boxplot(allgF1s,col=rep(c("green","blue"),8), main=maint)
  abline(v=seq(2,2*(ng-1),by=2)+0.5,lty=3)
}



