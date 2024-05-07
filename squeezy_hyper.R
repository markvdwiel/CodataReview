myminML.LA.ridgeGLM <- function (loglambdas, XXblocks, Y, sigmasq = 1, Xunpen = NULL, 
          intrcpt = TRUE, model, minlam = 0, opt.sigma = FALSE, penhyper=FALSE, targlam=1) 
{
  if (model != "linear") {
    opt.sigma <- FALSE
    sigmasq <- 1
  }
  n <- length(Y)
  lambdas <- exp(loglambdas) + minlam
  if (opt.sigma) {
    sigmasq <- lambdas[1]
    lambdas <- lambdas[-1]
  }
  Lam0inv <- squeezy:::.SigmaBlocks(XXblocks, lambdas)
  if (!is.null(Xunpen)) {
    fit <- IWLSridge(XXT = Lam0inv, Y = Y, model = model, 
                     intercept = intrcpt, X1 = Xunpen, maxItr = 500, eps = 10^-12)
  }
  else {
    fit <- IWLSridge(XXT = Lam0inv, Y = Y, model = model, 
                     intercept = intrcpt, maxItr = 500, eps = 10^-12)
  }
  eta <- fit$etas + fit$eta0
  if (model == "linear") {
    mu <- eta
    W <- diag(rep(1, n))
    Hpen <- Lam0inv - Lam0inv %*% solve(solve(W) + Lam0inv, 
                                        Lam0inv)
    t1 <- mvtnorm:::dmvnorm(c(Y), mean = eta, sigma = sigmasq * diag(rep(1, 
                                                               n)), log = TRUE)
  }
  else {
    mu <- 1/(1 + exp(-eta))
    W <- diag(c(mu) * c(1 - mu))
    Hpen <- Lam0inv - Lam0inv %*% solve(diag(1, n) + W %*% 
                                          Lam0inv, W %*% Lam0inv)
    t1 <- sum(Y * log(mu) + (1 - Y) * log(1 - mu))
  }
  t2 <- -1/2/sigmasq * sum(c(Y - mu) * eta)
  t3 <- 1/2 * squeezy:::.logdet(diag(rep(1, n)) - W %*% Hpen)
  
  if(!penhyper) toret <- -(t1 + t2 + t3) else {
    # taus <- sqrt(1/lambdas)
    # tautarg <- sqrt(1/targlam)
    # t4 <- log(2) + sum(dcauchy(taus,0,tautarg, log=TRUE))
    logtarg <- log(targlam)
    #t4 <- sum(dnorm(log(lambdas),logtarg,1,log=TRUE))
    t4 <- sum(log(1/2) - abs(log(lambdas)-logtarg))
    toret <- -(t1 + t2 + t3 + t4)
  }
  return(toret)
}

mydminML.LA.ridgeGLM <- function (loglambdas, XXblocks, Y, sigmasq = 1, Xunpen = NULL, 
          intrcpt = TRUE, model, minlam = 0, opt.sigma = FALSE, penhyper=FALSE,targlam=1) 
{
  if (model != "linear") 
    opt.sigma <- FALSE
  n <- length(Y)
  lambdas <- exp(loglambdas) + minlam
  if (opt.sigma) {
    sigmasq <- lambdas[1]
    lambdas <- lambdas[-1]
  }
  rhos <- log(lambdas)
  Lam0inv <- squeezy:::.SigmaBlocks(XXblocks, lambdas)
  if (!is.null(Xunpen)) {
    fit <- IWLSridge(XXT = Lam0inv, Y = Y, model = model, 
                     intercept = intrcpt, X1 = Xunpen, maxItr = 500, eps = 10^-12)
  }
  else {
    fit <- IWLSridge(XXT = Lam0inv, Y = Y, model = model, 
                     intercept = intrcpt, maxItr = 500, eps = 10^-12)
  }
  eta <- fit$etas + fit$eta0
  if (intrcpt) 
    Xunpen <- cbind(Xunpen, rep(1, n))
  if (model == "linear") {
    mu <- eta
    W <- diag(rep(1, n))
    dwdeta <- rep(0, n)
  }
  else if (model == "logistic") {
    mu <- 1/(1 + exp(-eta))
    W <- diag(c(mu) * c(1 - mu))
    dwdeta <- c(mu * (1 - mu) * (1 - 2 * mu))
  }
  else {
    stop(print(paste("Only model type linear and logistic supported")))
  }
  Hpen <- Lam0inv - Lam0inv %*% solve(diag(1, n) + W %*% Lam0inv, 
                                      W %*% Lam0inv)
  if (!is.null(Xunpen)) {
    WP1W <- diag(rep(1, n)) - Xunpen %*% solve(t(Xunpen) %*% 
                                                 W %*% Xunpen, t(Xunpen) %*% W)
  }
  else {
    WP1W <- diag(rep(1, n))
  }
  L2 <- WP1W %*% (diag(rep(1, n)) - t(solve(diag(1, n) + W %*% 
                                              WP1W %*% Lam0inv, W %*% WP1W %*% Lam0inv)))
  Hj <- lapply(1:length(XXblocks), function(i) {
    L2 %*% XXblocks[[i]]/lambdas[i]
  })
  if (is.null(Xunpen)) {
    Hres <- squeezy:::.Hpen(diag(W), Lam0inv)
    Hmat <- Hres$Hmat
  }
  else {
    Hres <- squeezy:::.Hunpen(diag(W), Lam0inv, Xunpen)
    Hmat <- Hres$Hmat
  }
  detadrho <- lapply(1:length(rhos), function(j) {
    -exp(rhos[j]) * (diag(rep(1, n)) - Hmat %*% W) %*% t(Hj[[j]]/lambdas[j]) %*% 
      c(Y - mu + W %*% eta)
  })
  dWdrho <- lapply(1:length(rhos), function(j) {
    dwdeta * c(detadrho[[j]])
  })
  t1 <- sapply(1:length(rhos), function(j) {
    -1/sigmasq * c((Y - mu)) %*% detadrho[[j]]
  })
  t2 <- sapply(1:length(rhos), function(j) {
    1/2/sigmasq * c(-W %*% eta + Y - mu) %*% detadrho[[j]]
  })
  t3 <- sapply(1:length(rhos), function(j) {
    -0.5 * sum(diag(solve(diag(1, n) + W %*% Lam0inv, W %*% 
                            XXblocks[[j]]/lambdas[j])))
  })
  t4 <- sapply(1:length(rhos), function(j) {
    0.5 * sum(diag(Hpen) * dWdrho[[j]])
  })
  
  if(!penhyper) toret <- (t1 + t2 + t3 + t4) else {
    # t5 <- -sum(((exp(-1/2 * loglambdas))^2/targlam)/
    #              (1+(exp(-1/2 * loglambdas)/targlam)^2))
    logtarg <- log(targlam)
    #t5 <- rhos-logtarg
    #t5 <- 2*(rhos-logtarg >=0)-1 
    #t5 <- rhos-logtarg
    t5 <- (rhos - logtarg)/sqrt((rhos - logtarg)^2 + 0.1)
    toret <- (t1 + t2 + t3 + t4 + t5)
  }
  if (!opt.sigma) return(toret)
  ts.1 <- n/2/sigmasq
  ts.2 <- -1/2/sigmasq^2 * t(Y - eta) %*% Y
 
  return(c(c(ts.1 + ts.2) * sigmasq, toret))
}

squeezyXXbl <- function (X, groupset, unpen = NULL) 
  {
    n <- dim(X)[1]
    p <- dim(X)[2]
    groupsets <- list(groupset)
    penfctr <- rep(1, p)
    if (length(unpen) > 0) {
      penfctr[unpen] <- 0
      if (any(unlist(groupsets) %in% unpen)) {
        warning("Unpenalised covariates removed from group set")
        for (i in 1:length(groupsets)) {
          for (j in 1:length(groupsets[[i]])) {
            if (all(groupsets[[i]][[j]] %in% unpen)) {
              groupsets[[i]][[j]] <- NULL
            }
            else {
              groupsets[[i]][[j]] <- groupsets[[i]][[j]][!(groupsets[[i]][[j]] %in% 
                                                             unpen)]
            }
          }
        }
      }
    }
    G <- sapply(groupsets, length)
    m <- length(G)
    indGrpsGlobal <- list(1:G[1])
    if (m > 1) {
      for (i in 2:m) {
        indGrpsGlobal[[i]] <- (sum(G[1:(i - 1)]) + 1):sum(G[1:i])
      }
    }
    Kg <- lapply(groupsets, function(x) (sapply(x, length)))
    i <- unlist(sapply(1:sum(G), function(x) {
      rep(x, unlist(Kg)[x])
    }))
    j <- unlist(unlist(groupsets))
    ind <- sparseMatrix(i, j, x = 1)
    Ik <- lapply(1:m, function(i) {
      x <- rep(0, sum(G))
      x[(sum(G[1:i - 1]) + 1):sum(G[1:i])] <- 1
      as.vector(x %*% ind)
    })
    Zt <- ind
    if (G[1] > 1) {
      Zt[1:G[1], ] <- t(t(ind[1:G[1], ])/apply(ind[1:G[1], 
      ], 2, sum))
    }
    if (m > 1) {
      for (i in 2:m) {
        if (G[i] > 1) {
          Zt[indGrpsGlobal[[i]], ] <- t(t(ind[indGrpsGlobal[[i]], 
          ])/apply(ind[indGrpsGlobal[[i]], ], 2, sum))
        }
      }
    }
    if (dim(Zt)[2] < p) 
      Zt <- cbind(Zt, matrix(rep(NaN, (p - dim(Zt)[2]) * sum(G)), 
                             c(sum(G), p - dim(Zt)[2])))
    Xxtnd <- do.call(cbind, lapply(groupsets[[1]], function(group) {
      t(t(X[, group])/sqrt(Ik[[1]][group]))
    }))
    Kg2 <- c(1, Kg[[1]])
    G2 <- length(Kg2) - 1
    groupxtnd <- lapply(2:length(Kg2), function(i) {
      sum(Kg2[1:(i - 1)]):(sum(Kg2[1:i]) - 1)
    })
    groupxtnd2 <- unlist(sapply(1:G2, function(x) {
      rep(x, Kg2[x + 1])
    }))
    Xunpen <- NULL
    if (sum((1:p) %in% unpen) > 0) 
      Xunpen <- X[, (1:p) %in% unpen]
    Xbl <- createXblocks(lapply(groupxtnd, function(ind) Xxtnd[, 
                                                               ind, drop = FALSE]))
    XXbl <- createXXblocks(lapply(groupxtnd, function(ind) Xxtnd[, 
                                                                 ind, drop = FALSE]))
   return(XXbl)
}
      

squeezy_hyper <- function (Y, X, groupset, alpha = 1, model = NULL, X2 = NULL, 
          Y2 = NULL, unpen = NULL, intrcpt = TRUE, method = c("ecpcEN", 
                                                              "MML", "MML.noDeriv", "CV"), fold = 10, compareMR = TRUE, 
          selectAIC = FALSE, fit.ecpc = NULL, lambdas = NULL, lambdaglobal = NULL, 
          lambdasinit = NULL, sigmasq = NULL, ecpcinit = TRUE, SANN = FALSE, 
          minlam = 10^-3, standardise_Y = NULL, reCV = NULL, opt.sigma = NULL, 
          resultsAICboth = FALSE, silent = FALSE, penhyper=FALSE, targlam=NULL) 
{
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (length(lambdas) == p && missing(groupset)) 
    groupset <- lapply(1:p, function(x) x)
  groupsets <- list(groupset)
  if (!is.null(X2)) 
    n2 <- dim(X2)[1]
  if (is.null(model)) {
    if (all(is.element(Y, c(0, 1))) || is.factor(Y)) {
      model <- "logistic"
    }
    else if (all(is.numeric(Y)) & !(is.matrix(Y) && dim(Y)[2] == 
                                    2)) {
      model <- "linear"
    }
    else {
      model <- "cox"
    }
  }
  if (length(method) == 4) {
    method <- "MML"
  }
  switch(method, ecpcEN = {
    if (is.null(fit.ecpc)) stop("provide ecpc fit results")
    if (is.null(lambdas)) lambdas <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
                                                          fit.ecpc$tauglobal)
    if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
    if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
    if (is.null(standardise_Y)) standardise_Y <- TRUE
    if (is.null(reCV)) reCV <- TRUE
    if (is.null(opt.sigma)) opt.sigma <- FALSE
  }, MML = {
    if (!is.null(fit.ecpc)) {
      if (is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
                                                                    fit.ecpc$tauglobal)
      if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
      if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
    }
    if (is.null(standardise_Y)) standardise_Y <- FALSE
    if (is.null(opt.sigma)) opt.sigma <- TRUE
    if (is.null(reCV)) {
      reCV <- FALSE
      if (model != "linear" | !opt.sigma) reCV <- TRUE
    }
  }, MML.noDeriv = {
    if (!is.null(fit.ecpc)) {
      if (is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
                                                                    fit.ecpc$tauglobal)
      if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
      if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
    }
    if (is.null(standardise_Y)) standardise_Y <- FALSE
    if (is.null(opt.sigma)) opt.sigma <- TRUE
    if (is.null(reCV)) {
      reCV <- FALSE
      if (!opt.sigma) reCV <- TRUE
    }
  }, CV = {
    if (!is.null(fit.ecpc)) {
      if (is.null(lambdasinit)) lambdasinit <- fit.ecpc$sigmahat/(fit.ecpc$gamma * 
                                                                    fit.ecpc$tauglobal)
      if (is.null(lambdaglobal)) lambdaglobal <- fit.ecpc$sigmahat/fit.ecpc$tauglobal
      if (is.null(sigmasq)) sigmasq <- fit.ecpc$sigmahat
    }
    if (is.null(standardise_Y)) standardise_Y <- FALSE
    if (is.null(opt.sigma)) opt.sigma <- FALSE
    if (is.null(reCV)) reCV <- TRUE
  })
  switch(model, linear = {
    fml <- "gaussian"
    sd_y <- sqrt(var(Y) * (n - 1)/n)[1]
    if (standardise_Y) {
      Y <- Y/sd_y
      sd_y_former <- sd_y
      sd_y <- 1
      if (!is.null(sigmasq)) sigmasq <- sigmasq/sd_y_former
    }
    if (method == "MML") minlam <- max(minlam, 10^-4 * var(Y))
  }, logistic = {
    fml <- "binomial"
    opt.sigma <- FALSE
    standardise_Y <- FALSE
    sd_y <- 1
    sd_y_former <- sd_y
    levelsY <- cbind(c(0, 1), c(0, 1))
    if (!all(is.element(Y, c(0, 1)))) {
      oldLevelsY <- levels(Y)
      levels(Y) <- c("0", "1")
      Y <- as.numeric(Y) - 1
      levelsY <- cbind(oldLevelsY, c(0, 1))
      colnames(levelsY) <- c("Old level names", "New level names")
      if (!is.null(Y2)) {
        levels(Y2) <- c("0", "1")
        Y2 <- as.numeric(Y2) - 1
      }
      if (!silent) print("Y is put in 0/1 format, see levelsY in output for new names")
    }
  }, cox = {
    fml <- "cox"
    opt.sigma <- FALSE
    standardise_Y <- FALSE
    sd_y <- 1
    sd_y_former <- sd_y
    intrcpt <- FALSE
  })
  penfctr <- rep(1, p)
  if (length(unpen) > 0) {
    penfctr[unpen] <- 0
    if (any(unlist(groupsets) %in% unpen)) {
      warning("Unpenalised covariates removed from group set")
      for (i in 1:length(groupsets)) {
        for (j in 1:length(groupsets[[i]])) {
          if (all(groupsets[[i]][[j]] %in% unpen)) {
            groupsets[[i]][[j]] <- NULL
          }
          else {
            groupsets[[i]][[j]] <- groupsets[[i]][[j]][!(groupsets[[i]][[j]] %in% 
                                                           unpen)]
          }
        }
      }
    }
  }
  G <- sapply(groupsets, length)
  m <- length(G)
  indGrpsGlobal <- list(1:G[1])
  if (m > 1) {
    for (i in 2:m) {
      indGrpsGlobal[[i]] <- (sum(G[1:(i - 1)]) + 1):sum(G[1:i])
    }
  }
  Kg <- lapply(groupsets, function(x) (sapply(x, length)))
  i <- unlist(sapply(1:sum(G), function(x) {
    rep(x, unlist(Kg)[x])
  }))
  j <- unlist(unlist(groupsets))
  ind <- sparseMatrix(i, j, x = 1)
  Ik <- lapply(1:m, function(i) {
    x <- rep(0, sum(G))
    x[(sum(G[1:i - 1]) + 1):sum(G[1:i])] <- 1
    as.vector(x %*% ind)
  })
  Zt <- ind
  if (G[1] > 1) {
    Zt[1:G[1], ] <- t(t(ind[1:G[1], ])/apply(ind[1:G[1], 
    ], 2, sum))
  }
  if (m > 1) {
    for (i in 2:m) {
      if (G[i] > 1) {
        Zt[indGrpsGlobal[[i]], ] <- t(t(ind[indGrpsGlobal[[i]], 
        ])/apply(ind[indGrpsGlobal[[i]], ], 2, sum))
      }
    }
  }
  if (dim(Zt)[2] < p) 
    Zt <- cbind(Zt, matrix(rep(NaN, (p - dim(Zt)[2]) * sum(G)), 
                           c(sum(G), p - dim(Zt)[2])))
  Xxtnd <- do.call(cbind, lapply(groupsets[[1]], function(group) {
    t(t(X[, group])/sqrt(Ik[[1]][group]))
  }))
  Kg2 <- c(1, Kg[[1]])
  G2 <- length(Kg2) - 1
  groupxtnd <- lapply(2:length(Kg2), function(i) {
    sum(Kg2[1:(i - 1)]):(sum(Kg2[1:i]) - 1)
  })
  groupxtnd2 <- unlist(sapply(1:G2, function(x) {
    rep(x, Kg2[x + 1])
  }))
  Xunpen <- NULL
  if (sum((1:p) %in% unpen) > 0) 
    Xunpen <- X[, (1:p) %in% unpen]
  Xbl <- createXblocks(lapply(groupxtnd, function(ind) Xxtnd[, 
                                                             ind, drop = FALSE]))
  XXbl <- createXXblocks(lapply(groupxtnd, function(ind) Xxtnd[, 
                                                               ind, drop = FALSE]))
  if (is.null(lambdaglobal)) {
    lambdaseq <- 10^c(-10:10)
    XXbl1 <- list(apply(simplify2array(XXbl), c(1, 2), sum))
    ML <- sapply(log(lambdaseq), function(lam) {
      temp <- try(myminML.LA.ridgeGLM(loglambdas = lam, XXblocks = XXbl1, 
                                    sigmasq = sd_y^2, Y = Y, Xunpen = Xunpen, intrcpt = intrcpt, 
                                    model = model, penhyper=FALSE), silent = TRUE)
      if (class(temp)[1] != "try-error") {
        return(temp)
      }
      else return(NaN)
    })
    if (all(is.nan(ML))) 
      stop("Error in estimating global lambda, try standardising data")
    lambdaseq <- lambdaseq[!is.nan(ML) & (ML != Inf)]
    ML <- ML[!is.nan(ML) & (ML != Inf)]
    se <- abs(ML[1] - rev(ML)[1])/100
    if (ML[1] < rev(ML)[1]) {
      lambda <- rev(lambdaseq)[which.min(rev(ML))]
      if (lambda == lambdaseq[1]) {
        lambda <- max(lambdaseq[ML <= (min(ML) + se)])
      }
    }
    else {
      lambda <- lambdaseq[which.min(ML)]
      if (lambda == rev(lambdaseq)[1]) {
        lambda <- min(lambdaseq[ML <= (min(ML) + se)])
      }
    }
  } else {
    lambda <- lambdaglobal
  }
  if (selectAIC | compareMR) {
    if (method == "CV") {
      leftout <- CVfolds(Y = Y, kfold = fold, nrepeat = 3, 
                         fixedfolds = FALSE)
      lambda1groupfit <- optLambdasWrap(penaltiesinit = lambda, 
                                        XXblocks = list(apply(simplify2array(XXbl), c(1, 
                                                                                      2), sum)), Y = Y, folds = leftout, X1 = Xunpen, 
                                        intercept = intrcpt, score = ifelse(model == 
                                                                              "linear", "mse", "loglik"), model = model, 
                                        maxItropt2 = 500, reltol = 10^-3)
      lambda <- lambda1groupfit$optpen
    }
    else if (method == "MML.noDeriv") {
      lambda1groupfit <- optLambdas_mgcv(penaltiesinit = lambda, 
                                         Y = Y, XXblocks = list(apply(simplify2array(XXbl), 
                                                                      c(1, 2), sum)), model = model, reltol = 0.001, 
                                         maxItropt = 500, tracescore = FALSE, fixedseed = FALSE, 
                                         optmethod = "Nelder-Mead", sigmasq = sigmasq)
      lambda <- lambda1groupfit$optpen
    }
    else if (method == "MML") {
      sigmahat <- sd_y
      if (!is.null(sigmasq) & model == "linear") 
        sigmahat <- sigmasq
      if (opt.sigma) {
        lambda1groupfit <- optim(par = c(log(sigmahat), 
                                         log(lambda)), fn = myminML.LA.ridgeGLM, gr = mydminML.LA.ridgeGLM, 
                                 method = "BFGS", minlam = minlam, XXblocks = list(apply(simplify2array(XXbl), 
                                                                                         c(1, 2), sum)), Y = Y, sigmasq = sigmahat, 
                                 model = model, intrcpt = intrcpt, Xunpen = Xunpen, 
                                 opt.sigma = opt.sigma, penhyper=penhyper)
        sigmasq <- exp(lambda1groupfit$par[1]) + minlam
        lambda <- exp(lambda1groupfit$par[-1]) + minlam
      }
      else {
        lambda1groupfit <- optim(par = log(lambda), fn = myminML.LA.ridgeGLM, 
                                 gr = mydminML.LA.ridgeGLM, method = "BFGS", minlam = minlam, 
                                 XXblocks = list(apply(simplify2array(XXbl), 
                                                       c(1, 2), sum)), Y = Y, sigmasq = sigmahat, 
                                 model = model, intrcpt = intrcpt, Xunpen = Xunpen, 
                                 opt.sigma = opt.sigma, penhyper=penhyper)
        lambda <- exp(lambda1groupfit$par) + minlam
      }
    }
  }
  sigmahat <- sd_y
  if (model == "linear") {
    if (!is.null(sigmasq)) 
      sigmahat <- sigmasq
    else {
      XXT1grp <- SigmaFromBlocks(XXblocks = list(apply(simplify2array(XXbl), 
                                                       c(1, 2), sum)), lambda)
      if (length(unpen) > 0 | intrcpt) {
        Xunpen2 <- Xunpen
        if (intrcpt) 
          Xunpen2 <- cbind(Xunpen, rep(1, n))
        if (intrcpt && length(unpen) == 0) {
          betaunpenML1grp <- sum(Y)/n
        }
        else {
          temp <- solve(XXT1grp + diag(rep(1, n)), Xunpen2)
          betaunpenML1grp <- solve(t(Xunpen2) %*% temp, 
                                   t(temp) %*% Y)
        }
        sigmahat <- c(t(Y - Xunpen2 %*% betaunpenML1grp) %*% 
                        solve(XXT1grp + diag(rep(1, n)), Y - Xunpen2 %*% 
                                betaunpenML1grp)/n)
      }
      else {
        sigmahat <- c(t(Y) %*% solve(XXT1grp + diag(rep(1, 
                                                        n)), Y)/n)
      }
    }
  }
  if (selectAIC) {
    sigmahat1group <- sigmahat
  }
  if (is.null(lambdas)) {
    if (is.null(lambdasinit) | !ecpcinit) {
      lambdasinit <- rep(lambda, G)
    }
    if (any(lambdasinit == Inf)) {
      lambdasinit[lambdasinit > 2 * lambda] <- 2 * lambda
    }
    if (method == "CV") {
      leftout <- CVfolds(Y = Y, kfold = fold, nrepeat = 3, 
                         fixedfolds = FALSE)
      jointlambdas <- optLambdasWrap(penaltiesinit = lambdasinit, 
                                     XXblocks = XXbl, Y = Y, folds = leftout, X1 = Xunpen, 
                                     intercept = intrcpt, score = ifelse(model == 
                                                                           "linear", "mse", "loglik"), model = model, 
                                     maxItropt2 = 500, reltol = 10^-3, traceCV = FALSE)
      lambdas <- jointlambdas$optpen
    }
    else if (method == "MML.noDeriv") {
      if (ecpcinit) {
        if (SANN) {
          jointlambdas <- optLambdas_mgcvWrap(penaltiesinit = lambdasinit, 
                                              XXblocks = XXbl, Y = Y, model = model, reltol = 1e-04, 
                                              maxItropt2 = 1000, tracescore = FALSE, fixedseed = FALSE, 
                                              optmethod2 = "Nelder-Mead", sigmasq = sigmahat, 
                                              opt.sigma = opt.sigma)
        }
        else {
          jointlambdas <- optLambdas_mgcv(penaltiesinit = lambdasinit, 
                                          XXblocks = XXbl, Y = Y, model = model, reltol = 1e-04, 
                                          maxItropt = 1000, tracescore = FALSE, fixedseed = FALSE, 
                                          optmethod = "Nelder-Mead", sigmasq = sigmahat, 
                                          opt.sigma = opt.sigma)
        }
      }
      else {
        if (SANN) {
          jointlambdas <- optLambdas_mgcvWrap(penaltiesinit = rep(lambda, 
                                                                  G), XXblocks = XXbl, Y = Y, model = model, 
                                              reltol = 1e-04, maxItropt2 = 1000, tracescore = FALSE, 
                                              fixedseed = FALSE, optmethod2 = "Nelder-Mead", 
                                              sigmasq = sigmahat, opt.sigma = opt.sigma)
        }
        else {
          jointlambdas <- optLambdas_mgcv(penaltiesinit = rep(lambda, 
                                                              G), XXblocks = XXbl, Y = Y, model = model, 
                                          reltol = 1e-04, maxItropt = 1000, tracescore = FALSE, 
                                          fixedseed = FALSE, optmethod = "Nelder-Mead", 
                                          sigmasq = sigmahat, opt.sigma = opt.sigma)
        }
      }
      if (opt.sigma) {
        sigmasq <- jointlambdas$optpen[1]
        lambdas <- jointlambdas$optpen[-1]
      }
      else {
        lambdas <- jointlambdas$optpen
      }
    }
    else if (method == "MML") {
      if (opt.sigma) {
        jointlambdas <- optim(par = c(log(sigmahat), 
                                      log(lambdasinit)), fn = myminML.LA.ridgeGLM, 
                              gr = mydminML.LA.ridgeGLM, method = "BFGS", minlam = minlam, 
                              XXblocks = XXbl, Y = Y, opt.sigma = opt.sigma, 
                              model = model, intrcpt = intrcpt, Xunpen = Xunpen, penhyper=penhyper, targlam = targlam)
        sigmasq <- exp(jointlambdas$par[1]) + minlam
        lambdas <- exp(jointlambdas$par[-1]) + minlam
      }
      else {
        jointlambdas <- optim(par = log(lambdasinit), 
                              fn = myminML.LA.ridgeGLM, gr = mydminML.LA.ridgeGLM, 
                              method = "BFGS", minlam = minlam, XXblocks = XXbl, 
                              Y = Y, sigmasq = sigmahat, model = model, intrcpt = intrcpt, 
                              Xunpen = Xunpen,penhyper=penhyper, targlam = targlam)
        lambdas <- exp(jointlambdas$par) + minlam
      }
    }
  }
  else {
    lambdasinit <- lambdas
  }
  sigmahat <- 1
  if (model == "linear") {
    if (!is.null(sigmasq)) 
      sigmahat <- sigmasq
    else {
      XXT <- SigmaFromBlocks(XXblocks = XXbl, lambdas)
      if (length(unpen) > 0 | intrcpt) {
        Xunpen2 <- Xunpen
        if (intrcpt) 
          Xunpen2 <- cbind(Xunpen, rep(1, n))
        if (intrcpt && length(unpen) == 0) {
          betaunpenML <- sum(Y)/n
        }
        else {
          temp <- solve(XXT + diag(rep(1, n)), Xunpen2)
          betaunpenML <- solve(t(Xunpen2) %*% temp, t(temp) %*% 
                                 Y)
        }
        sigmahat <- c(t(Y - Xunpen2 %*% betaunpenML) %*% 
                        solve(XXT + diag(rep(1, n)), Y - Xunpen2 %*% 
                                betaunpenML)/n)
      } else {
        sigmahat <- c(t(Y) %*% solve(XXT + diag(rep(1, 
                                                    n)), Y)/n)
      }
    }
  }
  
  MLinit <- myminML.LA.ridgeGLM(loglambdas = log(lambdasinit), 
                              opt.sigma = FALSE, sigmasq = sigmahat, XXblocks = XXbl, 
                              Y = Y, model = model, intrcpt = intrcpt, minlam = 0, penhyper=penhyper, targlam = targlam)
  MLfinal <- myminML.LA.ridgeGLM(loglambdas = log(lambdas), opt.sigma = FALSE, 
                               sigmasq = sigmahat, XXblocks = XXbl, Y = Y, model = model, 
                               intrcpt = intrcpt, minlam = 0, penhyper = penhyper, targlam = targlam)
  # MLfinal2 <- myminML.LA.ridgeGLM(loglambdas = log(lambdas), opt.sigma = FALSE, 
  #                     sigmasq = sigmahat, XXblocks = XXbl, Y = Y, model = model, 
  #                     intrcpt = intrcpt, minlam = 0, penhyper=TRUE, targlam = targlam)
  
  
  
  if (selectAIC) {
    lambda1group <- lambda
    if (sum((1:p) %in% unpen) > 0) {
      AIC1group <- mAIC.LA.ridgeGLM(log(lambda1group), 
                                    XXblocks = list(apply(simplify2array(XXbl), c(1, 
                                                                                  2), sum)), Y = Y, sigmasq = sigmahat1group, 
                                    Xunpen = Xunpen, intrcpt = intrcpt, model = model)
      AICmultigroup <- mAIC.LA.ridgeGLM(log(lambdas), XXblocks = XXbl, 
                                        Y = Y, sigmasq = sigmahat, Xunpen = Xunpen, intrcpt = intrcpt, 
                                        model = model)
    }
    else {
      AIC1group <- mAIC.LA.ridgeGLM(log(lambda1group), 
                                    XXblocks = list(apply(simplify2array(XXbl), c(1, 
                                                                                  2), sum)), Y = Y, sigmasq = sigmahat1group, 
                                    intrcpt = intrcpt, model = model)
      AICmultigroup <- mAIC.LA.ridgeGLM(log(lambdas), XXblocks = XXbl, 
                                        Y = Y, sigmasq = sigmahat, intrcpt = intrcpt, 
                                        model = model)
    }
    if (AIC1group <= AICmultigroup) {
      lambdasNotOptimalAIC <- lambdas
      sigmahatNotOptimalAIC <- sigmahat
      lambdas <- rep(lambda1group, G)
      sigmahat <- sigmahat1group
      modelbestAIC <- "onegroup"
    }
    else {
      lambdasNotOptimalAIC <- lambda1group
      sigmahatNotOptimalAIC <- sigmahat1group
      modelbestAIC <- "multigroup"
    }
  }
  tauglobal <- sigmahat/lambda
  gamma <- lambda/lambdas
  lambdap <- lambda/(as.vector(gamma %*% Zt))
  lambdap[lambdap < 0] <- Inf
  glmGR <- NA
  if (compareMR || alpha == 0) {
    XXT <- SigmaFromBlocks(XXbl, penalties = lambdas)
    if (sum((1:p) %in% unpen) > 0) {
      fit <- IWLSridge(XXT, Y = Y, model = model, intercept = intrcpt, 
                       X1 = X[, (1:p) %in% unpen])
    }
    else {
      fit <- IWLSridge(XXT, Y = Y, model = model, intercept = intrcpt)
    }
    betafit <- betasout(fit, Xblocks = Xbl, penalties = lambdas)
    a0MR <- 0
    if (intrcpt) 
      a0MR <- c(betafit[[1]][1])
    betaMR <- rep(0, p)
    betaMR[(1:p) %in% unpen] <- betafit[[1]][-1]
    for (i in 1:length(groupsets[[1]])) {
      betaMR[groupsets[[1]][[i]]] <- betaMR[groupsets[[1]][[i]]] + 
        betafit[[1 + i]]/sqrt(Ik[[1]][groupsets[[1]][[i]]])
    }
    rm(betafit)
  }
  pen <- which(!((1:p) %in% unpen))
  if (any(is.nan(sqrt(lambdap[pen])))) {
    browser()
  }
  if (alpha <= 1) {
    varFunc <- function(tauEN, alpha = alpha, tauR) {
      t2 <- -alpha/2/(1 - alpha)^(3/2) * sqrt(tauEN) * 
        exp(dnorm(alpha/2/sqrt(tauEN)/sqrt(1 - alpha), 
                  log = TRUE) - pnorm(-alpha/2/sqrt(tauEN)/sqrt(1 - 
                                                                  alpha), log.p = TRUE))
      varBeta <- tauEN/(1 - alpha) + alpha^2/4/(1 - alpha)^2 + 
        t2
      f <- varBeta - tauR
      return(f)
    }
    lamEN <- function(alpha, tauR) {
      if (alpha == 0) {
        lamEN <- sigmahat/tauR
      }
      else if (alpha == 1) {
        lamEN <- sigmahat * sqrt(8/tauR)
      }
      else if (tauR/sigmahat > 10^6) {
        lamEN <- sigmahat/tauR
        if (alpha > 0.2) {
          ub2 <- 10^6 * sigmahat
          lb2 <- sqrt(10^6 * sigmahat/8)
          temp <- try(uniroot(varFunc, c(0.9 * lb2, 1.1 * 
                                           ub2), alpha = alpha, tauR = 10^6 * sigmahat, 
                              tol = 10^-6))
          if (class(temp)[1] != "try-error") {
            if (temp$root < lamEN) 
              lamEN <- sigmahat/temp$root
          }
        }
      }
      else if (tauR/sigmahat < 10^-6) {
        lamEN <- sigmahat * sqrt(8/tauR)
        if (alpha < 0.8) {
          ub2 <- sqrt(10^-7 * sigmahat/8)
          lb2 <- 10^-7 * sigmahat
          temp <- try(uniroot(varFunc, c(0.9 * lb2, 1.1 * 
                                           ub2), alpha = alpha, tauR = 10^-7 * sigmahat, 
                              tol = 10^-6))
          if (class(temp)[1] != "try-error") {
            if (temp$root > lamEN) 
              lamEN <- temp$root
          }
        }
      }
      else {
        lb <- min(tauR, sqrt(tauR/8))
        ub <- max(tauR, sqrt(tauR/8))
        temp <- try(uniroot(varFunc, c(0.9 * lb, 1.1 * 
                                         ub), alpha = alpha, tauR = tauR, tol = 10^-6))
        if (class(temp)[1] == "try-error") {
          if (tauR < 0.5) 
            lamEN <- sigmahat * sqrt(8/tauR)
          if (tauR > 0.5) 
            lamEN <- sigmahat/tauR
        }
        else {
          lamEN <- sigmahat/temp$root
        }
      }
      return(lamEN)
    }
    alphat <- 1/(1 + 2 * sd_y * (1 - alpha)/alpha)
    lambdasEN <- sapply(sigmahat/lambdas, function(tau) lamEN(tauR = tau, 
                                                              alpha = alpha))
    tauEN <- sigmahat/lambdasEN
    lambdasENhat <- lambdasEN/2 * (alpha/sd_y + 2 * (1 - 
                                                       alpha))
    uniqueTaus <- unique(as.vector(c(sigmahat/lambdas) %*% 
                                     Zt))
    if (dim(X)[2] == dim(Xxtnd)[2]) {
      lambdap <- (as.vector(lambdasENhat %*% Zt))
      lambdapApprox <- (as.vector(lambdasEN %*% Zt))
    }
    else {
      lambdasENunique <- sapply(uniqueTaus, function(tau) lamEN(tauR = tau, 
                                                                alpha = alpha))
      lambdasENuniquehat <- lambdasENunique/2 * (alpha/sd_y + 
                                                   2 * (1 - alpha))
      indUnique <- sapply(as.vector(c(sigmahat/lambdas) %*% 
                                      Zt), function(x) which(uniqueTaus == x))
      lambdap <- lambdasENuniquehat[indUnique]
      lambdapApprox <- lambdasENunique[indUnique]
    }
    if (alpha == 0) {
      beta <- betaMR
      a0 <- a0MR
      lambdaMR_reCV <- lambdas
    }
    else {
      penfctr2 <- penfctr
      penfctr2[pen] <- lambdap[pen]
      not0 <- which(penfctr2 != Inf)
      lambdaEN <- sum(penfctr2[not0])/length(penfctr2[not0])
      penfctr2 <- penfctr2/lambdaEN
      if (model == "cox") {
        glmGR <- glmnet(X[, not0], Y, alpha = alphat, 
                        family = fml, standardize = FALSE, penalty.factor = penfctr2[not0], 
                        thresh = 10^-10)
      }
      else {
        glmGR <- glmnet(X[, not0], Y, alpha = alphat, 
                        family = fml, intercept = intrcpt, standardize = FALSE, 
                        penalty.factor = penfctr2[not0], thresh = 10^-10)
      }
      if (reCV) {
        glmGR.cv <- cv.glmnet(X[, not0], Y, alpha = alphat, 
                              family = fml, intercept = intrcpt, standardize = FALSE, 
                              penalty.factor = penfctr2[not0], thresh = 10^-10)
        sopt <- glmGR.cv$lambda.min
        tempItr <- 1
        while (glmGR.cv$lambda.min == glmGR.cv$lambda[1] & 
               tempItr <= 10) {
          glmGR.cv <- cv.glmnet(X[, not0], Y, alpha = alphat, 
                                family = fml, intercept = intrcpt, standardize = FALSE, 
                                penalty.factor = penfctr2[not0], thresh = 10^-10)
          sopt <- glmGR.cv$lambda.min
          tempItr <- tempItr + 1
        }
        lambdasEN_reCV <- sopt/lambdaEN * n/sd_y * lambdasEN
        lambdaMR_reCV <- sigmahat/sapply(sigmahat/lambdasEN_reCV, 
                                         varFunc, alpha = alpha, tauR = 0)
      }
      else {
        sopt <- lambdaEN/n * sd_y
        lambdaMR_reCV <- lambdas
      }
      temp <- coef(glmGR, s = sopt, exact = TRUE, x = X[, 
                                                        not0], y = Y, alpha = alphat, penalty.factor = penfctr2[not0], 
                   family = fml, intercept = intrcpt)
      beta <- rep(0, p)
      beta[not0] <- temp[-1]
      a0 <- temp[1]
    }
  }
  else {
    stop("alpha should be between 0 and 1")
  }
  if (standardise_Y) {
    beta <- beta * sd_y_former
    a0 <- a0 * sd_y_former
    if (compareMR) {
      betaMR <- betaMR * sd_y_former
      a0MR <- a0MR * sd_y_former
    }
  }
  if (!is.null(X2)) {
    if (compareMR) {
      if (model == "linear") {
        X2c <- cbind(X2, rep(1, n2))
        YpredMR <- X2c %*% c(betaMR, a0MR)
        MSEMR <- sum((YpredMR - Y2)^2)/n2
      }
      if (model == "logistic") {
        X2c <- cbind(X2, rep(1, n2))
        YpredMR <- 1/(1 + exp(-X2c %*% c(betaMR, a0MR)))
        MSEMR <- sum((YpredMR - Y2)^2)/n2
      }
      else if (model == "cox") {
        expXb <- exp(X %*% c(betaMR))
        h0 <- sapply(1:length(Y[, 1]), function(i) {
          Y[i, 2]/sum(expXb[Y[, 1] >= Y[i, 1]])
        })
        H0 <- sapply(Y2[, 1], function(Ti) {
          sum(h0[Y[, 1] <= Ti])
        })
        YpredMR <- H0 * exp(X2 %*% betaMR)
        MSEMR <- sum((YpredMR - Y2[, 2])^2)/n2
      }
    }
    if (model == "linear") {
      X2c <- cbind(X2, rep(1, n2))
      YpredApprox <- X2c %*% c(beta, a0)
      MSEApprox <- sum((YpredApprox - Y2)^2)/n2
    }
    if (model == "logistic") {
      X2c <- cbind(X2, rep(1, n2))
      YpredApprox <- 1/(1 + exp(-X2c %*% c(beta, a0)))
      MSEApprox <- sum((YpredApprox - Y2)^2)/n2
    }
    else if (model == "cox") {
      expXb <- exp(X %*% c(beta))
      h0 <- sapply(1:length(Y[, 1]), function(i) {
        Y[i, 2]/sum(expXb[Y[, 1] >= Y[i, 1]])
      })
      H0 <- sapply(Y2[, 1], function(Ti) {
        sum(h0[Y[, 1] <= Ti])
      })
      YpredApprox <- H0 * exp(X2 %*% beta)
      MSEApprox <- sum((YpredApprox - Y2[, 2])^2)/n2
    }
  }
  output <- list(betaApprox = beta, a0Approx = a0, lambdaApprox = lambdasEN, 
                 lambdapApprox = lambdapApprox, tauMR = sigmahat/lambdas, 
                 lambdaMR = lambdas, lambdaglobal = lambda, sigmahat = sigmahat, 
                 MLinit = MLinit, MLfinal = MLfinal, alpha = alpha, glmnet.fit = glmGR, 
                 lambdaMR_reCV = lambdaMR_reCV)
  if (compareMR) {
    output$betaMR <- betaMR
    output$a0MR <- a0MR
  }
  if (!is.null(X2)) {
    output$YpredApprox <- YpredApprox
    output$MSEApprox <- MSEApprox
    if (compareMR) {
      output$YpredMR <- YpredMR
      output$MSEMR <- MSEMR
    }
  }
  if (selectAIC) {
    output$AICmodels <- list(multigroup = list(lambdas = lambdas, 
                                               sigmahat = sigmahat, AIC = AICmultigroup), onegroup = list(lambdas = lambda1group, 
                                                                                                          sigmahat = sigmahat1group, AIC = AIC1group))
    if (modelbestAIC != "multigroup") {
      output$AICmodels$multigroup$lambdas <- lambdasNotOptimalAIC
      output$AICmodels$multigroup$sigmahat <- sigmahatNotOptimalAIC
    }
    if (resultsAICboth) {
      output$AICmodels$onegroup$fit <- squeezy(Y, X, groupset, 
                                               alpha = alpha, model = model, X2 = X2, Y2 = Y2, 
                                               unpen = unpen, intrcpt = intrcpt, method = "MML", 
                                               fold = fold, compareMR = compareMR, selectAIC = FALSE, 
                                               fit.ecpc = NULL, lambdas = rep(lambda1group, 
                                                                              G), lambdaglobal = lambda1group, sigmasq = sigmahat1group, 
                                               standardise_Y = standardise_Y, reCV = reCV, resultsAICboth = FALSE)
      output$AICmodels$multigroup$fit <- squeezy(Y, X, 
                                                 groupset, alpha = alpha, model = model, X2 = X2, 
                                                 Y2 = Y2, unpen = unpen, intrcpt = intrcpt, method = "MML", 
                                                 fold = fold, compareMR = compareMR, selectAIC = FALSE, 
                                                 fit.ecpc = NULL, lambdas = output$AICmodels$multigroup$lambdas, 
                                                 lambdaglobal = lambda1group, sigmasq = output$AICmodels$multigroup$sigmahat, 
                                                 standardise_Y = standardise_Y, reCV = reCV, resultsAICboth = FALSE)
    }
    output$modelbestAIC <- modelbestAIC
  }
  return(output)
}