assessPOE <- function(Y,gmm,gcc,X,test.locus,hap,f,method,inter){
  X <- as.matrix(X)
  hap <- as.matrix(hap)
  gmm <- as.matrix(gmm)
  gcc <- as.matrix(gcc)
  if(method == 'ROB-HAP'){ #data without NA or delete
    gmmcc <- cbind(gmm,gcc)
    keep <- !apply(is.na(gmmcc) , 1, any)
    Xt <- as.matrix(X[keep,])
    Yt <- Y[keep]
    gmmt <- as.matrix(gmm[keep,])
    gcct <- as.matrix(gcc[keep,])
    
    sp <- NULL
    for(i in 1:nrow(gmmt))
      sp <- c(sp,is.null(haplo(hap,gmmt[i,],gcct[i,])))
    Xtt <- as.matrix(Xt[!sp,])
    Ytt <- Yt[!sp]
    gmmtt <- as.matrix(gmmt[!sp,])
    gcctt <- as.matrix(gcct[!sp,])
    
    if(inter == FALSE)
      fit <- imprinting_robhap(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    if(inter == TRUE)
      fit <- imprinting_robhap.inter(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    return(list(est = fit$est,sd = fit$sd,est.log = fit$est.log,sd.log = fit$sd.log))
  }
  
  if(method == 'ROB-HAP-EM'){ #method incorporate NA
    keep <- is.na(gmm[,test.locus])
    Xt <- as.matrix(X[!keep,])
    Yt <- Y[!keep]
    gmmt <- as.matrix(gmm[!keep,])
    gcct <- as.matrix(gcc[!keep,])
    
    sp <- NULL
    for(i in 1:nrow(gmmt))
      sp <- c(sp,is.null(haplo.na(hap,gmmt[i,],gcct[i,])))
    Xtt <- as.matrix(Xt[!sp,])
    Ytt <- Yt[!sp]
    gmmtt <- as.matrix(gmmt[!sp,])
    gcctt <- as.matrix(gcct[!sp,])
    
    if(inter == FALSE)
      fit <- imprinting_robhap.na(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    if(inter == TRUE)
      fit <- imprinting_robhap.na.inter(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    return(list(est = fit$est,sd = fit$sd,est.log = fit$est.log,sd.log = fit$sd.log))
  }
  
  
  if(method == 'IND-HAP'){
    gmmcc <- cbind(gmm,gcc)
    keep <- !apply(is.na(gmmcc) , 1, any)
    Xt <- as.matrix(X[keep,])
    Yt <- Y[keep]
    gmmt <- as.matrix(gmm[keep,])
    gcct <- as.matrix(gcc[keep,])
    
    sp <- NULL
    for(i in 1:nrow(gmmt))
      sp <- c(sp,is.null(haplo(hap,gmmt[i,],gcct[i,])))
    Xtt <- as.matrix(Xt[!sp,])
    Ytt <- Yt[!sp]
    gmmtt <- as.matrix(gmmt[!sp,])
    gcctt <- as.matrix(gcct[!sp,])
    
    if(inter == FALSE)
      fit <- imprinting_indhap(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    if(inter == TRUE)
      fit <- imprinting_indhap.inter(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    return(list(est = fit$est,sd = fit$sd,est.log = fit$est.log,sd.log = fit$sd.log))
  }
  
  if(method == 'IND-HAP-EM'){
    keep <- is.na(gmm[,test.locus])
    Xt <- as.matrix(X[!keep,])
    Yt <- Y[!keep]
    gmmt <- as.matrix(gmm[!keep,])
    gcct <- as.matrix(gcc[!keep,])
    
    sp <- NULL
    for(i in 1:nrow(gmmt))
      sp <- c(sp,is.null(haplo.na(hap,gmmt[i,],gcct[i,])))
    Xtt <- as.matrix(Xt[!sp,])
    Ytt <- Yt[!sp]
    gmmtt <- as.matrix(gmmt[!sp,])
    gcctt <- as.matrix(gcct[!sp,])
    
    if(inter == FALSE)
      fit <- imprinting_indhap.na(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    if(inter == TRUE)
      fit <- imprinting_indhap.na.inter(Ytt,gmmtt,gcctt,Xtt,test.locus,hap,f)
    return(list(est = fit$est,sd = fit$sd,est.log = fit$est.log,sd.log = fit$sd.log))
  }
}

est.haplo <- function(gmm){
  keep <- !apply(is.na(gmm),1,any)
  gmmt <- as.matrix(gmm[keep,])
  
  library(haplo.stats)
  gmmt2 <- geno1to2(gmmt)
  haplo.raw <- haplo.em(gmmt2)
  haplo.rm <- haplo.raw$haplotype[haplo.raw$hap.prob > 0.01,]
  hapt <- matrix(as.numeric(as.matrix(haplo.rm)), nrow = nrow(haplo.rm), byrow = FALSE) - 1
  return(hapt)
}

summary.POE <- function(EST,SD){
  pvalue <- 2 * (1 - pnorm(abs(EST), mean = 0, sd = SD))
  return(round(c(EST,SD,EST - 1.96*SD,EST + 1.96*SD,pvalue),3))
}

#ROB-HAP & LOGIT-HAP
imprinting_robhap <- function(Y,gmm,gcc,X,loci,hap,f){
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  K <- ncol(hap)
  nX <- max(1,ncol(X))
  
  #initial value
  F <- distinguish(gmm,gcc,loci,hap)
  gm <- F$gm
  gc <- F$gc
  gcm <- F$gcm
  gcp <- F$gcp
  phi <- F$phi
  Z <- design.matrix(gm,gc,gcm,gcp,X,phi)
  fit <- glm(Y ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  sd.log <- as.vector(res[,2])
  
  if(ncol(gmm) > 1){
    library(haplo.stats)
    gmm2 <- geno1to2(gmm)
    haplo_raw <- haplo.em(gmm2)
    ppi_raw <- haplo_raw$hap.prob
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)),n_hap_raw,K,byrow = FALSE) - 1
  }
  else{
    p <- mean(gmm)/2
    ppi_raw <- c(1-p,p)
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(c(0,1),n_hap_raw,K)
  }
  
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1]
  ind_hh <- mul_whichrow(hh_raw,hap)
  ppi_ini <- ppi_raw[ind_hh]
  ppi_ini <- ppi_ini/sum(ppi_ini)
  theta0 <- mul_logitf(ppi_ini)
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.robust(c(beta0,theta),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.robust(c(beta,theta0),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.robust(para,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX)],mul_logistic(est[(5+nX):(3+nX+n_hap)])) 
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

#likelihood function for ROB-HAP
likeli0.robust <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  beta0 <- para0[1:(4+nX)]
  mu0 <- mul_logistic(para0[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(h in hh){
    u <- u + (log(P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.rob(h,gmm,gcc,hap,mu))  - log(P_g.rob(gmm[loci],theta))) * P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.robust <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.robust(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.robust(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.robust <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  for(h in hh)
    u <- u + P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.rob(h,gmm,gcc,hap,mu) / P_g.rob(gmm[loci],theta)
  return(u)
}

likelioriginal.robust <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.robust(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.robust(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

#some functions

H.robust <- function(gm,X,beta,theta){
  if(gm == 0)
    res <- P_Y1(1,0,0,0,X,beta) * (1-theta) + P_Y1(1,0,1,-1,X,beta) * theta
  else if(gm == 1)
    res <- 0.5* (P_Y1(1,1,0,0,X,beta) * (1-theta) + P_Y1(1,1,1,-1,X,beta) * theta + P_Y1(1,1,1,1,X,beta) * (1-theta) + P_Y1(1,1,2,0,X,beta) * theta)
  else if(gm == 2)
    res <- P_Y1(1,2,1,1,X,beta) * (1-theta) + P_Y1(1,2,2,0,X,beta) * theta
  return(res)
}

P_H.rob <- function(h,gmm,gcc,hap,mu){
  hm1 <- hap[h,]
  hm2 <- gmm - hm1
  hc2 <- gcc - hm1
  h2 <- num_hap(hm2,hap)
  h3 <- num_hap(hc2,hap)
  return(mu[h]*mu[h2]*mu[h3])
}

P_Y.rob <- function(Y,X,h,gmm,gcc,hap,loci,beta){
  hm1 <- hap[h,]
  tmp <- exp(sum(c(1,gmm[loci],gcc[loci],2*hm1[loci]-gcc[loci],X) * beta))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

P_g.rob <- function(g,theta){
  p <- c((1-theta) ^ 2,2 * (1-theta) * theta,(theta ^ 2))
  return(p[g+1])
}

num_hap <- function(h,hap){
  K <- nrow(hap)
  a <- 0
  for(k in 1:K)
    if(all(h == hap[k,]))
      a <- k
  return(a)
}

haplo <- function(hap,gmm,gcc){ #for one person to test hap
  h <- NULL
  K <- nrow(hap)
  for(k in 1:K){
    hm1 <- hap[k,] #hap for child
    hm2 <- gmm - hm1
    hc2 <- gcc - hm1
    r <- 0
    for (kk in 1:K)
      r <- r + all(hm2==hap[kk,]) + all(hc2==hap[kk,])
    if(r == 2)
      h <- c(h,k)
  }
  return(h)  #return the hap number mom gives child
}

P_Y1 <- function(Y,gm,gc,im,X,para){
  nX <- max(1,length(X))
  tmp <- exp(sum(c(1,gm,gc,im,X) * para[1:(4+nX)]))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

distinguish = function(gmm,gcc,loci,hap){
  l = nrow(gmm)
  K = nrow(hap)
  gm = gmm[,loci]
  gc = gcc[,loci]
  gcm = rep(3,l)
  gcp = rep(0,l)
  phi = rep(0,l)
  
  for(u in 1:l){
    
    if(gc[u]==0){
      gcm[u] = 0
      gcp[u] = 0
      phi[u] = 1
    }
    else if(gm[u]==0&gc[u]==1){
      gcm[u] = 0
      gcp[u] = 1
      phi[u] = 1
    }
    else if(gm[u]==1&gc[u]==1){
      for(k in 1:K){
        #mother
        hm1 = hap[k,]
        hm2 = gmm[u,] - hm1
        r = 0
        for (kk in 1:K){
          r = r + all(hm2==hap[kk,])
        }
        if(r==1){
          #pass hm1 to child
          hcm = hm1
          hcp = gcc[u,] - hcm
          t = 0
          for (kk in 1:K){
            t = t + all(hcp==hap[kk,])
          } 
          if(t==1&gcm[u]!=hcm[loci]){
            gcm[u] = hcm[loci]
            gcp[u] = hcp[loci]
            if(phi[u]==0) phi[u] = 1
            else phi[u] = phi[u] + 1
          }
        }
      }
    }
    else if(gm[u]==2&gc[u]==1){
      gcm[u] = 1
      gcp[u] = 0
      phi[u] = 1
    }
    else{
      gcm[u] = 1
      gcp[u] = 1
      phi[u] = 1
    }
  }
  return(list(gm = gm, gc = gc, gcm = gcm, gcp = gcp, phi = phi))
}

design.matrix <- function(gm,gc,gcm,gcp,X,phi){
  gcm[-which(phi==1)] <- 0
  gcp[-which(phi==1)] <- 0
  cbind(1,gm,gc,gcm-gcp,X)
}

design.matrix2 <- function(gm,gc,gcm,gcp,X,phi){
  gcm[-which(phi==1)] <- 0
  gcp[-which(phi==1)] <- 0
  cbind(1,gm,gc,gcm-gcp,X,gm * X,gc * X)
}

#ROB-HAP with NA
imprinting_robhap.na <- function(Y,gmm,gcc,X,loci,hap,f){
  #data pre-processing
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  K <- ncol(hap)
  nX <- max(1,ncol(X))
  
  group2 <- apply(is.na(cbind(gmm,gcc)),1,any)
  X1 <- as.matrix(X[!group2,])
  Y1 <- Y[!group2]
  gmm1 <- as.matrix(gmm[!group2,])
  gcc1 <- as.matrix(gcc[!group2,])
  X2 <- as.matrix(X[group2,])
  Y2 <- Y[group2]
  gmm2 <- as.matrix(gmm[group2,])
  gcc2 <- as.matrix(gcc[group2,])
  
  #primary value
  F <- distinguish(gmm1,gcc1,loci,hap)
  gm <- F$gm
  gc <- F$gc
  gcm <- F$gcm
  gcp <- F$gcp
  phi <- F$phi
  Z <- design.matrix(gm,gc,gcm,gcp,X1,phi)
  fit <- glm(Y1 ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  sd.log <- as.vector(res[,2])
  
  if(ncol(gmm1)>1){
    library(haplo.stats)
    gmm12 <- geno1to2(gmm1)
    haplo_raw <- haplo.em(gmm12)
    ppi_raw <- haplo_raw$hap.prob
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)),n_hap_raw,K,byrow = FALSE) - 1
  }
  else{
    p <- mean(gmm1)/2
    ppi_raw <- c(1-p,p)
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(c(0,1),n_hap_raw,K)
  }
  
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1]
  ind_hh <- mul_whichrow(hh_raw,hap)
  ppi_ini <- ppi_raw[ind_hh]
  ppi_ini <- ppi_ini/sum(ppi_ini)
  theta0 <- mul_logitf(ppi_ini)
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.robust.na(c(beta0,theta),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.robust(c(beta0,theta),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.robust.na(c(beta,theta0),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.robust(c(beta,theta0),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.robust.na(para,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n)-likelioriginal.robust(para,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX)],mul_logistic(est[(5+nX):((3+nX) + n_hap)])) 
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

######
#function
haplo.na <- function(hap,gmm,gcc){
  h <- NULL
  K <- nrow(hap)
  for(i in 1:K){
    for(j in 1:K){
      for(k in 1:K){
        hm1 <- hap[i,] #hap for child
        hm2 <- hap[j,]
        hc2 <- hap[k,]
        if(all((hm1 + hm2)[!is.na(gmm)] == gmm[!is.na(gmm)]) & all((hm1 + hc2)[!is.na(gcc)] == gcc[!is.na(gcc)]))
          h <- rbind(h,c(i,j,k))
      }
    }
  }
  return(h)
}

P_H.na <- function(h,gmm,gcc,hap,mu){
  return(mu[h[1]]*mu[h[2]]*mu[h[3]])
}

P_Y.na <- function(Y,X,h,gmm,gcc,hap,loci,beta){
  hm1 <- hap[h[1],]
  hm2 <- hap[h[2],]
  hc2 <- hap[h[3],]
  tmp <- exp(sum(c(1,hm1[loci]+hm2[loci],hm1[loci]+hc2[loci],hm1[loci]-hc2[loci],X) * beta))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

likeli0.robust.na <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  beta0 <- para0[1:(4+nX)]
  mu0 <- mul_logistic(para0[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + (log(P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.na(h,gmm,gcc,hap,mu))  - log(P_g.rob(gmm[loci],theta))) * P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.robust.na <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.robust.na(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.robust(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.robust.na <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.na(h,gmm,gcc,hap,mu) / P_g.rob(gmm[loci],theta)
  }
  return(u)
}

likelioriginal.robust.na <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.robust.na(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.robust(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

#IND-HAP
imprinting_indhap <- function(Y,gmm,gcc,X,loci,hap,f){
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  K <- ncol(hap)
  nX <- max(1,ncol(X))
  
  #primary value
  F <- distinguish(gmm,gcc,loci,hap)
  gm <- F$gm
  gc <- F$gc
  gcm <- F$gcm
  gcp <- F$gcp
  phi <- F$phi
  Z <- design.matrix(gm,gc,gcm,gcp,X,phi)
  fit <- glm(Y ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  sd.log <- as.vector(res[,2])
  
  if(ncol(gmm) > 1){
    library(haplo.stats)
    gmm2 <- geno1to2(gmm)
    haplo_raw <- haplo.em(gmm2)
    ppi_raw <- haplo_raw$hap.prob
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)),n_hap_raw,K,byrow = FALSE) - 1
  }
  else{
    p <- mean(gmm)/2
    ppi_raw <- c(1-p,p)
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(c(0,1),n_hap_raw,K)
  }
  
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1]
  ind_hh <- mul_whichrow(hh_raw,hap)
  ppi_ini <- ppi_raw[ind_hh]
  ppi_ini <- ppi_ini/sum(ppi_ini)
  theta0 <- mul_logitf(ppi_ini)
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.ind(c(beta0,theta),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.ind(c(beta,theta0),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.ind(para,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX)],mul_logistic(est[(5+nX):((3+nX) + n_hap)]))
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

#functions
likeli0.ind <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  beta0 <- para0[1:(4+nX)]
  mu0 <- mul_logistic(para0[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(h in hh){
    u <- u + (log(P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.rob(h,gmm,gcc,hap,mu))) * P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.ind <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.ind(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.ind(X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.ind <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  for(h in hh)
    u <- u + P_Y.rob(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.rob(h,gmm,gcc,hap,mu)
  return(u)
}

likelioriginal.ind <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.ind(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.ind(X[i,],beta,theta) - f)))
  }
  return(res)
}

H.ind <- function(X,beta,theta){
  res <- P_g.rob(0,theta) * (theta * P_Y1(1,0,1,-1,X,beta) + (1-theta) * P_Y1(1,0,0,0,X,beta))
  res <- res + 0.5*P_g.rob(1,theta) * ((theta * P_Y1(1,1,1,-1,X,beta) + (1-theta) * P_Y1(1,1,0,0,X,beta)) + (theta * P_Y1(1,1,2,0,X,beta) + (1-theta) * P_Y1(1,1,1,1,X,beta)))
  res <- res + P_g.rob(2,theta) * (theta * P_Y1(1,2,2,0,X,beta) + (1-theta) * P_Y1(1,2,1,1,X,beta))
  return(res)
}

#IND-HAP with NA
imprinting_indhap.na <- function(Y,gmm,gcc,X,loci,hap,f){
  #data pre-processing
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  K <- ncol(hap)
  nX <- max(1,ncol(X))
  
  group2 <- apply(is.na(cbind(gmm,gcc)) , 1, any)
  X1 <- as.matrix(X[!group2,])
  Y1 <- Y[!group2]
  gmm1 <- as.matrix(gmm[!group2,])
  gcc1 <- as.matrix(gcc[!group2,])
  X2 <- as.matrix(X[group2,])
  Y2 <- Y[group2]
  gmm2 <- as.matrix(gmm[group2,])
  gcc2 <- as.matrix(gcc[group2,])
  
  #primary value
  F <- distinguish(gmm1,gcc1,loci,hap)
  gm <- F$gm
  gc <- F$gc
  gcm <- F$gcm
  gcp <- F$gcp
  phi <- F$phi
  Z <- design.matrix(gm,gc,gcm,gcp,X1,phi)
  fit <- glm(Y1 ~ 0 + Z,family = binomial)
  res <- summary(fit)$coef
  est.log <- as.vector(res[,1])
  sd.log <- as.vector(res[,2])
  
  if(ncol(gmm1)>1){
    library(haplo.stats)
    gmm12 <- geno1to2(gmm1)
    haplo_raw <- haplo.em(gmm12)
    ppi_raw <- haplo_raw$hap.prob
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)),n_hap_raw,K,byrow = FALSE) - 1
  }
  else{
    p <- mean(gmm1)/2
    ppi_raw <- c(1-p,p)
    n_hap_raw <- length(ppi_raw)
    hh_raw <- matrix(c(0,1),n_hap_raw,K)
  }
  
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1]
  ind_hh <- mul_whichrow(hh_raw, hap)
  ppi_ini <- ppi_raw[ind_hh]
  ppi_ini <- ppi_ini/sum(ppi_ini)
  theta0 <- mul_logitf(ppi_ini)
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.ind.na(c(beta0,theta),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.ind(c(beta0,theta),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.ind.na(c(beta,theta0),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.ind(c(beta,theta0),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.ind.na(para,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n)-likelioriginal.ind(para,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX)],mul_logistic(est[(5+nX):((3+nX) + n_hap)])) 
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

likeli0.ind.na <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  beta0 <- para0[1:(4+nX)]
  mu0 <- mul_logistic(para0[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + (log(P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.na(h,gmm,gcc,hap,mu))) * P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.ind.na <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.ind.na(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.ind(X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.ind.na <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + P_Y.na(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.na(h,gmm,gcc,hap,mu)
    }
  return(u)
}

likelioriginal.ind.na <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX)]
  mu <- mul_logistic(para[(5+nX):((3+nX)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.ind.na(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.ind(X[i,],beta,theta) - f)))
  }
  return(res)
}

logistic <- function(x){
    return(1/(1+exp(-x)))
}

mul_logistic <- function(x){
    freq <- c(exp(x), 1)
    freq <- freq/sum(freq)
    return(freq)
}

logitf <- function(p){
    return(log(p/(1-p)))
}

mul_logitf <- function(p){
    n <- length(p)
    theta <- log(p[-n]/p[n])
    return(theta)
}

whichrow <- function(hh,h){
    ind <- which(apply(hh,1,function(x){all(x == h)}));
    if(length(ind) == 0)ind <- -1;
    return(ind);
}

mul_whichrow <- function(hh1, hh2){
    return(apply(hh2, 1, function(x){whichrow(hh1, x)}));
}

#methods with interaction
#ROB-HAP
imprinting_robhap.inter <- function(Y,gmm,gcc,X,loci,hap,f){
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  K <- ncol(hap)
  X <- as.matrix(X)
  nX <- max(1,ncol(X))
  
  #initial value
  F = distinguish(gmm, gcc, loci, hap)
  gm = F$gm
  gc = F$gc
  gcm = F$gcm
  gcp = F$gcp
  phi = F$phi
  Z = design.matrix2(gm,gc,gcm,gcp,X,phi)
  fit = glm(Y ~ 0 + Z,family = binomial)
  res = summary(fit)$coef
  est.log = as.vector(res[,1])
  sd.log = as.vector(res[,2])
  
  library(haplo.stats);
  gmm2 <- geno1to2(gmm);
  haplo_raw <- haplo.em(gmm2);
  ppi_raw <- haplo_raw$hap.prob;
  n_hap_raw <- length(ppi_raw);
  hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)), n_hap_raw, K, byrow = FALSE) - 1;
  
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1];
  ind_hh <- mul_whichrow(hh_raw, hap);
  ppi_ini <- ppi_raw[ind_hh];
  ppi_ini <- ppi_ini/sum(ppi_ini);
  theta0 <- mul_logitf(ppi_ini);
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.robust.inter(c(beta0,theta),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.robust.inter(c(beta,theta0),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.robust.inter(para,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX*3)],mul_logistic(est[(5+nX*3):(3+nX*3+n_hap)])) 
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX*3)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

####functions
likeli0.robust.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  beta0 <- para0[1:(4+nX*3)]
  mu0 <- mul_logistic(para0[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(h in hh){
    u <- u + (log(P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.rob(h,gmm,gcc,hap,mu))  - log(P_g.rob(gmm[loci],theta))) * P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.robust.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.robust.inter(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.inter(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.robust.inter <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  for(h in hh)
    u <- u + P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.rob(h,gmm,gcc,hap,mu) / P_g.rob(gmm[loci],theta)
  return(u)
}

likelioriginal.robust.inter <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.robust.inter(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.inter(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

P_Y.inter <- function(Y,X,h,gmm,gcc,hap,loci,beta){
  hm1 <- hap[h,]
  tmp <- exp(sum(c(1,gmm[loci],gcc[loci],2*hm1[loci]-gcc[loci],X,gmm[loci]*X,gcc[loci]*X) * beta))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

P_Y1.inter <- function(Y,gm,gc,im,X,para){
  nX <- max(1,length(X))
  tmp <- exp(sum(c(1,gm,gc,im,X,gm*X,gc*X) * para[1:(4+nX*3)]))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

H.inter <- function(gm,X,beta,theta){
  if(gm == 0)
    res <- P_Y1.inter(1,0,0,0,X,beta) * (1-theta) + P_Y1.inter(1,0,1,-1,X,beta) * theta
  else if(gm == 1)
    res <- 0.5* (P_Y1.inter(1,1,0,0,X,beta) * (1-theta) + P_Y1.inter(1,1,1,-1,X,beta) * theta + P_Y1.inter(1,1,1,1,X,beta) * (1-theta) + P_Y1.inter(1,1,2,0,X,beta) * theta)
  else if(gm == 2)
    res <- P_Y1.inter(1,2,1,1,X,beta) * (1-theta) + P_Y1.inter(1,2,2,0,X,beta) * theta
  return(res)
}

#ROB-HAP-EM
imprinting_robhap.na.inter <- function(Y,gmm,gcc,X,loci,hap,f){
  #data pre-processing
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  K <- ncol(hap)
  X <- as.matrix(X)
  nX <- max(1,ncol(X))
  
  group2 <- apply(is.na(cbind(gmm,gcc)) , 1, any)
  X1 <- as.matrix(X[!group2,])
  Y1 <- Y[!group2]
  gmm1 <- as.matrix(gmm[!group2,])
  gcc1 <- as.matrix(gcc[!group2,])
  X2 <- as.matrix(X[group2,])
  Y2 <- Y[group2]
  gmm2 <- as.matrix(gmm[group2,])
  gcc2 <- as.matrix(gcc[group2,])
  
  #primary value
  F = distinguish(gmm1, gcc1, loci, hap)
  gm = F$gm
  gc = F$gc
  gcm = F$gcm
  gcp = F$gcp
  phi = F$phi
  Z = design.matrix2(gm,gc,gcm,gcp,X1,phi)
  fit = glm(Y1 ~ 0 + Z,family = binomial)
  res = summary(fit)$coef
  est.log = as.vector(res[,1])
  sd.log = as.vector(res[,2])
  
  if(ncol(gmm1)>1)
  {
  library(haplo.stats);
  gmm12 <- geno1to2(gmm1);
  haplo_raw <- haplo.em(gmm12);
  ppi_raw <- haplo_raw$hap.prob;
  n_hap_raw <- length(ppi_raw);
  hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)), n_hap_raw, K, byrow = FALSE) - 1;
  }
  else
  {
    p = mean(gmm1)/2
    ppi_raw <- c(1-p,p)
    n_hap_raw <- length(ppi_raw);
    hh_raw = matrix(c(0,1), n_hap_raw, K)
  }
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1];
  ind_hh <- mul_whichrow(hh_raw, hap);
  ppi_ini <- ppi_raw[ind_hh];
  ppi_ini <- ppi_ini/sum(ppi_ini);
  theta0 <- mul_logitf(ppi_ini);
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.robust.na.inter(c(beta0,theta),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.robust.inter(c(beta0,theta),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.robust.na.inter(c(beta,theta0),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.robust.inter(c(beta,theta0),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.robust.na.inter(para,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n)-likelioriginal.robust.inter(para,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX*3)],mul_logistic(est[(5+nX*3):((3+nX*3) + n_hap)])) 
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX*3)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

likeli0.robust.na.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  beta0 <- para0[1:(4+nX*3)]
  mu0 <- mul_logistic(para0[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + (log(P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.na(h,gmm,gcc,hap,mu))  - log(P_g.rob(gmm[loci],theta))) * P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.robust.na.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.robust.na.inter(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.inter(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.robust.na.inter <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.na(h,gmm,gcc,hap,mu) / P_g.rob(gmm[loci],theta)
  }
  return(u)
}

likelioriginal.robust.na.inter <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.robust.na.inter(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.inter(gmm[i,loci],X[i,],beta,theta) - f)))
  }
  return(res)
}

P_Y.na.inter <- function(Y,X,h,gmm,gcc,hap,loci,beta){
  hm1 <- hap[h[1],]
  hm2 <- hap[h[2],]
  hc2 <- hap[h[3],]
  tmp <- exp(sum(c(1,hm1[loci]+hm2[loci],hm1[loci]+hc2[loci],hm1[loci]-hc2[loci],X,(hm1[loci]+hm2[loci])*X,(hm1[loci]+hc2[loci])*X) * beta))
  return(1-Y+(2*Y-1)*tmp/(1+tmp))
}

#IND-HAP
imprinting_indhap.inter <- function(Y,gmm,gcc,X,loci,hap,f){
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  X <- as.matrix(X)
  K <- ncol(hap)
  nX <- max(1,ncol(X))
  
  #primary value
  F = distinguish(gmm, gcc, loci, hap)
  gm = F$gm
  gc = F$gc
  gcm = F$gcm
  gcp = F$gcp
  phi = F$phi
  Z = design.matrix2(gm,gc,gcm,gcp,X,phi)
  fit = glm(Y ~ 0 + Z,family = binomial)
  res = summary(fit)$coef
  est.log = as.vector(res[,1])
  sd.log = as.vector(res[,2])
  
  library(haplo.stats);
  gmm2 <- geno1to2(gmm);
  haplo_raw <- haplo.em(gmm2);
  ppi_raw <- haplo_raw$hap.prob;
  n_hap_raw <- length(ppi_raw);
  hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)), n_hap_raw, K, byrow = FALSE) - 1;
  
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1];
  ind_hh <- mul_whichrow(hh_raw, hap);
  ppi_ini <- ppi_raw[ind_hh];
  ppi_ini <- ppi_ini/sum(ppi_ini);
  theta0 <- mul_logitf(ppi_ini);
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.ind.inter(c(beta0,theta),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.ind.inter(c(beta,theta0),para0,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.ind.inter(para,Y,X,gmm,gcc,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX*3)],mul_logistic(est[(5+nX*3):((3+nX*3) + n_hap)]))
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX*3)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

likeli0.ind.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  beta0 <- para0[1:(4+nX*3)]
  mu0 <- mul_logistic(para0[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(h in hh){
    u <- u + (log(P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.rob(h,gmm,gcc,hap,mu))) * P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.rob(h,gmm,gcc,hap,mu0) * P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.ind.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.ind.inter(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.ind.inter(X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.ind.inter <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo(hap,gmm,gcc)
  u <- 0
  for(h in hh)
    u <- u + P_Y.inter(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.rob(h,gmm,gcc,hap,mu)
  return(u)
}

likelioriginal.ind.inter <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.ind.inter(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.ind.inter(X[i,],beta,theta) - f)))
  }
  return(res)
}

H.ind.inter <- function(X,beta,theta){
  res <- P_g.rob(0,theta) * (theta * P_Y1.inter(1,0,1,-1,X,beta) + (1-theta) * P_Y1.inter(1,0,0,0,X,beta))
  res <- res + 0.5*P_g.rob(1,theta) * ((theta * P_Y1.inter(1,1,1,-1,X,beta) + (1-theta) * P_Y1.inter(1,1,0,0,X,beta)) + (theta * P_Y1.inter(1,1,2,0,X,beta) + (1-theta) * P_Y1.inter(1,1,1,1,X,beta)))
  res <- res + P_g.rob(2,theta) * (theta * P_Y1.inter(1,2,2,0,X,beta) + (1-theta) * P_Y1.inter(1,2,1,1,X,beta))
  return(res)
}

#IND-HAP-EM
imprinting_indhap.na.inter <- function(Y,gmm,gcc,X,loci,hap,f){
  #data pre-processing
  n1 <- sum(Y == 1)
  n0 <- sum(Y == 0)
  n <- n1+n0
  lambda <- n1 / (n*f) - n0 / (n * (1-f))
  K <- ncol(hap)
  X <- as.matrix(X)
  nX <- max(1,ncol(X))
  
  group2 <- apply(is.na(cbind(gmm,gcc)) , 1, any)
  X1 <- as.matrix(X[!group2,])
  Y1 <- Y[!group2]
  gmm1 <- gmm[!group2,]
  gcc1 <- gcc[!group2,]
  X2 <- as.matrix(X[group2,])
  Y2 <- Y[group2]
  gmm2 <- gmm[group2,]
  gcc2 <- gcc[group2,]
  
  #primary value
  F = distinguish(gmm1, gcc1, loci, hap)
  gm = F$gm
  gc = F$gc
  gcm = F$gcm
  gcp = F$gcp
  phi = F$phi
  Z = design.matrix2(gm,gc,gcm,gcp,X1,phi)
  fit = glm(Y1 ~ 0 + Z,family = binomial)
  res = summary(fit)$coef
  est.log = as.vector(res[,1])
  sd.log = as.vector(res[,2])
  
  library(haplo.stats);
  gmm12 <- geno1to2(gmm1);
  haplo_raw <- haplo.em(gmm12);
  ppi_raw <- haplo_raw$hap.prob;
  n_hap_raw <- length(ppi_raw);
  hh_raw <- matrix(as.numeric(as.matrix(haplo_raw$haplotype)), n_hap_raw, K, byrow = FALSE) - 1;
  
  # assuming true haplotypes in population
  n_hap <- dim(hap)[1];
  ind_hh <- mul_whichrow(hh_raw, hap);
  ppi_ini <- ppi_raw[ind_hh];
  ppi_ini <- ppi_ini/sum(ppi_ini);
  theta0 <- mul_logitf(ppi_ini);
  beta0 <- est.log
  para0 <- c(beta0,theta0)
  
  #EM
  iteration <- 5
  llik1 <- function(theta)
    -likeli.ind.na.inter(c(beta0,theta),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.ind.inter(c(beta0,theta),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  llik2 <- function(beta)
    -likeli.ind.na.inter(c(beta,theta0),para0,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n) -likeli.ind.inter(c(beta,theta0),para0,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  for(r in 1:iteration){
    fit1 <- optim(par = theta0,fn = llik1,method = 'L-BFGS-B')
    theta0 <- fit1$par
    para0 <- c(beta0,theta0)
    fit2 <- optim(par = beta0,fn = llik2,method = 'L-BFGS-B')
    beta0 <- fit2$par
    para0 <- c(beta0,theta0)
  }
  llik <- function(para)
    -likelioriginal.ind.na.inter(para,Y2,X2,gmm2,gcc2,hap,loci,f,lambda,n)-likelioriginal.ind.inter(para,Y1,X1,gmm1,gcc1,hap,loci,f,lambda,n)
  fit3 <- optim(par = para0,fn = llik,method = 'L-BFGS-B',hessian = TRUE)
  est <- fit3$par
  est <- c(est[1:(4+nX*3)],mul_logistic(est[(5+nX*3):((3+nX*3) + n_hap)])) 
  sd <- sqrt(diag(solve(fit3$hessian)))[1:(4+nX*3)]
  return(list(est = est,sd = sd,est.log = est.log,sd.log = sd.log))
}

likeli0.ind.na.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci){ #likelihood for one person
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  d <- 0
  K <- nrow(hap)
  nX <- max(1,length(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  beta0 <- para0[1:(4+nX*3)]
  mu0 <- mul_logistic(para0[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + (log(P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta)) + log(P_H.na(h,gmm,gcc,hap,mu))) * P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
    d <- d + P_H.na(h,gmm,gcc,hap,mu0) * P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta0)
  }
  return(u/d)
}

likeli.ind.na.inter <- function(para,para0,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + likeli0.ind.na.inter(para,para0,Y[i],X[i,],gmm[i,],gcc[i,],hap,loci)
    res <- res - log(n * (1 + lambda * (H.ind.inter(X[i,],beta,theta) - f)))
  }
  return(res)
}

likelioriginal0.ind.na.inter <- function(Y,X,gmm,gcc,hap,loci,beta,mu){
  theta <- sum(hap[,loci] * mu)
  hh <- haplo.na(hap,gmm,gcc)
  u <- 0
  for(i in 1:nrow(hh)){
    h <- hh[i,]
    u <- u + P_Y.na.inter(Y,X,h,gmm,gcc,hap,loci,beta) * P_H.na(h,gmm,gcc,hap,mu)
  }
  return(u)
}

likelioriginal.ind.na.inter <- function(para,Y,X,gmm,gcc,hap,loci,f,lambda,n){
  res <- 0
  K <- nrow(hap)
  nX <- max(1,ncol(X))
  beta <- para[1:(4+nX*3)]
  mu <- mul_logistic(para[(5+nX*3):((3+nX*3)+K)])
  theta <- sum(hap[,loci] * mu)
  for(i in 1:nrow(gmm)){
    res <- res + log(likelioriginal0.ind.na.inter(Y[i],X[i,],gmm[i,],gcc[i,],hap,loci,beta,mu))
    res <- res - log(n * (1 + lambda * (H.ind.inter(X[i,],beta,theta) - f)))
  }
  return(res)
}
