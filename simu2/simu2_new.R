

library(mvtnorm)
ff.Wei.reg = function(n,theta,nmu,sigma,b0,b1,b2,t,B=500,type=1,fix=0,seed=1){ 
  #nmu = multivariate mean ; sigma = multivariate covariance
  # fix=1: known shape ; 0 = unknown
  logL.c = function(para,x,xx,type,fix){
    # if fix = 0, para = c(theta,beta0,beta1*I(type=1),beta2*I(type(1,2)))
    if (fix == 0){
      if(type ==1){ mu = para[2] + colSums( t(xx) * c(para[3],para[4]) ) }
      else if(type ==2) { mu=para[2] + xx[,1] * para[3] }
      else if(type ==3) { mu=para[2] + xx[,2] * para[3] }
      else if(type ==4) { mu=rep( para[2] , n) }
    }
    else if (fix == 1){
      if(type ==1){ mu = para[1] + colSums( t(xx) * c(para[2],para[3]) ) }
      else if(type ==2) { mu=para[1] + xx[,1] * para[2] }
      else if(type ==3) { mu=para[1] + xx[,2] * para[2] }
      else if(type ==4) { mu=rep(x = para,n) }
    }
    di= as.numeric(x <=t); r = sum(di)
    sig=ifelse(fix == 1,yes = 1/theta, no = 1/para[1])
    if(sig<0){sig=0.001}
    l=0
    z=(log(x)-mu)/sig
    zt = (log(t)-mu)/sig
    for(i in 1:length(x)){
      l=l+di[i]*(-log(sig)- log(x[i])+z[i]-exp(z[i]) ) -
        (1-di[i])*exp((zt[i]))
    }
    return(l)
  }
  logL.t=function(para,x,xx,type,fix){
    # if fix = 0, para = c(theta,beta0,beta1*I(type=1),beta2*I(type(1,2)))
    if (fix == 0){
      if(type ==1){ mu = para[2] + colSums( t(xx) * c(para[3],para[4]) ) }
      else if(type ==2) { mu=para[2] + xx[,1] * para[3] }
      else if(type ==3) { mu=para[2] + xx[,2] * para[3] }
      else if(type ==4) { mu=rep( para[2] , n) }
    }
    else if (fix == 1){
      if(type ==1){ mu = para[1] + colSums( t(xx) * c(para[2],para[3]) ) }
      else if(type ==2) { mu=para[1] + xx[,1] * para[2] }
      else if(type ==3) { mu=para[1] + xx[,2] * para[2] }
      else if(type ==4) { mu=rep(para,n) }
    }
    di= as.numeric(x <=t); r = sum(di)
    ifelse(fix == 1, sig<- 1/theta, sig <- 1/para[1])
    if(sig<0){sig=0.001}
    l=0
    y = x[di==1]
    mean = mu[di==1]
    # for (i in 1:length(y)){
    #   l = l + log(dexp(x = y[i], rate = 1/exp(mean[i]))) - log(pexp(q = t, rate = 1/exp(mean[i])))
    # }
    z=(log(y)-mean)/sig
    zt = (log(t)-mean)/sig
    for(i in 1:length(y)){
      #l = l+di[i] * dexp(x = x[i], rate = 1/exp(mu[i]) )/ pexp(q = t, rate = 1/exp(mu[i]))
      l=l+(-log(sig)- log(y[i])+z[i]-exp(z[i]) - log(1-exp(-exp(zt[i]))) )
      #l=l+(-log(sig)- log(y[i])+z[i]-exp(z[i]) )
    }
    return(l)
  } # truncated Weibull log L
  
  #  sigma.matrix = matrix(c(0.1,p*0.1,p*0.1,0.1),2,2)
  set.seed(seed) ; bigX = rmvnorm( n, mean = nmu, sigma = sigma )
  r.mu = exp( b0 + bigX %*% c(b1,b2) )
  true.mu = exp( b0 + sum( nmu * c(b1, b2) ) + 1/2* ( c(b1,b2) %*% sigma %*% c(b1,b2)) )
  ifelse( type== 4, cutoff <- 10, cutoff <- 2.5 )
  simu = replicate( B, expr={
    x=0
    for(i in 1:n){ x[i] = rexp(1,rate = 1/r.mu[i]^theta)^(1/theta) }# data
    di= as.numeric(x <=t); r = sum(di)
    if( type == 1 ) { ifelse(fix==1, para <- c(1,1,1), para <- c(theta,1,1,1) ) }
    else if ( type ==2 || type ==3) { 
      ifelse(fix==1, para <- c(1,1), para <- c(theta,1,1) ) 
    }
    else if ( type == 4 ){ ifelse(fix==1, para <- b0, para <- c(theta,1) ) }
    
    ifelse ( type == 4 && fix == 1 ,
             op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                              control=list("fnscale"=-1),hessian=T,
                              method =  "Brent" , lower = 0.5, upper = 5) ,
             op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                              control=list("fnscale"=-1),hessian=T) )
    
    ifelse( fix==1,mle.sha.t <- theta, mle.sha.t <- op.trun$par[1] )
    ifelse( fix==1,beta.hat.t <- op.trun$par, beta.hat.t <- op.trun$par[-1] )
    
    if (type == 1 ) {  adj.X = cbind(1, bigX) }
    else if (type == 2) {  adj.X = cbind(1, bigX[,1]) }
    else if (type == 3) {  adj.X = cbind(1, bigX[,2]) }
    else if (type == 4) {  adj.X = 1 }
    #mle.sca.t = sum( exp( adj.X %*% beta.hat.t)* di )/ r 
    mle.sca.t = mean( exp( adj.X %*% beta.hat.t) )
    #mle.sca.t = exp( beta.hat.t[1] + 1/2* (beta.hat.t[-1] %*% sigma %*% beta.hat.t[-1]) )
    # sample.mean = sum( exp( adj.X %*% beta.hat.t)* di )/ r 
    # ifelse(sample.mean <= t/2, mle.sca.t <- sample.mean, mle.sca.t <- t/2) 
    ratio = mle.sca.t/true.mu * gamma(1+1/mle.sha.t)/gamma(1+1/theta)# determine if discarding?
    while(ratio >=cutoff  || ratio <= 1/cutoff || op.trun$convergence != 0){
      x=0
      for(i in 1:n){x[i] = rexp(1,rate = 1/r.mu[i]^theta)^(1/theta)}# data
      di= as.numeric(x <=t); r = sum(di)
      if( type == 1 ) { ifelse(fix==1, para <- c(1,1,1), para <- c(theta,1,1,1) ) }
      else if ( type ==2 || type ==3) { 
        ifelse(fix==1,para <- c(1,1), para <- c(theta,1,1) ) 
      }
      else if ( type == 4 ){ ifelse(fix==1,para <- b0, para <- c(theta,1) ) }
      
      ifelse ( type == 4 && fix == 1 ,
               op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                                control=list("fnscale"=-1),hessian=T,
                                method =  "Brent" , lower = 0.5, upper = 5) ,
               op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                                control=list("fnscale"=-1),hessian=T) )
      
      ifelse( fix==1,mle.sha.t <- theta, mle.sha.t <- op.trun$par[1] )
      ifelse( fix==1,beta.hat.t <- op.trun$par,beta.hat.t <- op.trun$par[-1] )
      # adj.X = cbind(1,bigX) ;
      if (type == 1 ) {  adj.X = cbind(1, bigX) }
      else if (type == 2) {  adj.X = cbind(1, bigX[,1]) }
      else if (type == 3) {  adj.X = cbind(1, bigX[,2]) }
      else if (type == 4) {  adj.X = 1 }
      mle.sca.t = mean( exp( adj.X %*% beta.hat.t) )
      #mle.sca.t = exp( beta.hat.t[1] + 1/2* (beta.hat.t[-1] %*% sigma %*% beta.hat.t[-1]) )
      # sample.mean = sum( exp( adj.X %*% beta.hat.t)* di )/ r 
      # ifelse(sample.mean <= t/2, mle.sca.t <- sample.mean, mle.sca.t <- t/2) 
      ratio = mle.sca.t/true.mu * gamma(1+1/mle.sha.t)/gamma(1+1/theta)# determine if discarding?
    }
    ifelse ( type == 4 && fix == 1 ,
             op.cen <- optim(para,logL.c,x=x,xx=bigX,type=type,fix=fix,
                             control=list("fnscale"=-1),hessian=T,
                             method =  "Brent" , lower = 0.5, upper = 5) ,
             op.cen <- optim(para,logL.c,x=x,xx=bigX,type=type,fix=fix,
                             control=list("fnscale"=-1),hessian=T) )
    ifelse( fix==1,mle.sha.c <-theta, mle.sha.c <- op.cen$par[1] )
    ifelse( fix==1,beta.hat.c <- op.cen$par,beta.hat.c <- op.cen$par[-1] )
    
    mle.sca.c = mean( exp( adj.X %*% beta.hat.c ) )  
    
    ## cp(Wald for mu only) 
    
    if (fix == 0){
      Fisher.c = solve(-op.cen$hessian)
      Fisher.t = solve(-op.trun$hessian)
      if ( type == 1 ){
        c.var.mu = apply(X = adj.X, MARGIN = 1, FUN = function(x){
          exp(x %*% beta.hat.c)^2 * (x %*% Fisher.c %*% x)
        })
        t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
          exp(x[-4] %*% beta.hat.t)^2 * (x[-4] %*% Fisher.t %*% x[-4]) * x[4]
        })
        c.CI = mle.sca.c + c(-1,1) * qnorm(0.975) * sqrt(sum(c.var.mu))/n
        t.CI = mle.sca.t + c(-1,1) * qnorm(0.975) * sqrt(sum(t.var.mu))/r
      }
      else if ( type == 2 || type == 3){
        c.var.mu = apply(X = adj.X, MARGIN = 1, FUN = function(x){
          exp(x %*% beta.hat.c)^2 * (x %*% Fisher.c %*% x)
        })
        t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
          exp(x[-3] %*% beta.hat.t)^2 * (x[-3] %*% Fisher.t %*% x[-3]) * x[3]
        })
        c.CI = mle.sca.c + c(-1,1) * qnorm(0.975) * sqrt(sum(c.var.mu))/n
        t.CI = mle.sca.t + c(-1,1) * qnorm(0.975) * sqrt(sum(t.var.mu))/r
      }  
      else if ( type == 4){
        c.var.mu = exp(beta.hat.c)^2 * Fisher.c
        t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
          exp(x[-2] * beta.hat.t)^2 * (x[-2]^2 * Fisher.t) * x[2]
        })
        c.CI = mle.sca.c + c(-1,1) * qnorm(0.975) * sqrt(c.var.mu/n)
        t.CI = mle.sca.t + c(-1,1) * qnorm(0.975) * sqrt(sum(t.var.mu))/r
      }
      sigma.c.CI =  mle.sha.c + c(-1,1) * qnorm(0.975) * sqrt(Fisher.c[2,2])
      sigma.t.CI =  mle.sha.t + c(-1,1) * qnorm(0.975) * sqrt(Fisher.t[2,2])
    }
    
    else if( fix == 1 ){
      Fisher.c = solve(-op.cen$hessian)
      Fisher.t = solve(-op.trun$hessian)
      if ( type == 1 ){
        c.var.mu = apply(X = adj.X, MARGIN = 1, FUN = function(x){
          exp(2 * x %*% beta.hat.c) * (x %*% Fisher.c %*% x) 
        })
        t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
          #exp(2 * x[-4] %*% beta.hat.t) * (x[-4] %*% Fisher.t %*% x[-4]) *x[4]
          exp(2 * x[-4] %*% beta.hat.t) * (x[-4] %*% Fisher.t %*% x[-4])
        })
        sd.cen = sqrt( mean(c.var.mu) ) ; sd.trun = sqrt( mean(t.var.mu) )
        c.CI.l =  mean( exp(adj.X %*% beta.hat.c) ) - qnorm(0.975) * sd.cen
        c.CI.u =  mean( exp(adj.X %*% beta.hat.c) ) + qnorm(0.975) * sd.cen
        t.CI.l =  mle.sca.t - qnorm(0.975) * sd.trun
        t.CI.u =  mle.sca.t + qnorm(0.975) * sd.trun
      }
      else if ( type == 2 || type == 3){
        c.var.mu = apply(X = adj.X, MARGIN = 1, FUN = function(x){
          exp(x %*% beta.hat.c)^2 * (x %*% Fisher.c %*% x)
        })
        t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
          exp(x[-3] %*% beta.hat.t)^2 * (x[-3] %*% Fisher.t %*% x[-3])
          #exp(x[-3] %*% beta.hat.t)^2 * (x[-3] %*% Fisher.t %*% x[-3]) *x[3]
        })
        sd.cen = sqrt( mean(c.var.mu) ) ; sd.trun = sqrt( mean(t.var.mu) )
        c.CI.l =  mean( exp(adj.X %*% beta.hat.c) ) - qnorm(0.975) * sd.cen
        c.CI.u =  mean( exp(adj.X %*% beta.hat.c) ) + qnorm(0.975) * sd.cen
        t.CI.l =  mle.sca.t - qnorm(0.975) * sd.trun
        t.CI.u =  mle.sca.t + qnorm(0.975) * sd.trun
      }  
      else if ( type == 4){
        c.var.mu = exp(beta.hat.c)^2 * Fisher.c
        t.var.mu = exp(beta.hat.t)^2 * Fisher.t
        # t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
        #   exp(x[-2] * beta.hat.t)^2 * (x[-2]^2 * Fisher.t) * x[2]
        # })
        sd.cen = sqrt( c.var.mu) ; sd.trun = sqrt(t.var.mu) 
        c.CI.l = exp(adj.X %*% beta.hat.c) - qnorm(0.975) * sqrt(c.var.mu) 
        c.CI.u = exp(adj.X %*% beta.hat.c) + qnorm(0.975) * sqrt(c.var.mu) 
        t.CI.l = mle.sca.t - qnorm(0.975) * sqrt(t.var.mu)
        t.CI.u = mle.sca.t + qnorm(0.975) * sqrt(t.var.mu)
      }
    }
    
    cp.c = as.numeric( true.mu >= c.CI.l && true.mu <= c.CI.u )
    cp.t = as.numeric( true.mu >= t.CI.l && true.mu <= t.CI.u )
    
    
    
    if ( type == 1 ){
      return(c("cen.b0"=beta.hat.c[1],"cen.b1"=beta.hat.c[2],"cen.b2"=beta.hat.c[3],
               "cen.mu" =mle.sca.c,
               "cen.theta"=mle.sha.c,
               "trun.b0"=beta.hat.t[1],"trun.b1"=beta.hat.t[2],"trun.b2"=beta.hat.t[3],
               "trun.mu" =mle.sca.t,
               "trun.theta"=mle.sha.t,
               "count"=r,
               "cp.censored" = cp.c,
               "cp.truncated" = cp.t, 
               "cen.sd" = sd.cen, "trun.sd"= sd.trun) )
    }
    else if ( type == 2 || type ==3 ){
      return(c("cen.b0"=beta.hat.c[1],"cen.b1"=beta.hat.c[2],
               "cen.mu" =mle.sca.c,
               "cen.theta"=mle.sha.c,
               "trun.b0"=beta.hat.t[1],"trun.b1"=beta.hat.t[2],
               "trun.mu" =mle.sca.t,
               "trun.theta"=mle.sha.t,
               "count"=r,
               "cp.censored" = cp.c,
               "cp.truncated" = cp.t,
               "cen.sd" = sd.cen, "trun.sd"= sd.trun) )
    }
    else if (type == 4){
      return(c("cen.b0"=beta.hat.c[1],
               "cen.mu" =mle.sca.c,
               "cen.theta"=mle.sha.c,
               "trun.b0"=beta.hat.t[1],
               "trun.mu" =mle.sca.t,
               "trun.theta"=mle.sha.t,
               "count"=r,
               "cp.censored" = cp.c,
               "cp.truncated" = cp.t,
               "cen.sd" = sd.cen, "trun.sd"= sd.trun) )
    }
    
  })
  
  read.mu = exp( cbind(1, bigX) %*% c(b0,b1,b2) )
  if ( type == 1 ){
    est.c.b0 = mean(simu[1,])
    est.c.b1 = mean(simu[2,])
    est.c.b2 = mean(simu[3,])
    est.c.beta = c(est.c.b0,est.c.b1,est.c.b2)
    
    #est.c.mu = exp( rowMeans( cbind(1,bigX) %*% simu[1:3,] ) )  
    #est.c.mu = mean( exp( cbind(1, bigX) %*% est.c.beta))
    est.c.mu = exp(est.c.beta[1] + 1/2* (est.c.beta[-1] %*%sigma %*% est.c.beta[-1]) )
    #est.c.mu = mean( simu[4,] ) 
    rmse.c.mu = sd(simu[4,])
    a.c.sd = mean(simu[14,] )
    est.c.theta = mean( simu[5,] )
    mle.rate.censored = mean( exp(-(t/est.c.mu)^est.c.theta) )
    
    est.t.b0 = mean(simu[6,])
    est.t.b1 = mean(simu[7,])
    est.t.b2 = mean(simu[8,])
    est.t.beta = c(est.t.b0,est.t.b1,est.t.b2)
    
    #est.t.mu =  exp( rowMeans( cbind(1,bigX) %*% simu[6:8,] ) ) 
    #est.t.mu = mean( simu[9,] ) 
    #est.t.mu = mean( exp( cbind(1, bigX) %*% est.t.beta) )
    est.t.mu = exp(est.t.beta[1] + 1/2* (est.t.beta[-1] %*% sigma %*% est.t.beta[-1]) )
    rmse.t.mu = sd( simu[9,] )
    a.t.sd = mean(simu[15,] )
    est.t.theta = mean(simu[10,])
    mle.rate.truncated =  exp(-(t/est.t.mu)^est.t.theta) 
    
    mean.bias.c = est.c.mu*gamma(1+1/est.c.theta) - true.mu * gamma(1+1/theta) 
    mean.bias.t = est.t.mu*gamma(1+1/est.t.theta) - true.mu * gamma(1+1/theta)  
    mean.rate = 1-mean(simu[11,]/n)
    
    WaldCP.c = mean( simu[12,] )
    WaldCP.t = mean( simu[13,] )
  }
  
  else if (type ==2 || type ==3) {
    est.c.b0 = mean(simu[1,])
    est.c.b1 = mean(simu[2,])
    est.c.beta = c(est.c.b0,est.c.b1)
    if (type == 2) { adjX = bigX[,1]}
    else { adjX = bigX[,2] }
   # est.c.mu = exp( rowMeans( cbind(1,adjX) %*% simu[1:2,] ) ) 
    
    est.c.mu = mean( simu[3,] ) 
    #est.c.mu = mean( exp( cbind(1, adjX) %*% est.c.beta) )
    rmse.c.mu = sd( simu[3,] )
    a.c.sd = mean(simu[12,] )
    est.c.theta = mean( simu[4,] )
    mle.rate.censored =  exp(-(t/est.c.mu)^est.c.theta)
    
    est.t.b0 = median(simu[5,])
    est.t.b1 = median(simu[6,])
    est.t.beta = c(est.t.b0,est.t.b1)
    
   # est.t.mu =  exp( rowMeans( cbind(1,adjX) %*% simu[5:6,] ) ) 
    est.t.mu = median( simu[7,])
    #est.t.mu = mean( exp( cbind(1, adjX) %*% est.t.beta) )
    rmse.t.mu = sd(simu[7,])
    a.t.sd = mean(simu[13,] )
    est.t.theta = mean( simu[8,] )
    mle.rate.truncated =  exp(-(t/est.t.mu)^est.t.theta)
    
    mean.bias.c = est.c.mu*gamma(1+1/est.c.theta) - true.mu * gamma(1+1/theta) 
    mean.bias.t = est.t.mu*gamma(1+1/est.t.theta) - true.mu * gamma(1+1/theta) 
    mean.rate = 1-mean(simu[9,]/n)
    
    WaldCP.c = mean( simu[10,] )
    WaldCP.t = mean( simu[11,] )
  }
  else if (type ==4 ) {
    est.c.b0 = mean(simu[1,])
    est.c.beta = est.c.b0
    
   # est.c.mu =  exp( est.c.b0 )
    est.c.mu = mean( simu[2, ] ) ; rmse.c.mu = sqrt( mean((simu[2,] - est.c.mu)^2) )
    a.c.sd = mean(simu[10,] )
    est.c.theta = mean(simu[3,])
    mle.rate.censored = mean( exp(-(t/est.c.mu )^est.c.theta) )
    
    est.t.b0 = mean(simu[4,])
    est.t.beta = est.t.b0
    
   # est.t.mu =  exp( est.t.b0 )
    est.t.mu = mean( simu[5,] ) ; rmse.t.mu = sqrt( mean((simu[5,] - est.t.mu)^2) )
    a.t.sd = mean(simu[11,] )
    est.t.theta = mean( simu[6,] )
    mle.rate.truncated = mean( exp(-(t/est.t.mu )^est.t.theta ) )
    
    mean.bias.c = est.c.mu*gamma(1+1/est.c.theta) - true.mu * gamma(1+1/theta) 
    mean.bias.t = est.t.mu*gamma(1+1/est.t.theta) - true.mu * gamma(1+1/theta) 
    mean.rate = 1-mean(simu[7,]/n)
    
    WaldCP.c = mean( simu[8,] )
    WaldCP.t = mean( simu[9,] )
  }
  return(list("bigX" = bigX, "true.mu" = true.mu,
              "Censored.para"= c(est.c.beta, est.c.theta),
              "Truncated.para"= c(est.t.beta, est.t.theta),
              "Censored.sesd" = c(rmse.c.mu, a.c.sd),
              "Truncated.sesd" = c(rmse.t.mu, a.t.sd),
              "Censored.mean.bias"= c(mean.bias.c),
              "Truncated.mean.bias"= c(mean.bias.t),
              "est.rate.censored" = mle.rate.censored,
              "est.rate.truncated" = mle.rate.truncated,
              "mean.rate" = mean.rate,
              "CP" = c(WaldCP.c, WaldCP.t)))
}


test = ff.Wei.reg(n = 1000,theta = 1,nmu = c(0,0), 
                  sigma = matrix(c(1, 0.5, 0.5,1), 2,2) ,b0 = 1, b1 = 0.8,b2 = 0.8,t = 3,fix = 1, type = 2)
#######################
t = c(3, 5, 10)
beta1 = c(0.5, 0.9)
rho = c(-0.9, -0.7, -.5, -0.2 ,0, 0.2, 0.5, 0.7, .9)

combind=expand.grid(Time = t,beta1=beta1, rho=rho )

generate.data = function(combind,type){
  comb = as.matrix(combind)
  sigma = matrix(c(1, comb[3], comb[3],1), 2,2)
  true.mu = exp( 1 + 1/2* ( c(comb[2], comb[2]) %*% sigma %*% c(comb[2], comb[2])) )
  
  result = ff.Wei.reg(n = 1000,theta = 1,b0=1,b1=comb[2],b2=comb[2],t=comb[1],
                      nmu = c(0,0), sigma = sigma, type = type, fix=1)
  c.bias = result$Censored.mean.bias
  c.rmse = result$Censored.sesd[1]; c.ase = result$Censored.sesd[2]
  t.bias = result$Truncated.mean.bias
  t.rmse = result$Truncated.sesd[1]; t.ase = result$Truncated.sesd[2]
  c.rate = result$est.rate.censored
  t.rate = result$est.rate.truncated
  m.rate = result$mean.rate 
  true.rate = exp(-comb[1]/true.mu)
  CP.c = result$CP[1] ; CP.t = result$CP[2]
  
  final1 = data.frame( "Mean.bias" = c.bias, "true.mean"=true.mu,
                       "Mean.rmse"=c.rmse,"Mean.ase"=c.ase,
                       "Rate" = c.rate, "m.rate" = m.rate, "true.rate" = true.rate, 
                       "Time"= comb[1], "rho" = comb[3],
                       "beta1" = comb[2], "CP" = CP.c, "Scenario" = "censored")
  final2 = data.frame( "Mean.bias" = t.bias, "true.mean"=true.mu,
                       "Mean.rmse"=t.rmse,"Mean.ase"=t.ase,
                       "Rate" =t.rate, "m.rate" = m.rate, "true.rate" = true.rate, 
                       "Time"= comb[1], "rho" = comb[3],
                       "beta1" = comb[2], "CP" = CP.t, "Scenario" = "truncated")
  final = rbind( final1, final2)
  return(final)
}

library(parallel)  # do paraelle computing
myCoreNums <- detectCores()
cl <- makeCluster(myCoreNums-1)
clusterExport(cl, c("combind","ff.Wei.reg"))               
clusterEvalQ(cl, c(library(mvtnorm))) 

## parallel computing
ptm <- proc.time()
all_wei.reg1 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data, type=1)
all_wei.reg1 = as.data.frame(do.call(rbind, all_wei.reg1)) 
all_wei.reg2 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data, type=2)
all_wei.reg2 = as.data.frame(do.call(rbind, all_wei.reg2)) 
all_wei.reg3 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data, type=4)
all_wei.reg3 = as.data.frame(do.call(rbind, all_wei.reg3))
parApplyTime <- proc.time() - ptm # 1315.36 s 
stopCluster(cl)

 

ptm <- proc.time()
all_wei.reg1 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data, type=1)
all_wei.reg1 = as.data.frame(do.call(rbind, all_wei.reg1)) 
parApplyTime <- proc.time() - ptm # 1315.36 s 
stopCluster(cl)


#----------------------
aaaaa= ftable(xtabs(round( cbind(Mean.bias, Mean.rmse, Mean.ase, Rate, CP) ,4) ~ beta1+Time+rho, 
                    all_wei.reg3[all_wei.reg3$Scenario=="truncated" & all_wei.reg3$rho %in% c(-0.9, -0.5, 0, 0.5, 0.9),]))
aaaaa
library(memisc)
toLatex(aaaaa, digits = 4)


# 
# write.csv(x = all_wei.reg1, file = "reg_simu_1.csv",row.names = F)
# 
# write.csv(x = all_wei.reg2, file = "reg_simu_2.csv",row.names = F)
# 
# write.csv(x = all_wei.reg3, file = "reg_simu_3.csv",row.names = F)


#----------plot ------------

library(ggplot2)
library(ggnewscale)
ggplot(all_wei.reg3, aes(x= rho, group=Scenario, shape= Scenario)) + 
  geom_line(aes(y= 100*Mean.bias/true.mean, color=Scenario), size = 1.05) +
  scale_color_manual(name = '',values=c(censored='black', truncated="blue"))+
  new_scale_color() +
  geom_point(aes(y=100*Mean.bias/true.mean, color=Scenario, shape=Scenario), size=2.5) +
  scale_color_manual(values=c(censored='red', truncated="yellow1"),
                     guide = guide_legend(override.aes = list(shape = c(16,17))))+
  scale_shape_manual(values = c(censored = 16, truncated = 17)) + 
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red') +
  #  geom_point(aes(y=m.rate), shape=21, color="black",size=4) +
  scale_linetype_manual(name = "Reference : ", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red")))+
  # scale_size_manual(values=c(1.5, 1.5)) +
  guides( shape="none") +
  facet_grid(beta1 ~ Time,labeller= label_bquote( rows = beta[1]: .(beta1), cols= T[p]: .(Time)),
             scales='free') + 
  labs(color="Trendline: ",
       x= expression(paste(rho," : correlation")), 
       y= expression(paste("Relative Bias ", "( ", 100, "(",hat(mu),"/",mu - 1, " )%", ")")) ) +
  theme(legend.position="bottom", legend.text = element_text(size=10),
        strip.text = element_text(size = 15, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(face = "bold",hjust = 0.5,size=20))


ggplot(all_wei.reg3, aes(x= rho,group= Scenario, shape= Scenario)) + 
  geom_line(aes(y= 100*(Rate/true.rate-1), color=Scenario), size = 1.05) +
  scale_color_manual(name = '',values=c(censored='black', truncated="blue"))+
  new_scale_color() +
  geom_point(aes(y=100*(Rate/true.rate-1), color=Scenario, shape=Scenario),size=2.5) +
  scale_color_manual(values=c(censored='red', truncated="yellow1"),
                     guide = guide_legend(override.aes = list(shape = c(16,17))))+
  scale_shape_manual(values = c(censored = 16, truncated = 17)) + 
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red') +
  #  geom_point(aes(y=m.rate), shape=21, color="black",size=4) +
  scale_linetype_manual(name = "Reference : ", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red")))+
  # scale_size_manual(values=c(1.5, 1.5)) +
  guides( shape="none") +
  # guides(linetype= "none") +
  facet_grid(beta1 ~ Time,labeller= label_bquote( rows = beta[1]: .(beta1), cols= T[p]: .(Time)),
             scales='free') +
  labs(shape="Points: ", color="Trendline: ",
       x= expression(paste(rho," : correlation")), 
       y= expression(paste("Relative Bias ", "( ", 100, "(",hat(S), 
                           "(",T[p], ")", "/",S ,"(",T[p], ")", - 1, " )%", ")"))  ) +
  theme(legend.position="bottom", legend.text = element_text(size=10),
        strip.text = element_text(size = 15, face="bold"),
        axis.text=element_text(size=10),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(face = "bold",hjust = 0.5,size=20))
