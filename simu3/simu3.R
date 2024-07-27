# simu 3-1

ff.Wei.uni = function(n, mu, theta, k, t, B=500,fix=0){ 
  #nmu = multivariate mean ; sigma = multivariate covariance
  # fix=1: known shape ; 0 = unknown
  logL.c = function(para,x){
    ifelse(fix==1, mu <- log(para[1]), mu <- log(para[2]) )
    di= as.numeric(x <=t); r = sum(di)
    sig=ifelse(fix == 1,yes = 1/theta, no = 1/para[1])
    if(sig<0){sig=0.001}
    l=0
    z=(log(x)-mu)/sig
    zt = (log(t)-mu)/sig
    for(i in 1:length(x)){
      l=l+di[i]*( -log(sig)- log(x[i])+z[i]-exp(z[i]) ) -
        (1-di[i])*exp((zt))
    }
    return(l)
  }
  logL.t=function(para,x){
    ifelse(fix==1, mu <- log(para[1]), mu <- log(para[2]) )
    di= as.numeric(x <=t); r = sum(di)
    sig=ifelse(fix == 1,yes = 1/theta, no = 1/para[1])
    if(sig<0){sig=0.001}
    l=0
    z=(log(x)-mu)/sig
    zt = (log(t)-mu)/sig
    for(i in 1:length(x)){
      l=l+di[i]*( -log(sig)- log(x[i])+z[i]-exp(z[i]) - log(1-exp(-exp(zt))) )
    }
    return(l)
  } # truncated Weibull log L
  
  r.mu = rgamma(n,shape = k,scale = mu/k) #var=mu^2/k given Mu samplings
  simu = replicate( B, expr={
    x=0
    for(i in 1:n){ x[i] = rexp(1,rate = 1/r.mu[i]^theta)^(1/theta) }# data
    di= as.numeric(x <=t); r = sum(di)
    ifelse( fix==1, para <- c(mu), para <- c(theta,mu) )
    ifelse (  fix == 1 ,
              op.trun <- optim(para,logL.t,x=x,
                               control=list("fnscale"=-1),hessian=T,
                               method =  "Brent" , lower = 0.1, upper = 10) ,
              op.trun <- optim(para,logL.t,x=x,
                               control=list("fnscale"=-1),hessian=T) )
    
    ifelse( fix==1,mle.sha.t <- theta, mle.sha.t <- op.trun$par[1] )
    ifelse( fix==1,mle.sca.t <- op.trun$par, mle.sca.t <- op.trun$par[-1] )
    
    ratio = mle.sca.t/mu * gamma(1+1/mle.sha.t)/gamma(1+1/theta)# determine if discarding?
    while(ratio >=3  || ratio <= 1/3 || op.trun$convergence != 0){
      x=0
      for(i in 1:n){ x[i] = rexp(1,rate = 1/r.mu[i]^theta)^(1/theta) }# data
      di= as.numeric(x <=t); r = sum(di)
      ifelse( fix==1, para <- c(mu), para <- c(theta,mu) )
      ifelse (  fix == 1 ,
                op.trun <- optim(para,logL.t,x=x,
                                 control=list("fnscale"=-1),hessian=T,
                                 method =  "Brent" , lower = 0.1, upper = 5) ,
                op.trun <- optim(para,logL.t,x=x,
                                 control=list("fnscale"=-1),hessian=T) )
      
      ifelse( fix==1,mle.sha.t <- theta, mle.sha.t <- op.trun$par[1] )
      ifelse( fix==1,mle.sca.t <- op.trun$par, mle.sca.t <- op.trun$par[-1] )
      
      ratio = mle.sca.t/mu * gamma(1+1/mle.sha.t)/gamma(1+1/theta)# determine if discarding?
    }
    ifelse (  fix == 1 ,
              op.cen <- optim(para,logL.c,x=x,
                              control=list("fnscale"=-1),hessian=T,
                              method =  "Brent" , lower = 0.1, upper = 5) ,
              op.cen <- optim(para,logL.c,x=x,
                              control=list("fnscale"=-1),hessian=T)  )
    ifelse( fix==1,mle.sha.c <- theta, mle.sha.c <- op.cen$par[1] )
    ifelse( fix==1,mle.sca.c <- op.cen$par, mle.sca.c <- op.cen$par[-1] )
    
    ## cp( Wald for mu)
    
    if (fix == 1){
      Fisher.c = solve(-op.cen$hessian)
      Fisher.t = solve(-op.trun$hessian)
      ### censored part
      # glambda2=mle.sca.c*gamma(1+1/mle.sha.c)*digamma(1+1/mle.sha.c)*(-1/mle.sha.c^2)
      # glambda = matrix( c( glambda2, gamma(1+1/mle.sha.c) ),nrow = 1,ncol = 2)
      # V1 = glambda %*% Fisher.c %*% t(glambda)
      sd.c.mu = sqrt(Fisher.c) ; sd.c.theta= 0
      ub = mle.sca.c + qnorm(0.975)*sqrt(Fisher.c)
      lb = mle.sca.c - qnorm(0.975)*sqrt(Fisher.c)
      cp.c = as.numeric(mu <= ub && mu >= lb)
      
      ### truncated part
      # glambda21=mle.sca.t*gamma(1+1/mle.sha.t)*digamma(1+1/mle.sha.t)*(-1/mle.sha.t^2)
      # glambda1 = matrix(c( glambda21, gamma(1+1/mle.sha.t) ),nrow = 1,ncol = 2)
      # V2 = glambda1 %*% Fisher.t %*% t(glambda1)
      sd.t.mu = sqrt(Fisher.t) ; sd.t.theta= 0
      ub1 = mle.sca.t + qnorm(0.975)*sqrt(Fisher.t)
      lb1 = mle.sca.t - qnorm(0.975)*sqrt(Fisher.t)
      cp.t = as.numeric(mu<= ub1 && mu >= lb1)
    }
    
    else if( fix == 0 ){
      Fisher.c = solve(-op.cen$hessian)
      Fisher.t = solve(-op.trun$hessian)
      sd.c.mu = sqrt(Fisher.c)[2,2] ; sd.c.theta= sqrt(Fisher.c)[1,1]
      sd.t.mu = sqrt(Fisher.t)[2,2] ; sd.t.theta= sqrt(Fisher.t)[1,1]
      glambda2=mle.sca.c*gamma(1+1/mle.sha.c)*digamma(1+1/mle.sha.c)*(-1/mle.sha.c^2)
      glambda = matrix(c(gamma(1+1/mle.sha.c),glambda2),nrow = 1,ncol = 2)
      # censored part
      V1 = glambda %*% Fisher.c %*% t(glambda)
      ub = mle.sca.c*gamma(1+1/mle.sha.c) + qnorm(0.975)*sqrt(V1)
      lb = mle.sca.c*gamma(1+1/mle.sha.c) - qnorm(0.975)*sqrt(V1)
      cp.c = as.numeric(mu*gamma(1+1/theta)<= ub && mu*gamma(1+1/theta) >= lb)
      # truncated part
      glambda21=mle.sca.t*gamma(1+1/mle.sha.t)*digamma(1+1/mle.sha.t)*(-1/mle.sha.t^2)
      glambda1 = matrix(c(gamma(1+1/mle.sha.t),glambda21),nrow = 1,ncol = 2)
      V2 = glambda1 %*% Fisher.t %*% t(glambda1)
      ub1 = mle.sca.t*gamma(1+1/mle.sha.t) + qnorm(0.975)*sqrt(V2)
      lb1 = mle.sca.t*gamma(1+1/mle.sha.t) - qnorm(0.975)*sqrt(V2)
      cp.t = as.numeric(mu*gamma(1+1/theta)<= ub1 && mu*gamma(1+1/theta) >= lb1)
    }
    
    
    return(c("cen.mu" =mle.sca.c, "cen.theta"=mle.sha.c,
             "trun.mu" =mle.sca.t,"trun.theta"=mle.sha.t,
             "count"=r,
             "cp.censored" = cp.c, "cp.truncated" = cp.t, 
             "sd.c.mu" = sd.c.mu, "sd.c.theta"=sd.c.theta,
             "sd.t.mu" = sd.t.mu, "sd.t.theta"=sd.t.theta) )
  })
  
  
  ## censored results
  est.c.lambda = mean( simu[1,] ) ; est.c.theta = mean( simu[2,] ) 
  rmse.c.mu = sqrt( mean((simu[1,] - est.c.lambda)^2) )
  rmse.c.theta = sqrt( mean((simu[2,] - est.c.theta)^2) )
  sd.c.mu = mean(simu[8,]) ; sd.c.theta = mean(simu[9,])
  est.mean.c = mean(simu[1,] * gamma(1+ 1/ simu[2,]) )
  mle.rate.censored = mean( exp(-( t/simu[1, ] )^simu[2,] ) )
  
  ## truncated results
  est.t.lambda = median( simu[3,] ) ; est.t.theta = median( simu[4,] ) 
  rmse.t.mu = sqrt( mean((simu[3,] - est.t.lambda)^2) )
  rmse.t.theta = sqrt( mean((simu[4,] - est.t.theta)^2) )
  sd.t.mu = mean(simu[10,]) ; sd.t.theta = mean(simu[11,])
  est.mean.t = median(simu[3,] * gamma(1+ 1/ simu[4,]) )
  mle.rate.truncated = mean( exp(-( t/simu[3, ] )^simu[4,] ) )
  
  mean.bias.c = est.mean.c - mu * gamma(1+1/theta) 
  mean.bias.t = est.mean.t - mu * gamma(1+1/theta)
  mean.rate = 1-mean(simu[5,]/n)
  
  WaldCP.c = mean( simu[6,] )
  WaldCP.t = mean( simu[7,] )
  
  return(list("Censored.para"= c(est.c.lambda, est.c.theta, rmse.c.mu, rmse.c.theta),
              "Truncated.para"= c(est.t.lambda, est.t.theta, rmse.t.mu, rmse.t.theta),
              "Censored.sd" = c(sd.c.mu, sd.c.theta),
              "Truncated.sd"=c(sd.t.mu, sd.t.theta),
              "Censored.mean.bias"= mean.bias.c,
              "Truncated.mean.bias"= mean.bias.t,
              "est.rate.censored" = mle.rate.censored,
              "est.rate.truncated" = mle.rate.truncated,
              "mean.rate" = mean.rate,
              "CP" = c(WaldCP.c, WaldCP.t)))
  
}

generate.data = function(comb,kappa,t=1, fix=1){
  c.sca = c.sha = t.sca = t.sha = 0
  r.c.sca = r.c.sha = r.t.sca = r.t.sha = 0
  sd.c.sca = sd.c.sha=sd.t.sca=sd.t.sha=0
  c.bias =t.bias =0;c.rate =t.rate =all.rate= 0
  CP.c=CP.t=0
  for (i in 1:length(kappa)){
    bias = ff.Wei.uni(n=comb[3],mu = comb[1],theta=comb[2],k = kappa[i],t=t, fix=fix) 
    c.sca[i] = bias$Censored.para[1] ; c.sha[i] = bias$Censored.para[2]
    r.c.sca[i] = bias$Censored.para[3] ; r.c.sha[i] = bias$Censored.para[4]
    sd.c.sca[i] = bias$Censored.sd[1] ; sd.c.sha[i] = bias$Censored.sd[2]
    t.sca[i] = bias$Truncated.para[1] ; t.sha[i] = bias$Truncated.para[2]
    r.t.sca[i] = bias$Truncated.para[3] ; r.t.sha[i] = bias$Truncated.para[4]
    sd.t.sca[i] = bias$Truncated.sd[1] ; sd.t.sha[i] = bias$Truncated.sd[2]
    c.bias[i] = bias$Censored.mean.bias
    #  c.asy[i] = bias$Censored.asy.mean.bias
    t.bias[i] = bias$Truncated.mean.bias
    #  t.asy[i] = bias$Truncated.asy.mean.bias
    c.rate[i] = bias$est.rate.censored 
    t.rate[i] = bias$est.rate.truncated 
    all.rate[i] = bias$mean.rate
    CP.c[i] = bias$CP[1] ; CP.t[i] = bias$CP[2]
  }
  true.rate = exp(-(t/comb[1])^(comb[2]) )
  final1 = data.frame( "p.scale"=comb[1], "p.theta"=comb[2], "kappa" = kappa, 
                       "scale" = c.sca , "r.scale" =r.c.sca, "sd.scale"=sd.c.sca,
                       "shape" = c.sha , "r.shape" =r.c.sha, "sd.shape"=sd.c.sha,
                       "Mean.bias" = c.bias, "Rate" = c.rate, 
                       "m.rate" = all.rate, "true.rate" = true.rate,
                       "CP" = CP.c, "Scenario" = "censored")
  final2 = data.frame( "p.scale"=comb[1], "p.theta"=comb[2], "kappa" = kappa, 
                       "scale" = t.sca , "r.scale" =r.t.sca, "sd.scale"=sd.t.sca,
                       "shape" = t.sha , "r.shape" =r.t.sha, "sd.shape"=sd.t.sha,
                       "Mean.bias" = t.bias, "Rate" = t.rate, 
                       "m.rate" = all.rate, "true.rate" = true.rate,
                       "CP" = CP.t, "Scenario" = "truncated")
  final = rbind( final1, final2)
  return(final)
}

## settings:
kappa = 10^(c(seq(from=log(1.01,base = 10),to=log(10,base = 10),length=10),
              seq(from=log(10.1,base = 10),to=log(100,base = 10),length=20)))
lambda = c(0.5,1,1.5) ; theta=c(0.5,1,1.5)
n=1000

combind=expand.grid(lambda = lambda,theta=theta, n=n )

# paraelle computing
library(parallel)  
myCoreNums <- detectCores()
cl <- makeCluster(myCoreNums-1)
clusterExport(cl, c("combind","ff.Wei.uni","kappa","generate.data"))               

ptm <- proc.time()
simu3 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data,kappa=kappa, fix=0)
simu3 = as.data.frame(do.call(rbind, simu3)) 
parApplyTime <- proc.time() - ptm # 1978 s 
stopCluster(cl)

## Save 
setwd("C:\\Users\\stat_835\\Desktop\\模擬\\simu3")
write.csv(simu3, file = "simu3_uni_result.csv", row.names = F)
simu3= read.csv("simu3_uni_result.csv",header=T)

## convert to latex table
aaaaa = ftable(xtabs(cbind( scale, r.scale, sd.scale, Mean.bias, CP ,Rate) ~ p.scale+ kappa, 
                     data=simu3, subset = simu3$kappa %in% c(100, 10, 1.01) & simu3$Scenario == "truncated"))
library(memisc)
toLatex(aaaaa, digits = 4)

################## plot ###############################
ggplot(simu3, aes(x= kappa ,group= Scenario ) ) + 
  geom_smooth(aes( y=Mean.bias/(p.scale*gamma(1+1/p.theta))*100,color=Scenario),se = FALSE,method = loess) +
  # geom_smooth(aes(y=m.rate,color="moment"),linewidth=1.2,se = FALSE,method = loess) +
  scale_x_continuous(trans='log10',limits=c(1,100)) +
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red') +
  scale_linetype_manual(name = "Reference : ", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red"))) +
  #  geom_text(aes(x=1.2,y=true.rate-0.01, label=round(true.rate,4)),size=5,colour="red") +
  scale_color_manual(values=c(censored='black', truncated="blue"))+
  facet_grid( rows = vars(p.theta), cols = vars(p.scale),
              labeller = label_bquote( cols = mu: .(p.scale), rows = theta: .(p.theta)),
              scales = "free_y") +
  labs(color = "Linetype : ",
       x= expression(paste(log[10], "(",kappa, ")"," labelled with ", kappa)), 
       y= expression(paste("Relative Bias ", "( ", 100, "(",hat(mu),"/",mu - 1, " )%", ")")) ) +
  theme(legend.position="bottom",
        legend.text = element_text(size=10),
        strip.text = element_text(size = 15, face="bold"),
        axis.text=element_text(size=10, face="bold"),
        axis.title=element_text(size=14,face="bold"))


ggplot(simu3, aes(x= kappa ,group= Scenario )) + 
  geom_smooth(aes( y=100*(Rate/true.rate-1),color=Scenario),se = FALSE,method = loess) +
  #  geom_smooth(aes(y=m.rate,color="moment"),linewidth=1.2,se = FALSE,method = loess) +
  scale_x_continuous(trans='log10',limits=c(1,100)) +
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red')+
  scale_linetype_manual(name = "Reference :", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red"))) +
  #  geom_text(aes(x=1.2,y=true.rate*0.95, label=round(true.rate,4)),size=3,colour="red") +
  scale_color_manual(values=c(censored='black', truncated="blue"))+
  facet_grid( rows = vars(p.theta), cols = vars(p.scale),
              labeller = label_bquote( cols = mu: .(p.scale), rows = theta: .(p.theta)),
              scales = "free_y") +
  labs(color = "Linetype : ",
       x= expression(paste(log[10], "(",kappa, ")"," labelled with ", kappa)), 
       y= expression(paste("Relative Bias ", "( ", 100, "(",hat(S), 
                           "(",T[p], ")", "/",S ,"(",T[p], ")", - 1, " )%", ")"))  ) +
  theme(legend.position="bottom",
        legend.text = element_text(size=10),
        strip.text = element_text(size = 15, face="bold"),
        axis.text=element_text(size=10, face="bold"),
        axis.title=element_text(size=14,face="bold"))


# for simu 3-2 : mixed analysis; 

ff.Wei.mix = function(n,theta,nmu,sigma,b0,b1,k,t,B=500,type=1,fix=0,seed=1){ 
  #nmu = multivariate mean ; sigma = multivariate covariance
  # fix=1: known shape ; 0 = unknown
  logL.c = function(para,x,xx,type,fix){
    # if fix = 0, para = c(theta,beta0,beta1*I(type=1),beta2*I(type(1,2)))
    
    if (fix == 0){
      if(type ==1){ mu=para[2] + xx[,1] * para[3] }
      else if(type ==2) { mu=rep( para[2] , n) }
    }
    else if (fix == 1){
      if(type ==1){ mu=para[1] + xx[,1] * para[2] }
      else if(type ==2) { mu=rep(x = para,n) }
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
      if(type ==1){ mu=para[2] + xx[,1] * para[3] }
      else if(type ==2) { mu=rep( para[2] , n) }
    }
    else if (fix == 1){
      if(type ==1){ mu=para[1] + xx[,1] * para[2] }
      else if(type ==2) { mu=rep(x = para, n) }
    }
    di= as.numeric(x <=t); r = sum(di)
    sig=ifelse(fix == 1,yes = 1/theta, no = 1/para[1])
    if(sig<0){sig=0.001}
    l=0
    z=(log(x)-mu)/sig
    zt = (log(t)-mu)/sig
    for(i in 1:length(x)){
      l=l+di[i]*(-log(sig)- log(x[i])+z[i]-exp(z[i]) - log(1-exp(-exp(zt[i]))))
    }
    return(l)
  } # truncated Weibull log L
  
  #  sigma.matrix = matrix(c(0.1,p*0.1,p*0.1,0.1),2,2)
  set.seed(seed) ; bigX = matrix( rnorm( n, mean = nmu, sd = sigma) ,ncol = 1)
  r.mu = rgamma(n,shape = k,scale = 1/k) * exp( b0 + bigX %*% b1 )
  true.mu = exp( b0 + nmu*b1 + 1/2* ( b1*sigma)^2 )
  ifelse( type== 2, cutoff <- 10, cutoff <- 3 )
  simu = replicate( B, expr={
    x=0
    for(i in 1:n){ x[i] = rexp(1,rate = 1/r.mu[i]^theta)^(1/theta) }# data
    di= as.numeric(x <=t); r = sum(di)
    if( type == 1 ) { ifelse(fix==1,para <- c(b0,b1), para <- c(theta,b0,b1) )}
    else if ( type == 2 ){ ifelse(fix==1,para <- b0, para <- c(theta, b0) ) }
    
    ifelse ( type == 2 && fix == 1 ,
             op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                              control=list("fnscale"=-1),hessian=T,
                              method =  "Brent" , lower = 0.5, upper = 5) ,
             op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                              control=list("fnscale"=-1),hessian=T) )
    
    ifelse( fix==1, mle.sha.t <- theta, mle.sha.t <- op.trun$par[1] )
    ifelse( fix==1, beta.hat.t <- op.trun$par, beta.hat.t <- op.trun$par[-1] )
    
    if (type == 1 ) {  adj.X = cbind(1, bigX) }
    else if (type == 2) {  adj.X = 1 }
    mle.sca.t = sum( exp( adj.X %*% beta.hat.t)* di )/ r 
    ratio = mle.sca.t/true.mu * gamma(1+1/mle.sha.t)/gamma(1+1/theta)# determine if discarding?
    while(ratio >=cutoff  || ratio <= 1/cutoff || op.trun$convergence != 0){
      x=0
      for(i in 1:n){ x[i] = rexp(1,rate = 1/r.mu[i]^theta)^(1/theta) }# data
      di= as.numeric(x <=t); r = sum(di)
      if( type == 1 ) { ifelse(fix==1,para <- c(1,1), para <- c(theta,1,1) )}
      else if ( type == 2 ){ ifelse(fix==1,para <- b0, para <- c(theta,1) ) }
      
      ifelse ( type == 2 && fix == 1 ,
               op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                                control=list("fnscale"=-1),hessian=T,
                                method =  "Brent" , lower = 0.5, upper = 5) ,
               op.trun <- optim(para,logL.t,x=x,xx=bigX,type=type,fix=fix,
                                control=list("fnscale"=-1),hessian=T) )
      
      ifelse( fix==1,mle.sha.t <- theta, mle.sha.t <- op.trun$par[1] )
      ifelse( fix==1,beta.hat.t <- op.trun$par,beta.hat.t <- op.trun$par[-1] )
      
      if (type == 1 ) {  adj.X = cbind(1, bigX) }
      else if (type == 2) {  adj.X = 1 }
      mle.sca.t = sum( exp( adj.X %*% beta.hat.t)* di )/ r 
      ratio = mle.sca.t/true.mu * gamma(1+1/mle.sha.t)/gamma(1+1/theta)# determine if discarding?
    }
    ifelse ( type == 2 && fix == 1 ,
             op.cen <- optim(para,logL.c,x=x,xx=bigX,type=type,fix=fix,
                             control=list("fnscale"=-1),hessian=T,
                             method =  "Brent" , lower = 0.5, upper = 10) ,
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
        c.var.mu = apply(X = cbind(1, bigX), MARGIN = 1, FUN = function(x){
          (exp(x %*% beta.hat.c) * gamma(1+1/mle.sha.c))^2 *
            (c(digamma(1+1/mle.sha.c)*(-1/mle.sha.c^2) ,x) %*% 
               Fisher.c %*% c( digamma(1+1/mle.sha.c)*(-1/mle.sha.c^2), x) )
        })
        t.var.mu = apply(X = cbind(1, bigX, di), MARGIN = 1, FUN = function(x){
          (exp(x[-3] %*% beta.hat.t) * gamma(1+1/mle.sha.t))^2 *
            (c(digamma(1+1/mle.sha.t)*(-1/mle.sha.t^2), x[-3]) %*% 
               Fisher.t %*% c(digamma(1+1/mle.sha.t)*(-1/mle.sha.t^2), x[-3]))
        })
        sd.cen = sqrt( mean(c.var.mu) ) ; sd.trun = sqrt( sum(t.var.mu)/n )
        c.CI.u = mle.sca.c*gamma(1+1/mle.sha.c) + qnorm(0.975)*sqrt(sum(c.var.mu)/n)
        c.CI.l = mle.sca.c*gamma(1+1/mle.sha.c) - qnorm(0.975)*sqrt(sum(c.var.mu)/n)
        t.CI.u = mle.sca.t*gamma(1+1/mle.sha.t) + qnorm(0.975)*sqrt(sum(t.var.mu)/n)
        t.CI.l = mle.sca.t*gamma(1+1/mle.sha.t) - qnorm(0.975)*sqrt(sum(t.var.mu)/n)  
      }
      else if ( type == 2){
        # censored part
        c.var.mu = exp(2*beta.hat.c) * gamma(1+1/mle.sha.c)^2 *
           (c(digamma(1+1/mle.sha.c)*(-1/mle.sha.c^2) ,1) %*% 
              Fisher.c %*% c( digamma(1+1/mle.sha.c)*(-1/mle.sha.c^2), 1) )
        
        c.CI.u = mle.sca.c*gamma(1+1/mle.sha.c) + qnorm(0.975)*sqrt(c.var.mu)
        c.CI.l = mle.sca.c*gamma(1+1/mle.sha.c) - qnorm(0.975)*sqrt(c.var.mu)
        # truncated part
        t.var.mu = exp( beta.hat.t*2) * gamma(1+1/mle.sha.t)^2 *
          (c(digamma(1+1/mle.sha.t)*(-1/mle.sha.t^2) ,1) %*% 
             Fisher.t %*% c( digamma(1+1/mle.sha.t)*(-1/mle.sha.t^2), 1) )
        t.CI.u = mle.sca.t*gamma(1+1/mle.sha.t) + qnorm(0.975)*sqrt(t.var.mu)
        t.CI.l = mle.sca.t*gamma(1+1/mle.sha.t) - qnorm(0.975)*sqrt(t.var.mu)
        
        sd.cen = sqrt( c.var.mu ) ; sd.trun = sqrt( t.var.mu )
      }
    }
    
    else if( fix == 1 ){
      Fisher.c = solve(-op.cen$hessian)
      Fisher.t = solve(-op.trun$hessian)
      if ( type == 1 ){
        c.var.mu = apply(X = adj.X, MARGIN = 1, FUN = function(x){
          exp(x %*% beta.hat.c)^2 * (x %*% Fisher.c %*% x)
        })
        t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
          exp(x[-3] %*% beta.hat.t)^2 * (x[-3] %*% Fisher.t %*% x[-3]) * x[3]
        })
        c.CI.l =  mean( exp(adj.X %*% beta.hat.c) - qnorm(0.975) * sqrt(c.var.mu) )
        c.CI.u =  mean( exp(adj.X %*% beta.hat.c) + qnorm(0.975) * sqrt(c.var.mu) )
        t.CI.l =  mean( exp(adj.X %*% beta.hat.t) - qnorm(0.975) * sqrt(t.var.mu) )
        t.CI.u =  mean( exp(adj.X %*% beta.hat.t) + qnorm(0.975) * sqrt(t.var.mu) )
      }
      else if ( type == 2){
        c.var.mu = exp(beta.hat.c)^2 * Fisher.c
        t.var.mu = apply(X = cbind(adj.X, di), MARGIN = 1, FUN = function(x){
          exp(x[-2] * beta.hat.t)^2 * (x[-2]^2 * Fisher.t) * x[2]
        })
        c.CI.l =  mean( exp(adj.X %*% beta.hat.c) - qnorm(0.975) * sqrt(c.var.mu) )
        c.CI.u =  mean( exp(adj.X %*% beta.hat.c) + qnorm(0.975) * sqrt(c.var.mu) )
        t.CI.l =  mean( exp(adj.X %*% beta.hat.t) - qnorm(0.975) * sqrt(t.var.mu) )
        t.CI.u =  mean( exp(adj.X %*% beta.hat.t) + qnorm(0.975) * sqrt(t.var.mu) )
      }
    }
    cp.c = as.numeric( true.mu*gamma(1+1/theta) >= c.CI.l && true.mu*gamma(1+1/theta) <= c.CI.u )
    cp.t = as.numeric( true.mu*gamma(1+1/theta) >= t.CI.l && true.mu*gamma(1+1/theta) <= t.CI.u )
    
    if ( type == 1 ){
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
    else if (type == 2){
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
  
  read.mu = exp( cbind(1, bigX) %*% c(b0,b1) )
  if ( type == 1 ){
    est.c.b0 = mean(simu[1,])
    est.c.b1 = mean(simu[2,])
    est.c.beta = c(est.c.b0,est.c.b1)
    # if (type == 2) { adjX = bigX[,1]}
    # else { adjX = bigX[,2] }
    est.c.mu = mean( exp( cbind(1, bigX) %*% est.c.beta) )
    #est.c.mu = rowMeans( exp( cbind(1,bigX) %*% simu[1:2,] ) ) 
 #   est.c.mu =  mean(simu[3,])
    est.c.theta = mean(simu[4,])
    rmse.c.mean = sd(simu[3,]*gamma(1+1/simu[4,]))
    a.c.sd = mean(simu[12,] )
    mle.rate.censored = exp(-(t/est.c.mu)^est.c.theta) 
    
    est.t.b0 = mean(simu[5,])
    est.t.b1 = mean(simu[6,])
    est.t.beta = c(est.t.b0,est.t.b1)
    
    #est.t.mu = apply( X = cbind(1, bigX) %*% simu[5:6,], MARGIN = 1, FUN = function(x) median(exp(x)) )  
    #est.t.mu = median( simu[7,])
    est.t.mu = mean( exp( cbind(1, bigX) %*% est.t.beta ) )
    est.t.theta = mean(simu[8,])
    rmse.t.mean = sd(simu[7,]*gamma(1+1/simu[8,]))
    a.t.sd = mean(simu[13,] )
    mle.rate.truncated =  exp(-(t/est.t.mu)^est.t.theta)
    
    mean.bias.c = est.c.mu * gamma(1+1/est.c.theta) - true.mu * gamma(1+1/theta)  
    mean.bias.t = est.t.mu * gamma(1+1/est.t.theta) - true.mu * gamma(1+1/theta) 
    mean.rate = 1-mean(simu[9,]/n)
    
    WaldCP.c = mean( simu[10,] )
    WaldCP.t = mean( simu[11,] )
  }
 
  else if (type ==2 ) {
    est.c.b0 = mean(simu[1,])
    est.c.beta = est.c.b0
    
    est.c.mu =  exp(est.c.b0)
    est.c.theta = mean(simu[3,])
    rmse.c.mean = sd(simu[2,]*gamma(1+1/simu[3,]))
    a.c.sd = mean(simu[10,] )
    mle.rate.censored = exp(-(t/est.c.mu )^est.c.theta) 
    
    est.t.b0 = mean(simu[4,])
    est.t.beta = est.t.b0
    
    est.t.mu = exp(est.t.b0)
    est.t.theta = mean(simu[6,])
    rmse.t.mean = sd(simu[5,]*gamma(1+1/simu[6,]))
    a.t.sd = mean(simu[11,] )
    mle.rate.truncated = exp(-(t/est.t.mu )^est.t.theta) 
    
    mean.bias.c =  est.c.mu * gamma(1+1/est.c.theta) - true.mu * gamma(1+1/theta) 
    mean.bias.t =  est.t.mu * gamma(1+1/est.t.theta) - true.mu * gamma(1+1/theta) 
    mean.rate = 1-mean(simu[7,]/n)
    
    WaldCP.c = mean( simu[8,] )
    WaldCP.t = mean( simu[9,] )
  }
  return(list("bigX" = bigX,
              "Censored.para"= c(est.c.beta, est.c.theta),
              "Truncated.para"= c(est.t.beta, est.t.theta),
              "Censored.sesd" = c(rmse.c.mean, a.c.sd),
              "Truncated.sesd" = c(rmse.t.mean, a.t.sd),
              "Censored.mean.bias"= mean.bias.c,
              "Truncated.mean.bias"= mean.bias.t,
              "est.rate.censored" = mle.rate.censored,
              "est.rate.truncated" = mle.rate.truncated,
              "mean.rate" = mean.rate,
              "CP" = c(WaldCP.c, WaldCP.t)))
  
}
# setting:
t = c(3, 5, 10)
beta1 = c(0.5, 0.9)
theta = c(0.5,1,1.5)
library(tidyr)  # rearrang the order of combination
combind=crossing(Time = t,beta1=beta1, theta= theta)
n=1000
kappa = 10^(c(seq(from=log(1.01,base = 10),to=log(10,base = 10),length=10),
              seq(from=log(10.1,base = 10),to=log(100,base = 10),length=20)))

generate.data = function(combind,type,kappa){
  comb = as.matrix(combind)
#  sigma = 1
  true.mu = exp( 1 + 1/2* comb[2]^2 * 1^2 )
#  c.sca = c.sha = t.sca = t.sha = 0
#  r.c.sca = r.c.sha = r.t.sca = r.t.sha = 0
  c.bias = t.bias =0;c.rate =t.rate =all.rate= true.rate=0
  c.rmse = t.rmse = c.ase = t.ase = 0
  CP.c=CP.t=0
  for (i in 1:length(kappa)){
    result = ff.Wei.mix(n = 1000,theta = comb[3],b0=1,b1=comb[2],t=comb[1],
                        nmu = 0, sigma = 1, type = type, fix=0,k = kappa[i])
    # c.sca[i] = bias$Censored.para[1] ; c.sha[i] = bias$Censored.para[2]
    # r.c.sca[i] = bias$Censored.para[3] ; r.c.sha[i] = bias$Censored.para[4]
    # t.sca[i] = bias$Truncated.para[1] ; t.sha[i] = bias$Truncated.para[2]
    # r.t.sca[i] = bias$Truncated.para[3] ; r.t.sha[i] = bias$Truncated.para[4]
    c.bias[i] = result$Censored.mean.bias
    c.rmse[i] = result$Censored.sesd[1]; c.ase[i] = result$Censored.sesd[2]
    #  c.asy[i] = bias$Censored.asy.mean.bias
    t.bias[i] = result$Truncated.mean.bias
    t.rmse[i] = result$Truncated.sesd[1]; t.ase[i] = result$Truncated.sesd[2]
    #  t.asy[i] = bias$Truncated.asy.mean.bias
    c.rate[i] = result$est.rate.censored
    t.rate[i] = result$est.rate.truncated
    all.rate[i] = result$mean.rate
    CP.c[i] = result$CP[1] ; CP.t[i] = result$CP[2]
  }
  true.rate = exp(-(comb[1]/true.mu)^(comb[3]) )
  true.mean = true.mu* gamma(1+1/comb[3])
  final1 = data.frame( "p.theta"=comb[3], "Time"= comb[1], "beta" = comb[2], "kappa" = kappa, 
                       "Mean.bias" = c.bias, "Rate" = c.rate, 
                       "Mean.rmse"=c.rmse,"Mean.ase"=c.ase,
                       "m.rate" = all.rate, "true.rate" = true.rate,"true.mean"=true.mean,
                       "CP" = CP.c, "Scenario" = "censored")
  final2 = data.frame( "p.theta"=comb[3], "Time"= comb[1], "beta" = comb[2], "kappa" = kappa,
                       "Mean.bias" = t.bias, "Rate" = t.rate, 
                       "Mean.rmse"=t.rmse,"Mean.ase"=t.ase,
                       "m.rate" = all.rate, "true.rate" = true.rate,"true.mean"=true.mean,
                       "CP" = CP.t, "Scenario" = "truncated")
  final = rbind( final1, final2)
  
  return(final)
}
#eqwe = generate.data(combind = combind[2,], type = 1,kappa = kappa[1:10])
# test = apply(combind[28,],1,generate.data, type=4)
# test = generate.data(combind[30,], type=4)

library(parallel)  # do paraelle computing
myCoreNums <- detectCores()
cl <- makeCluster(myCoreNums-1)
clusterExport(cl, c("combind","ff.Wei.mix","kappa","generate.data","n"))               
#clusterEvalQ(cl, c(library(mvtnorm))) 

# apply(combind,1,generate.data, type=1)
ptm <- proc.time()
simu3_1 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data,kappa=kappa,type=1)
simu3_1 = as.data.frame(do.call(rbind, simu3_1)) 
simu3_2 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data,kappa=kappa,type=2)
simu3_2 = as.data.frame(do.call(rbind, simu3_2)) 
parApplyTime <- proc.time() - ptm # 1315.36 s 
stopCluster(cl)

write.csv(simu3_1, file = "simu32_1_result.csv", row.names = F)
write.csv(simu3_2, file = "simu32_2_result.csv", row.names = F)
#############latex###############

aaaaa = ftable(xtabs(cbind( Mean.bias, Mean.rmse, Mean.ase, Rate, CP) ~ p.theta + Time +beta + kappa, 
                     data=simu3_2, subset = simu3_2$kappa %in% c(100, 10, 1.01) & simu3_2$Scenario == "truncated"))
aaaaa
library(memisc)
toLatex(aaaaa, digits = 4)

#############plot################

simu3_1 = read.csv("simu32_1_result.csv",header = T)
simu3_2 = read.csv("simu32_2_result.csv",header = T)
# true.mu = apply(X = simu3_1[, c(1,3)], MARGIN = 1, FUN = function(x){
#   mu = 1 ; AA = x[2]^2
#   exp(mu + 1/2*AA)* gamma(1+1/x[1])
# })
# simu3_1 = cbind(simu3_1, "true.mean"=true.mu)
# simu3_2 = cbind(simu3_2, "true.mean"=true.mu)

library(ggplot2)
ggplot(simu3_2, aes(x= kappa ,group= Scenario ) ) + 
  geom_smooth(aes( y=100*Mean.bias/true.mean,color=Scenario),se = FALSE,method = loess) +
  # geom_smooth(aes(y=m.rate,color="moment"),linewidth=1.2,se = FALSE,method = loess) +
  scale_x_continuous(trans='log10',limits=c(1,100)) +
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red') +
  scale_linetype_manual(name = "Reference : ", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red"))) +
  #  geom_text(aes(x=1.2,y=true.rate-0.01, label=round(true.rate,4)),size=5,colour="red") +
  scale_color_manual(values=c(censored='black', truncated="blue"))+
  facet_grid( Time + beta ~ p.theta,
              labeller = label_bquote(
                cols = theta:.(p.theta), 
                rows = atop(T[p]:.(Time), beta[1]:.(beta))) ,
              scales = "free_y") +
  labs(shape="Points: ", color="Trendline: ",
       x= expression(paste(log[10], "(",kappa, ")"," labelled with ", kappa)), 
       y= expression(paste("Relative Bias ", "( ", 100, "(",hat(mu),"/",mu - 1, " )%", ")")) ) +
  theme(legend.position="bottom",
        legend.text = element_text(size=10),
        strip.text = element_text(size = 15, face="bold"),
        strip.text.y = element_text(angle = 0, size = 12),
        axis.text=element_text(size=10, face="bold"),
        axis.title=element_text(size=14,face="bold"))

## rate
#library(faraway)
ggplot(simu3_2, aes(x= kappa ,group= Scenario )) + 
  geom_smooth(aes( y=Rate/true.rate-1,color=Scenario),se = FALSE,method = loess) +
#  geom_smooth(aes(y=m.rate,color="moment"),linewidth=1.2,se = FALSE,method = loess) +
  scale_x_continuous(trans='log10',limits=c(1,100)) +
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red') +
  scale_linetype_manual(name = "Reference : ", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red"))) +
#  geom_text(aes(x=1.2,y=true.rate*0.95, label=round(true.rate,4)),size=2.5,colour="red") +
  scale_color_manual(values=c(censored='black', truncated="blue"))+
  facet_grid( Time + beta ~ p.theta,
              labeller = label_bquote(
                cols = theta:.(p.theta), 
                rows = atop(T[p]:.(Time), beta[1]:.(beta))) ,
              scales = "free_y") +
  labs(shape="Points: ", color="Linetype: ",
       x= expression(paste(log[10], "(",kappa, ")"," labelled with ", kappa)), 
       y= expression(paste("Relative Bias ", "( ", hat(S), 
                           "(",T[p], ")", "/",S ,"(",T[p], ")", - 1, ")"))  ) +
  theme(legend.position="bottom",
        legend.text = element_text(size=10),
        strip.text = element_text(size = 15, face="bold"),
        strip.text.y = element_text(angle = 0, size = 12),
        axis.text=element_text(size=10, face="bold"),
        axis.title=element_text(size=14,face="bold"))
