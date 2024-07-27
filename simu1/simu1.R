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


############ simulation 1 ################

kappa = 10^(c(seq(from=log(1.01,base = 10),to=log(10,base = 10),length=10),
              seq(from=log(10.1,base = 10),to=log(100,base = 10),length=20)))
lambda = c(0.5,1,1.5) ; theta=1
n=1000

combind=expand.grid(lambda = lambda,theta=theta, n=n )

library(parallel)  # do paraelle computing
myCoreNums <- detectCores()
cl <- makeCluster(myCoreNums-1)
clusterExport(cl, c("combind","ff.Wei.uni","kappa","generate.data"))               
clusterEvalQ(cl, ) 

# apply(combind,1,generate.data, type=1)
ptm <- proc.time()
simu1 = parApply(cl, X = combind, MARGIN = 1 ,FUN = generate.data,kappa=kappa)
simu1 = as.data.frame(do.call(rbind, simu1)) 
parApplyTime <- proc.time() - ptm # 1315.36 s 
stopCluster(cl)
# simu1 =cbind( rep(c(0.5,1,1.5),each=100), simu1)
# colnames(simu1)[1] = "p.scale"
write.csv(simu1, file = "simu1_result.csv", row.names = F)

############plot ###############

simu1 = read.csv("simu1_result.csv",header = T)

# simu1 : 
## mean
library(ggplot2)
ggplot(simu1, aes(x= kappa ,group= Scenario )) + 
  geom_smooth(aes( y=Mean.bias/p.scale*100,color=Scenario),se = FALSE,method = loess) +
 # geom_smooth(aes(y=m.rate,color="moment"),linewidth=1.2,se = FALSE,method = loess) +
  scale_x_continuous(trans='log10',limits=c(1,100)) +
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red') +
  scale_linetype_manual(name = "Reference : ", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red"))) +
#  geom_text(aes(x=1.2,y=true.rate-0.01, label=round(true.rate,4)),size=5,colour="red") +
  scale_color_manual(values=c(censored='black', truncated="blue"))+
  facet_grid( rows = vars(p.scale), 
              labeller = label_bquote( rows = mu: .(p.scale)),
              scales = "free_y") +
  labs(color = "Linetype : ",
       x= expression(paste(log[10], "(",kappa, ")"," labelled with ", kappa)), 
       y= expression(paste("Relative Bias ", "( ", 100, "(",hat(mu),"/",mu - 1, " )%", ")")) ) +
  theme(legend.position="bottom",
        legend.text = element_text(size=10),
        strip.text = element_text(size = 15, face="bold"),
        axis.text=element_text(size=10, face="bold"),
        axis.title=element_text(size=14,face="bold"))

## rate
ggplot(simu1, aes(x= kappa ,group= Scenario )) + 
  geom_smooth(aes( y=100*(Rate/true.rate-1),color=Scenario),se = FALSE,method = loess) +
#  geom_smooth(aes(y=m.rate,color="moment"),linewidth=1.2,se = FALSE,method = loess) +
  scale_x_continuous(trans='log10',limits=c(1,100)) +
  geom_hline(aes(yintercept= 0, linetype = "Horizontal 0"), colour= 'red')+
  scale_linetype_manual(name = "Reference :", values = 2, 
                        guide = guide_legend(override.aes = list(color = "red"))) +
#  geom_text(aes(x=1.2,y=true.rate*0.95, label=round(true.rate,4)),size=3,colour="red") +
  scale_color_manual(values=c(censored='black', truncated="blue"))+
  facet_grid( rows = vars(p.scale), 
              labeller = label_bquote( rows = mu: .(p.scale)),
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




