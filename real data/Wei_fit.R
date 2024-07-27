
# pre-processing : 
data = read.csv("heart_failure_clinical_records_dataset.csv",header=T)
data = data[,c(5,8,12,13)]  # Taking the index we interested
data$ejection_fraction = data$ejection_fraction/100 


# exponential fit:

logL.c = function(para,y,xx,type, di,t){
  # if fix = 0, para = c(theta,beta0,beta1*I(type=1),beta2*I(type(1,2)))
  n=length(y)
  if(type ==1){ mu = para[1] + colSums( t(xx) * c(para[2],para[3]) ) }
  else if(type ==2) { mu=para[1] + xx[,1] * para[2] }
  else if(type ==3) { mu=para[1] + xx[,2] * para[2] }
  else if(type ==4) { mu=rep(x = para,n) }
  
  r = sum(di)
  sig=1
  l=0
  z=(log(y)-mu)/sig
  zt = (log(t)-mu)/sig
  for(i in 1:length(y)){
    l=l+di[i]*(-log(sig)- log(y[i])+z[i]-exp(z[i]) ) -
      (1-di[i])*exp((zt[i]))
  }
  return(l)
}# X1 => ^(-2) ; X2 => ^(-1)

## for Model1
op.cen1 <- optim(c(10,0,0),logL.c,y=data$time,
                 xx=data[,1:2],type=1,di=data$DEATH_EVENT, hessian=T)
solve( op.cen1$hessian )
mu1 =  cbind(1, data[,1]^(-2), data[,2]^(-1)) %*% op.cen1$par
#mean(cbind(1,(1-data[,1]^(-2))/2,data[,2]) %*% op.cen1$par)
scale1 =  mean( exp(mu1) ) 
scale1 # sd= 11.6956
exp(-(241/scale1))


# #for model 2

op.cen2 <- optim(c(10,10),logL.c,y=data$time,
                 xx=data[,1:2],type=2,di=data$DEATH_EVENT, hessian=T)
solve( op.cen2$hessian )
mu2 = cbind(1, data[,1]^(-2)) %*% op.cen2$par
#mean(exp(cbind(1, (1-as.matrix(data[,1])^(-2))/2) %*% op.cen2$par))
scale2 =  mean( exp(mu2) )
scale2 # sd= 3.8674
#2*(op.cen2$value) + 2*2 # AIC for model 2

# #for model 3

op.cen3 <- optim(c(10,10),logL.c,y=data$time,
                 xx=data[,1:2],type=3,di=data$DEATH_EVENT,hessian=T)
solve( op.cen3$hessian )
mu3 = cbind(1, data[,2]^(-1)) %*% op.cen3$par
#mean(exp(cbind(1, as.matrix(data[,2])) %*% op.cen3$par))
scale3 = mean( exp(mu3) )
scale3 #sd=11.3490
exp(-(241/scale3))
2*(op.cen3$value) + 2*2# AIC for model 3


## for null model
op.cen4 <- optim(c(50),logL.c,y=data$time,
                 xx=data[,1:2],type=4,di=data$DEATH_EVENT,hessian=T)

mu4 = op.cen4$par
scale4 = exp(mu4)
scale4  #sd=41.4214
2*(op.cen4$value) + 2*1 # AIC for mode4



# Weibull fit

logL.W = function(para,y,xx,type, di){
  n=length(y)
  para.length = length(para)
  desigh.X = cbind(xx[,1]^(-2) , xx[,2]^(-1))
  if(type ==1){ mu = para[2] + colSums( t(desigh.X) * para[3:para.length] ) }
  else if(type ==2) { mu=para[2] + desigh.X[,1] * para[3] }
  else if(type ==3) { mu=para[2] + desigh.X[,2] * para[3] }
  else if(type ==4) { mu=rep( para[2] , n) }
  
  r = sum(di)
  sig=1/para[1]
  if(sig<=0) {sig =0.01}
  l=0
  z=(log(y)-mu)/sig
  for(i in 1:length(y)){
    l=l+di[i]*(-log(sig)- log(y[i])+z[i]-exp(z[i]) ) -
      (1-di[i])*exp((z[i]))
  }
  return(-l)
} # X1 => ^(-2) ; X2 => ^(-1)

# for Model1
op.cen1w <- optim(c(1,10,0,0),logL.W,y=data$time,
                 xx=data[,1:2],type=1,di=data$DEATH_EVENT, hessian=T)
var1w= solve( op.cen1w$hessian )
mu1w =  cbind(1, data[,1]^(-2), data[,2]^(-1)) %*% op.cen1w$par[-1]
#mean(cbind(1,(1-data[,1]^(-2))/2,data[,2]) %*% op.cen1$par)
scale1w =  exp(mu1w) 

mean1w = mean( scale1w* gamma(1+1/op.cen1w$par[1]) )
delta = cbind(digamma(1+1/op.cen1w$par[1])*(-1/op.cen1w$par[1]^2)*scale1w* gamma(1+1/op.cen1w$par[1]), 
              scale1w* gamma(1+1/op.cen1w$par[1]), 
              data[,1]^(-2)*scale1w* gamma(1+1/op.cen1w$par[1]), 
              data[,2]^(-1)*scale1w* gamma(1+1/op.cen1w$par[1]) )
sd.mean = apply(X = delta, MARGIN = 1, 
                FUN = function(x){x %*% var1w %*% x } )  
sd1w = sqrt(mean(sd.mean)/299)
mean1w
sd1w
#2*(op.cen1w$value) + 2*4 # AIC for model 1

# for model 2

op.cen2w <- optim(c(1,10,0),logL.W,y=data$time,
                 xx=data[,1:2],type=2,di=data$DEATH_EVENT, hessian=T)
var2w=solve( op.cen2w$hessian )
mu2w = cbind(1, data[,1]^(-2)) %*% op.cen2w$par[-1]
scale2w =  exp(mu2w) 

mean2w = mean( scale2w* gamma(1+1/op.cen2w$par[1]) )
delta = cbind(digamma(1+1/op.cen1w$par[1])*(-1/op.cen1w$par[1]^2)*scale2w* gamma(1+1/op.cen1w$par[1]), 
              scale2w* gamma(1+1/op.cen2w$par[1]), 
              data[,1]^(-2)*scale2w* gamma(1+1/op.cen2w$par[1]))
sd.mean = apply(X = delta, MARGIN = 1, 
                FUN = function(x){x %*% var2w %*% x } )  
sd2w = sqrt(mean(sd.mean)/299)
mean2w
sd2w

# for model 3

op.cen3w <- optim(c(1,10,0),logL.W,y=data$time,
                 xx=data[,1:2],type=3,di=data$DEATH_EVENT,hessian=T)
var3w= solve( op.cen3w$hessian )
mu3w = cbind(1, data[,2]^(-1)) %*% op.cen3w$par[-1]
scale3w =  exp(mu3w) 

mean3w = mean( scale3w* gamma(1+1/op.cen3w$par[1]) )
delta = cbind(digamma(1+1/op.cen3w$par[1])*(-1/op.cen3w$par[1]^2)*scale3w* gamma(1+1/op.cen3w$par[1]), 
              scale3w* gamma(1+1/op.cen3w$par[1]), 
              data[,2]^(-1)*scale3w* gamma(1+1/op.cen3w$par[1]))
sd.mean = apply(X = delta, MARGIN = 1, 
                FUN = function(x){x %*% var3w %*% x } )  
sd3w = sqrt(mean(sd.mean)/299)
mean3w
sd3w


# for null model
op.cen4w <- optim(c(1,10),logL.W,y=data$time,
                 xx=data[,1:2],type=4,di=data$DEATH_EVENT,hessian=T)
var4w=solve( op.cen4w$hessian )
mu4w = op.cen4w$par[2]
scale4w = exp(mu4w) 
mean4w = scale4w* gamma(1+1/op.cen4w$par[1])
delta = cbind(digamma(1+1/op.cen4w$par[1])*(-1/op.cen4w$par[1]^2)*scale4w* gamma(1+1/op.cen4w$par[1]), 
              scale4w* gamma(1+1/op.cen4w$par[1]))
sd.mean = delta %*% var4w %*% t(delta) 
sd4w = sqrt(sd.mean)
mean4w
sd4w
