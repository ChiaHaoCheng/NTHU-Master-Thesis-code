## boxcox

logL.box0 = function(para,y,xx, di){
  # if fix = 0, para = c(theta,beta0,beta1*I(type=1),beta2*I(type(1,2)))
  lambda= para[3]
  if(lambda==0){
    mu=para[1] + log(xx)* para[2]
  }
  else{
    mu=para[1] + (xx^lambda-1)/lambda * para[2]
  }
  r = sum(di)
  sig=1
  l=0
  z=(log(y)-mu)/sig
  for(i in 1:length(y)){
    l=l+di[i]*(-log(sig)- log(y[i])+z[i]-exp(z[i]) ) -
      (1-di[i])*exp((z[i]))
  }
  return(-l)
}

profile.result = optim(c(5, 1, 1),logL.box0,y=data$time,
      xx=data[,1],di=data$DEATH_EVENT, hessian=T, 
      control=list("maxit"=1e4))

aa =  profile.result$value

logL.lambda = function(para,lambda, y, xx, di){
  if(lambda==0){
    mu=para[1] + log(xx)* para[2]
  }
  else{
    mu=para[1] + (xx^lambda-1)/lambda * para[2]
  }
    r = sum(di)
    sig=1
    l=0
    z=(log(y)-mu)/sig
    for(i in 1:length(y)){
      l=l+di[i]*(-log(sig)- log(y[i])+z[i]-exp(z[i]) ) -
        (1-di[i])*exp((z[i]))
    }
    return(-l)
} 


lambda = seq(-3,3,length.out=100)
profile.value = 0
for (i in 1:100){
  optim = optim(c(mean(data$time), 1), fn = logL.lambda, y=data$time,
                xx = data$ejection_fraction, di =data$DEATH_EVENT,
                lambda=lambda[i], control=c("maxit"= 1e4))
  profile.value[i] = exp( aa - optim$value )
}

plot(lambda, profile.value, main = "ejection_fraction")
abline(h = exp(-qchisq(0.95,df = 1)/2), col="red")
## ejection_fraction : lambda = -2 order transformation



  profile.result2 = optim(c(mean(data$time), 1,1),logL.box0,y=data$time,
                        xx=data[,2],di=data$DEATH_EVENT, hessian=T,
                         control= list("maxit"=1e4))
  
  bb =  profile.result2$value
  lambda = seq(-2,2,length.out=50)
  profile.value2 = 0
  for (i in 1:50){
    optim = optim(c(mean(data$time), 1), fn = logL.lambda, y=data$time,
                  xx = data[,2], di =data$DEATH_EVENT, lambda=lambda[i],
                  control= list("maxit"=1e4))
    profile.value2[i] = exp( bb - optim$value )
  }
# 
 plot(x=lambda, y = profile.value2, main = "serum_creatinine")
 abline(h = exp(-qchisq(0.95,df = 1)/2), col="red")

## serum_creatinine : lambda = -1 order transformation