# see https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html to learn more
dldat <- read.csv("~/GitHub/NonNegativeLASSO-PortfolioReplication/Data/dldat.csv", row.names=1)
require(glmnet)
a=colSums(is.na(dldat))
sum(a!=0)
dldat=dldat[,a==0]

plot(ts(dldat[,1]))

x=dldat[,-1]
x=as.matrix(x)
y=(dldat[,1])


mod=cv.glmnet(x[1:250,] ,y[1:250],lower.limits=0 )

plot(mod)
plot(glmnet(x ,y,lower.limits=0 ),xvar="lambda")

sum(coef(mod,s=0.001)==0)
sum(coef(mod,s=0.001)!=0)
dim(coef(mod,s=0.001))



#R^2
r2=mod$glmnet.fit$dev.ratio[which(mod$glmnet.fit$lambda >= 0.001)][length(mod$glmnet.fit$dev.ratio[which(mod$glmnet.fit$lambda >= 0.001)])]

1-mod$cvm/var(y)

beta_time=matrix(NA,494,558)
mod$lambda[mod$lambda>=0.001][length(mod$lambda[mod$lambda>=0.001])]



# Predict Y_n ~ X_n
i=1
beta_time[,1]=(coef(mod,s=0.001)[,1])
sum(beta_time[,1]==0)
for( i in 1:558){
  x=as.matrix(dldat[i:(i+249),-1])
  y=as.matrix(dldat[i:(i+249),1])
  mod=cv.glmnet(x ,y,lower.limits=0 )
  beta_time[,i]=(coef(mod,s=0.001)[,1])

}

sum(abs(beta_time[,1]-beta_time[,2]))/494

quantile(beta_time[,1],probs =0.91)
qq=1:557
for (i in 1:553){
  
  
  qq[i]=sum(abs(beta_time[beta_time[,i],i]-beta_time[beta_time[,i],i+5]))/length(beta_time[,i])
  
}

mean(qq)
sum(abs(beta_time[beta_time[,1]<0,1]-beta_time[beta_time[,1]<0,2]))/449


#write.csv(beta_time,"ts_beta.csv")

#Predict Y_n+1 ~ X_n


mod=cv.glmnet(x[11:260,] ,y[12:261],lower.limits=0 )

plot(mod)
mod$lambda.1se
sum(coef(mod,s=0.0001)==0)
sum(coef(mod,s=mod$lambda.min)!=0)
dim(coef(mod,s=0.0001))


a=coef(mod,s=mod$lambda.min)
a@Dimnames[[1]][a@i+1]

beta_time=matrix(NA,494,558)




i=1
beta_time_n1[,1]=(coef(mod,s=0.0001)[,1])
sum(beta_time[,1]==0)
for( i in 1:557){
  x=as.matrix(dldat[i:(i+249),-1])
  y=as.matrix(dldat[i+1:(i+250),1])
  mod=cv.glmnet(x ,y,lower.limits=0 )
  beta_time_n1[,i]=(coef(mod,s=0.001)[,1])
}



# R2

r2 <- mod$glmnet.fit$dev.ratio[which(mod$glmnet.fit$lambda == 0.0001)]


library(readr)
ts_beta <- read_csv("GitHub/NonNegativeLASSO-PortfolioReplication/Data/ts_beta.csv")















# KALMAN FILTER


# New Kalman y ~ 2 series
  
  
dldat <- read.csv("~/GitHub/NonNegativeLASSO-PortfolioReplication/Data/dldat.csv", row.names=1)
x=dldat[1:250,c(2,3)]
x=as.matrix(x)
y=(dldat[1:250,1])

num = length(y)


A = array(1, dim=c(1 ,3 , num))

for (i in 1: num ) {
  A[1,2:3,i]= x[i,]
}
input = rep(1, num )
mu0 = matrix(0, 3, 1)
Sigma0 = diag(c(.1 , .1, 1) , 3)

Linn = function ( para ) {
  Phi = diag (0 ,3)
  Phi [1 ,1]= para [1]  # beta0
  Phi [2,2]= para[2]    # beta1
  Phi [3,3]=para[3]     # beta2
  cQ1= para [4]         
  cQ2= para [5] 
  cQ3= para[6]
  cQ = diag (0 ,3)
  cQ [1 ,1]= cQ1
  cQ [2 ,2]= cQ2
  cQ [3 ,3]= cQ3
  cR = para [7] 
  drift1=para[8]
  drift2=para[9]
  drift3=para[10]
  ups=matrix(0,3,1)
  ups[1,1]=drift1
  ups[2,1]=drift2
  ups[3,1]=drift3
  kf = Kfilter1(num , y , A, mu0 , Sigma0 , Phi ,Ups = ups,Gam = 0, cQ , cR,input)
  return (kf$ like )
}

init.par = c(0.8,.8,.8   ,.1,.1,.1,  .5,  0.1,0.1,0.1  ) 

require(astsa)

est = optim( init.par , Linn ,NULL , method = 'L-BFGS-B', hessian =TRUE)
SE = sqrt ( diag ( solve ( est $ hessian )))
u = cbind ( estimate = est $par , SE)
rownames (u)=c('phi_beta0','phi_beta1','phi_beta2','sig_beta0','sig_beta1','sig_beta2','sig_y','dift0','drift1','drift2')
u



Phi = diag (0 ,3)
Phi [1 ,1]= est$par[1]  # beta0
Phi [2,2]= est$par[2]    # beta1
Phi [3,3]=est$par[3]     # beta2
cQ1= est$par[4]         
cQ2= est$par [5] 
cQ3= est$par[6]
cQ = diag (0 ,3)
cQ [1 ,1]= cQ1
cQ [2 ,2]= cQ2
cQ [3 ,3]= cQ3
cR = est$par [7] 
drift1=est$par[8]
drift2=est$par[9]
drift3=est$par[10]
ups=matrix(0,3,1)
ups[1,1]=drift1
ups[2,1]=drift2
ups[3,1]=drift3
kf = Kfilter1(num , y , A, mu0 , Sigma0 , Phi ,Ups = ups,Gam = 0, cQ , cR,input)


res=1:250
for (i in 1:250){
  res[i]=A[,,i]%*%kf$xp[,,i]
  
}



tsplot(y[200:250])
lines(res[200:250],col=2)





# Kalman extension more regressors


library(glmnet)


dldat <- read.csv("~/GitHub/NonNegativeLASSO-PortfolioReplication/Data/dldat.csv", row.names=1)
a=colSums(is.na(dldat))
sum(a!=0)
dldat=dldat[,a==0]

x=dldat[,-1]
x=as.matrix(x)
y=dldat[,1]
mod=cv.glmnet(x[559:808,] ,y[559:808],lower.limits=0 )
m=coef(mod,s = 0.001)

index=m@i+1
mask=(m[,1]>0)


num_reg=sum(mask)

x=dldat[559:808,mask]
x=x[,c(2:num_reg)]
x=as.matrix(x)
y=(dldat[559:808,1])

num = length(y)


A = array(1, dim=c(1 ,num_reg , num))

for (i in 1: num ) {
  A[1,2:num_reg,i]= x[i,]
}
input = rep(1, num )
mu0 = matrix(0, num_reg, 1)
Sigma0 = diag(0.1 , num_reg)
Sigma0[num_reg,num_reg]=1


Linn = function ( para ) {
  Phi = diag (0 ,num_reg)
  for (i in 1:num_reg){
    Phi[i,i]=para[i]
  }
  
  cQ = diag (0 ,num_reg)
  for (i in 1:num_reg){
    cQ[i,i]=para[i+num_reg]
  }
  cR = para [num_reg*2+1] 
  
  ups=matrix(0,num_reg,1)
  for (i in 1:num_reg){
    ups[i,1]=para[num_reg*2+1+i]
  }
  
  kf = Kfilter1(num , y , A, mu0 , Sigma0 , Phi ,Ups = ups,Gam = 0, cQ , cR,input)
  return ((kf$ like )/100)
}


init.par = c(rep(0.7,num_reg)   ,rep(0.1,num_reg),  .5,  rep(0.1,num_reg)  )

require(astsa)

est = optim( init.par , Linn ,NULL , method = 'L-BFGS-B', hessian =TRUE,control = list(trace=1,REPORT=10))

SE = sqrt ( diag ( solve ( est $ hessian )))
u = cbind ( estimate = est $par , SE)


Phi = diag (0 ,num_reg)
for (i in 1:num_reg){
  Phi[i,i]=est$par[i]
}

cQ = diag (0 ,num_reg)
for (i in 1:num_reg){
  cQ[i,i]=est$par[i+num_reg]
}
cR = est$par [num_reg*2+1] 

ups=matrix(0,num_reg,1)
for (i in 1:num_reg){
  ups[i,1]=est$par[num_reg*2+1+i]
}
kf = Kfilter1(num , y , A, mu0 , Sigma0 , Phi ,Ups = ups,Gam = 0, cQ , cR,input)


res=1:250
for (i in 1:250){
  res[i]=A[,,i]%*%kf$xp[,,i]
  
}



tsplot(y[200:250])
lines(res[200:250],col=2)





# Kalman filter all Lasso-positive series



library(readr)
ts_beta <- read_csv("GitHub/NonNegativeLASSO-PortfolioReplication/Data/ts_beta.csv")
ts_beta=ts_beta[,-1]
t558=(ts_beta[,558])



mask=(t558[,1]>0)
dldat <- read.csv("~/GitHub/NonNegativeLASSO-PortfolioReplication/Data/dldat.csv", row.names=1)
a=colSums(is.na(dldat))
sum(a!=0)
dldat=dldat[,a==0]

x=dldat[559:808,mask]
x=as.matrix(x)
y=(dldat[559:808,1])

num = length(y)


A = array(1, dim=c(1 ,48 , num))

for (i in 1: num ) {
  A[1,2:48,i]= x[i,]
}

input = rep(1, num )
mu0 = matrix(0, 48, 1) #put last available value later
Sigma0 = diag(c(.1 , .1, 1) , 48)



Linn = function ( para ) {
  Phi = diag (0 ,48)
  for (i in 1:48){
    Phi[i,i]=para[i]
  }
  
  cQ = diag (0 ,48)
  for (i in 1:48){
    cQ[i,i]=para[i+48]
  }
  cR = para [97] 
  
  ups=matrix(0,48,1)
  for (i in 1:48){
    ups[i,1]=para[97+i]
  }
  
  kf = Kfilter1(num , y , A, mu0 , Sigma0 , Phi ,Ups = ups,Gam = 0, cQ , cR,input)
  return (kf$ like )
}

init.par = c(rep(0.8,48)   ,rep(0.1,48),  .5,  rep(0.1,48)  ) 
est = optim( init.par , Linn ,NULL , method = 'L-BFGS-B', hessian =TRUE)

SE = sqrt ( diag ( solve ( est $ hessian )))
u = cbind ( estimate = est $par , SE)


Phi = diag (0 ,48)
for (i in 1:48){
  Phi[i,i]=est$par[i]
}

cQ = diag (0 ,48)
for (i in 1:48){
  cQ[i,i]=est$par[i+48]
}
cR = est$par [97] 

ups=matrix(0,48,1)
for (i in 1:48){
  ups[i,1]=est$par[97+i]
}
kf = Kfilter1(num , y , A, mu0 , Sigma0 , Phi ,Ups = ups,Gam = 0, cQ , cR,input)

res=1:250
for (i in 1:250){
  res[i]=A[,,i]%*%kf$xp[,,i]
  
}



lines(y[200:250])
tsplot(res[200:250],col=2)

write.csv(est$par,"param1.csv")

