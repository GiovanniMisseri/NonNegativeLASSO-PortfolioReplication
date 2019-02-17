# see https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html to learn more
require(glmnet)
dldat <- read.csv("~/GitHub/NonNegativeLASSO-PortfolioReplication/Data/dldat.csv", row.names=1)
a=colSums(is.na(dldat))
sum(a!=0)
dldat=dldat[,a==0]
plot(ts(dldat[,1]))
x=dldat[,-1]
x=as.matrix(x)
y=(dldat[,1])
library(readr)
ts_beta <- read_csv("GitHub/NonNegativeLASSO-PortfolioReplication/Data/ts_beta.csv")
ts_beta=ts_beta[,-1]
require(nnls)
library(readr)
require(astsa)
require(tictoc)









##  Non negative  Lasso Portfolio replication   ##



#example on first window
mod=cv.glmnet(x[1:250,] ,y[1:250],lower.limits=0 ) #non negative lasso on first window

plot(mod)
abline(v=-6.907755,lty=2,col=4)
title("Lasso accuracy given lambda",line = 2.5)
plot(glmnet(x ,y,lower.limits=0 ),xvar="lambda")

sum(coef(mod,s=0.001)==0)
sum(coef(mod,s=0.001)!=0)




beta_time=matrix(NA,494,558)


# Predict Y_n ~ X_n
#cycle to get all the beta-lasso for all the windows (mooving window)
i=1
beta_time[,1]=(coef(mod,s=0.001)[,1])
sum(beta_time[,1]==0)
for( i in 1:558){
  x1=as.matrix(dldat[i:(i+249),-1])
  y1=as.matrix(dldat[i:(i+249),1])
  mod=cv.glmnet(x1 ,y1,lower.limits=0 )
  beta_time[,i]=(coef(mod,s=0.001)[,1])

}

# mean absolute change

mean_abs=1:250
for (j in 1:250){
  
  step=j
  qq=1:(558-step)
  for (i in 1:(558-step)){ 
    qq[i]=sum((beta_time[,i]-beta_time[,(i+step)]))/length(beta_time[,i])
    
  }
  
  mean_abs[j]=mean(abs(qq))

}

plot(mean_abs) # as the lag grows the beta change, due to time dependence of y
#write.csv(beta_time,"ts_beta.csv")



mean(colSums(ts_beta>0)) #on average 51 assets selected with lasso (lambda=0.001)





## Non negative OLS

#example of the model in the first window
mask=ts_beta[,1]>0 #mask to take only the variables lasso gives positive value
mask=mask[-1]

mod=nnls(cbind(rep(1,250),x[1:250,mask]),y[1:250])#non negative ols

#plot the real time series agains fitted OLS and fitted lasso

tsplot(y[200:250],lwd=1)
lines(mod$fitted[200:250],col=2,lwd=1)
pred_ols=mod$fitted[200:250]

co=as.matrix(ts_beta[c(TRUE,mask),1])
lines(cbind(rep(1,length(200:250)),x[200:250,mask])%*%co,col=4)
pred_lasso=cbind(rep(1,length(200:250)),x[200:250,mask])%*%co
# black=original series
#red=ols fit
#blue=lasso fit

#as we can notice the result is similar, ols seems to approximate sligtly better y


#cycle to get R2 and absolute error for each window in the dataset
R2=matrix(NA,nrow = 558,ncol = 2)
Absol=matrix(NA,nrow = 558,ncol = 2)

i=1
for (i in 1:558){
  mask=ts_beta[,i]>0
  mask=mask[-1]
  mod=nnls(cbind(rep(1,250),x[i:(i+249),mask]),y[i:(i+249)])
  
  pred_ols=mod$fitted
  
  co=as.matrix(ts_beta[c(TRUE,mask),i])
  pred_lasso=cbind(rep(1,250),x[i:(i+249),mask])%*%co
  
  R2[i,1]=1-sum((y[i:(i+249)]-pred_ols)^2)/sum((y[i:(i+249)]-mean(y[i:(i+249)]))^2)
  R2[i,2]=1-sum((y[i:(i+249)]-pred_lasso)^2)/sum((y[i:(i+249)]-mean(y[i:(i+249)]))^2)
  
  a=(abs(y[i:(i+249)]-pred_ols)/abs(y[i:(i+249)]))
  Absol[i,1]=mean(a[a!="Inf" & a!="NaN"])
  a=(abs(y[i:(i+249)]-pred_lasso)/abs(y[i:(i+249)]))
  Absol[i,2]=mean(a[a!="Inf" & a!="NaN"])
}

sum(R2[,2]>R2[,1]) #the R squared using OLS is always higher, similar result for asbolute error
sum(Absol[,2]<Absol[,1]) # OLS seems to have better performances but difference not as much drastic as before



## Portfolio simulation passive approach ##
mask=ts_beta[,1]>0
mask=mask[-1]
x=as.matrix(x)

mod=nnls(cbind(rep(1,250),x[1:250,mask]),y[1:250]) #Non negative ols on the lasso-selected assets

co=as.matrix(ts_beta[c(TRUE,mask),1])
pred_lasso=cbind(rep(1,length(1:808)),x[1:808,mask])%*%co # lasso prediction for the whole period


oo=exp(cumsum(y[251:808])) # correspond to percentage value of the sp500 index with respect time 250
oo[558] # final sp500 value (111% initial value)

tsplot(oo,ylim=c(0.88,1.192),lwd=2,col="blue",main = "OLS Portfolio")
lines(exp(cumsum(pred_lasso[251:808])),col="gray40",lty=2)


co=mod$x
pred_ols=cbind(rep(1,length(1:808)),x[1:808,mask])%*%co
lines(exp(cumsum(pred_ols[251:808])),col="orange",lty=1,lwd=2)


lines(exp(cumsum(pred_lasso[251:808])),col="gray40",lty=2)

legend( "topleft",legend = c("initial buget","Sp500 index","OLS portfolio","Lasso portfolio"),col = c(1,"blue","orange","gray40"),lty = c(1,1,1,2),lwd=c(2,2,2,1))
abline(h=1)



# elastic net

mod=glmnet(x[1:250,] ,y[1:250],lower.limits=0,alpha = 0.5 )
ss=0.0025
sum(coef(mod,s=ss)!=0)
a=coef(mod,s=ss)[,1]
pred_elas=cbind(rep(1,length(1:808)),x[1:808,])%*%a
require(astsa)
ap=exp(cumsum(pred_elas[251:808]))
ap[558]      
oo=exp(cumsum(y[251:808]))
oo[558]
tsplot(oo,ylim=c(0.88,1.19),lwd=2,col="blue",ylab="",main="Elastic net Portfolio")
lines(exp(cumsum(pred_elas[251:808])),col="orange",lwd=2)
legend( "topleft",legend = c("initial buget","Sp500 index","Elastic portfolio","Lasso portfolio"),col = c(1,"blue","orange","gray40"),lty = c(1,1,1,2),lwd=c(2,2,2,1))
abline(h=1)
# 
mod=glmnet(x[1:250,] ,y[1:250],lower.limits=0,alpha = 1)
sum(coef(mod,s=0.001)!=0)
a=coef(mod,s=0.001)[,1]
pred_elas=cbind(rep(1,length(1:808)),x[1:808,])%*%a

lines(exp(cumsum(pred_elas[251:808])),col="gray40",lwd=1,lty=2)




# KALMAN FILTER


# New Kalman y ~ 2 series

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
  return (kf$ like/100 )
}

init.par = c(0.8,.8,.8   ,.1,.1,.1,  .5,  0.1,0.1,0.1  ) 

tic("optim")
est = optim( init.par , Linn ,NULL , method = 'L-BFGS-B', hessian =T)
toc()

tic("nlminb")
est1= nlminb(init.par,Linn,control = list(trace=5))
toc()
SE = sqrt ( diag ( solve ( est $ hessian )))
u = cbind ( estimate = est $par , SE)
rownames (u)=c('phi_beta0','phi_beta1','phi_beta2','sig_beta0','sig_beta1','sig_beta2','sig_y','dift0','drift1','drift2')
u


est=est1
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




## Dynamic regression model ##

mod=cv.glmnet(x[1:250,] ,y[1:250],lower.limits=0 ) #changing lambda change the number of coefficient, with lambda=0.001 the dynamic regression model is done on the same variable lasso has been done
m=coef(mod,s = 0.001)

index=m@i+1
m@Dimnames[[1]][m@i+1] #name of the considered assets

mask=(m[,1]>0)

num_reg=sum(mask)

x=dldat[1:808,mask]
x=x[,c(2:num_reg)]
x=as.matrix(x)
y=(dldat[1:808,1])

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
  return ((kf$ like ) )
}


init.par = c(rep(0.75,num_reg)   ,rep(0.1,num_reg),  .5,  rep(0.1,num_reg)  )


Linn(init.par)
require(astsa)
require(tictoc)

tic("nlminb")
est=nlminb(init.par,Linn,control = list(trace=10)) #nlminb doesn't give the hessian so we dont have estimate for the variance, but with the full model optim took 6 hours to run, nlminb took 1.2 hour
toc()


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


res=1:808
for (i in 1:808){
  res[i]=A[,,i]%*%kf$xp[,,i]
  
}



tsplot(res[1:808],lwd=1,col=2,main="Kalman filter prediction",ylab = "")
lines(y[1:808],col=1,lwd=2)
legend("topleft",legend = c("Daily log-return","Kalman prediction"),col=c(1,2),lty=c(1,1),lwd=c(2,1))







#Dynamic mod
sum(rowSums(ts_beta)==0)
post=ts_beta[rowSums(ts_beta)!=0,]
tsplot(as.numeric(post[2,]),col=2,main="Beta Lasso time series", ylim=c(0,0.091),ylab = "")
for (i in 3:558){
  lines(as.numeric(post[i,]),col=i)
}

u=1:558
for (i in 1:558){
  u[i]=quantile(t(post[,i]),probs = 0.95)
  
}
mean(u)
0.05*216





# revenue of dynamic regression model - plot

oo=exp(cumsum(y[251:808])) 
oo[808]
tsplot(oo,ylim=c(0.8,1.38),lwd=2,col="blue",ylab="",main="Dynamic Portfolio - 47 assets")
lines(exp(cumsum(res[251:808])),col="orange",lwd=2)
legend( "topleft",legend = c("initial buget","Sp500 index","Dynamic portfolio","Ols portfolio"),col = c(1,"blue","orange","gray40"),lty = c(1,1,1,2),lwd=c(2,2,2,1))
abline(h=1)



