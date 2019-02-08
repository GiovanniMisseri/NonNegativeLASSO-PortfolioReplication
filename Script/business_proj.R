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

x=dldat[1:250,-1]
x=as.matrix(x)
y=(dldat[1:250,1])

z = x[,c(7,10)]
num = length(y)
A = array(z, dim=c(1 ,2 , num))
input = matrix (1,num ,1)






#b=array(1,dim=c(2,1))


#mu0 = array(0,dim=c(2,1))
#Sigma0 = diag(10,2,2)

#Phi =diag(1,2)
#Ups = (1- phi )%*%b
#Gam=0
#Theta=diag(1,2)
#cQ=array(1,dim=c(1,2))

#cR=1
#S=array(0,dim=c(2,1))


# Function to Calculate Likelihood
Linn = function(para, y.data){ # pass data also
  phi = diag(0,2); phi[1,1]=para[1] ; phi[2,2]=para[2]
  gam = 0
  b = array(0,dim=c(2,1)) ; b[1,1]=para[3]; b[2,1]=para[4]
  Ups = (1- phi )%*%b
  cQ = array(0,dim=c(1,2)); cQ[1,1]=para[5]; cQ[1,2]=para[6]
  cR = para[7]
  kf = Kfilter2(num ,y.data ,A,mu0 , Sigma0 ,phi ,Ups ,gam ,diag(1,2),cQ ,cR ,array(0,dim=c(2,1)), input )
  return (kf$like )
}



require(astsa)
mu0 = array(0,dim=c(2,1)); Sigma0 = diag(10,2,2)
# initial values of the parameters
init.par = c( phi1=.5,phi2=.5 , b1=1 ,b2=1, cQ1=.01,cQ2=.01 , cR1= 0.1)# maximize - loglikelihood
est = optim( init.par , Linn , NULL , y.data =y, method ="L-BFGS-B",hessian =T,control = list(maxit = 6, temp = 2000, trace = TRUE,
                                                                                              REPORT = 1) )
SE = sqrt( diag( solve( est$hessian )))


phi = est$par [1]; alpha = est$par[2]
b = est$par[3]; Ups = (1- phi)*b
cQ = est$par[4]; cR = est$par[5]
round ( cbind ( estimate = est$par , SE), 3)



########################




x=dldat[1:250,-1]
x=as.matrix(x)
y=(dldat[1:250,1])

z = x[,c(1,2)]
num = length(y)
A = array(z, dim=c(1 ,2 , num))
input = matrix (1,num ,1)
b=array(1,dim=c(2,1))


mu0 = array(0,dim=c(2,1))
Sigma0 = diag(10,2,2)

Phi =diag(1,2)
Ups = (1- phi )%*%b
Gam=0
Theta=diag(1,2)
cQ=array(1,dim=c(1,2))

cR=1
S=array(0,dim=c(2,1))


  
  
  Q = t(cQ) %*% cQ
  R = t(cR) %*% cR
  Phi = as.matrix(Phi)
  pdim = nrow(Phi)
  y = as.matrix(y)
  qdim = ncol(y)
  rdim = ncol(as.matrix(input))
  if (max(abs(Ups)) == 0) 
    Ups = matrix(0, pdim, rdim)
  if (max(abs(Gam)) == 0) 
    Gam = matrix(0, qdim, rdim)
  ut = matrix(input, num, rdim)
  xp = array(NA, dim = c(pdim, 1, num))
  Pp = array(NA, dim = c(pdim, pdim, num))
  xf = array(NA, dim = c(pdim, 1, num))
  Pf = array(NA, dim = c(pdim, pdim, num))
  Gain = array(NA, dim = c(pdim, qdim, num))
  innov = array(NA, dim = c(qdim, 1, num))
  sig = array(NA, dim = c(qdim, qdim, num))
  like = 0
  xp[, , 1] = Phi %*% mu0 + Ups %*% as.matrix(ut[1, ], rdim)
  Pp[, , 1] = Phi %*% Sigma0 %*% t(Phi) + Theta %*% Q %*% 
    t(Theta)
  for (i in 1:num) {
    B = matrix(A[, , i], nrow = qdim, ncol = pdim)
    innov[, , i] = y[i, ] - B %*% xp[, , i] - Gam %*% as.matrix(ut[i, 
                                                                   ], rdim)
    sigma = B %*% Pp[, , i] %*% t(B) + R
    sigma = (t(sigma) + sigma)/2
    sig[, , i] = sigma
    siginv = solve(sigma)
    Gain[, , i] = (Phi %*% Pp[, , i] %*% t(B) + Theta %*% 
                     S) %*% siginv
    K = as.matrix(Gain[, , i], nrow = qdim, ncol = pdim)
    xf[, , i] = xp[, , i] + Pp[, , i] %*% t(B) %*% siginv %*% 
      innov[, , i]
    Pf[, , i] = Pp[, , i] - Pp[, , i] %*% t(B) %*% siginv %*% 
      B %*% Pp[, , i]
    sigma = matrix(sigma, nrow = qdim, ncol = qdim)
    like = like + log(det(sigma)) + t(innov[, , i]) %*% 
      siginv %*% innov[, , i]
    if (i == num) 
      break
    xp[, , i + 1] = Phi %*% xp[, , i] + Ups %*% as.matrix(ut[i + 
                                                               1, ], rdim) + K %*% innov[, , i]
    Pp[, , i + 1] = Phi %*% Pp[, , i] %*% t(Phi) + Theta %*% 
      Q %*% t(Theta) - K %*% sig[, , i] %*% t(K)
  }
  like = 0.5 * like
  list(xp = xp, Pp = Pp, xf = xf, Pf = Pf, K = Gain, like = like, 
       innov = innov, sig = sig)












