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

