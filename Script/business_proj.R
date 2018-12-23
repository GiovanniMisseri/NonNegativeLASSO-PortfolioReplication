# see https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html to learn more
dldat <- read.csv("~/GitHub/NonNegativeLASSO-PortfolioReplication/Data/dldat.csv", row.names=1)
require(glmnet)
a=colSums(is.na(dldat))
sum(a!=0)
dldat=dldat[,a==0]


x=dldat[,-1]
x=as.matrix(x)
y=(dldat[,1])


mod=cv.glmnet(x ,y,lower.limits=0 )

plot(mod)
mod$lambda.1se
sum(coef(mod,s=0.0001)==0)
sum(coef(mod,s=0.0001)!=0)
dim(coef(mod,s=0.0001))

beta_time=matrix(NA,494,807)

i=1
beta_time[,1]=(coef(mod,s=0.001)[,1])
sum(beta_time[,1]==0)
for( i in 1:558){
  x=as.matrix(dldat[i:(i+249),-1])
  y=as.matrix(dldat[i:(i+249),1])
  mod=cv.glmnet(x ,y,lower.limits=0 )
  beta_time[,i]=(coef(mod,s=0.001)[,1])

}

#write.csv(beta_time,"ts_beta.csv")
