library(readxl)
dat <- read_excel("C:/Users/luigimissri/Desktop/FILE DA SALVARE/Data Science/II anno/Business economics and financial data/Project/SP500_10feb2017_components.xlsx", 
                  sheet = "SPX Components", skip = 2)
library(readxl)
dat1 <- read_excel("C:/Users/luigimissri/Desktop/FILE DA SALVARE/Data Science/II anno/Business economics and financial data/Project/SP500_10feb2017_components.xlsx", 
                   sheet = "SPX Index", skip = 2)
names=read_excel("C:/Users/luigimissri/Desktop/FILE DA SALVARE/Data Science/II anno/Business economics and financial data/Project/SP500_10feb2017_components.xlsx", 
                 sheet = "Names",col_names = F)

tot=merge(dat1,dat,by = "CURRENCY")
colnames(tot)=names[1,]

# write.csv(tot,"dat_tot.csv")


dat=dat_tot[10700:11508,c(-1,-2)]
date=tot$Date[10701:11508]
lodat=log(dat)
dldat=lodat[2:809,]-lodat[1:808,]
#write.csv(date,"date.csv")
#write.csv(dldat,"dldat.csv")

u=1:506
u[colSums(is.na(dldat))>0]
