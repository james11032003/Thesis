source("NetworkVectorAutoregression.R")
source("DebiasLasso.R")
library(BigVAR)
library(gplots)
library(statnet)


source("OGA.R")
source("HDIC.R")
source("Trim.R")
source("threeS_loc.R")
source("threeS.R")

n=500
p=20


set.seed(88)


du=runif(p,-0.35,-0.3)

dia=du

A=diag(dia)

C=combn(p,2)

for(i in 1:ncol(C)){
  cu1=C[1,i]
  cu2=C[2,i]
  
  A[cu1,cu2]=rbinom(1,1,0.06)*(2*rbinom(1,1,0.8)-1)*runif(1,0.4,0.45)
  
  cl1=C[2,i]
  cl2=C[1,i]
  
  A[cl1,cl2]=rbinom(1,1,0.06)*(2*rbinom(1,1,0.8)-1)*runif(1,0.4,0.45)
}

eig=eigen(A)$values
max(abs(eig))


#####################################################
r=500


rbeta0.1=numeric(r)
rbeta1.1=numeric(r)
rbeta2.1=numeric(r)

rbeta0.2=numeric(r)
rbeta1.2=numeric(r)
rbeta2.2=numeric(r)

rbeta0.3=numeric(r)
rbeta1.3=numeric(r)
rbeta2.3=numeric(r)


precition1=numeric(r)
precition2=numeric(r)
precition3=numeric(r)

recall1=numeric(r)
recall2=numeric(r)
recall3=numeric(r)

NU.1=matrix(0,p,r)
APP.NU.1=matrix(0,p,r)

NU.2=matrix(0,p,r)
APP.NU.2=matrix(0,p,r)

NU.3=matrix(0,p,r)
APP.NU.3=matrix(0,p,r)


final.VAR=matrix(0,p,p)
final.ADA=matrix(0,p,p)
final.three=matrix(0,p,p)


set.seed(88)
for(w in 1:r){
  
  X=NULL
  y=NULL
  x0=matrix(rnorm(p),ncol = 1)
  y1=A%*%x0+rnorm(p)
  y=rbind(y,t(y1))
  X=rbind(X,t(x0))
  for(i in 2:1000){
    x0=y[i-1,]
    X=rbind(X,x0)
    y1=A%*%x0+rnorm(p)
    y=rbind(y,t(y1))
  }
  data=matrix(rbind(X,y[1000,]),ncol=p)
  

  data=data[(nrow(data)-((n+1+11))+1):(nrow(data)),]
  
  
  Y.train=data[1:(n+1),]
  
  X.train=Y.train[1:n,]
  y.train=Y.train[2:(n+1),]
  

  
  ################################################ VAR debiased
  A.VAR=matrix(0,ncol = p,nrow = p)
  pval.VAR=matrix(0,ncol = p,nrow = p)
  
  
  M=M.proxy(X.train)
  for(j in 1:p){
    yj=y.train[,j]
    th = SSLasso(X.train,yj,M,verbose = T,alpha = 0.05
                 ,lambda = 10*sqrt(2*log(p)/n)
    )
    A.VAR[j,]=th$unb
    pval.VAR[j,]=th$pvals
  }
  
  adjacency = FDR(A = pval.VAR,alpha = 0.01)$Adjacency
  
  final.VAR=final.VAR+adjacency
  
  Table(adjacency,A)
  
  precition1[w]=Table(adjacency,A)[1]
  recall1[w]=Table(adjacency,A)[2]
  

  
  ############################################# adaptive lasso
  A.VAR.ada=matrix(0,ncol = p,nrow = p)
  
  for(j in 1:p){
    yj=y.train[,j]
    th2=Adalasso(X.train,yj,intercept = F,method = "lasso",lambda=seq(1.3*sqrt(2*log(p)/n),1.6*sqrt(2*log(p)/n),length.out = 100))$AdaLasso
    A.VAR.ada[j,]=th2
  }
  
  final.ADA=final.ADA+((A.VAR.ada!=0)*1)
  
  Table(A.VAR.ada,A)
  precition2[w]=Table(A.VAR.ada,A)[1]
  recall2[w]=Table(A.VAR.ada,A)[2]

  
  ############################################# threes
  A.three=matrix(0,ncol = p,nrow = p)
  
  for(j in 1:p){
    yj=y.train[,j]
    th3=threeS_loc(yj,X.train,2,5,3.01)$trim
    A.three[j,th3]=1
  }
  
  final.three=final.three+A.three
  
  Table(A.three,A)
  precition3[w]=Table(A.three,A)[1]
  recall3[w]=Table(A.three,A)[2]
  #############################################
  
  NAR1=NVA(Y.train,adjacency)
  
  beta0.1=NAR1$Coefficients[1,1]
  beta1.1=NAR1$Coefficients[2,1]
  beta2.1=NAR1$Coefficients[3,1]
  
  w.1=Matrix.G(adjacency,beta1=  beta1.1,beta2 = beta2.1)$w
  G.1=Matrix.G(adjacency,beta1=  beta1.1,beta2 = beta2.1)$G
  
  NU.1[,w]=Influential.power(w.1,G.1)$nu
  APP.NU.1[,w]=Influential.power(w.1,G.1,beta1 = beta1.1,beta2 = beta2.1)$approx.nu
  
  
  rbeta0.1[w]=beta0.1
  rbeta1.1[w]=beta1.1
  rbeta2.1[w]=beta2.1
  ############################################################
  NAR2=NVA(Y.train,((A.VAR.ada!=0)*1))
  beta0.2=NAR2$Coefficients[1,1]
  beta1.2=NAR2$Coefficients[2,1]
  beta2.2=NAR2$Coefficients[3,1]
  
  w.2=Matrix.G(((A.VAR.ada!=0)*1),beta1=  beta1.2,beta2 = beta2.2)$w
  G.2=Matrix.G(((A.VAR.ada!=0)*1),beta1=  beta1.2,beta2 = beta2.2)$G
  
  NU.2[,w]=Influential.power(w.2,G.2)$nu
  APP.NU.2[,w]=Influential.power(w.2,G.2,beta1 = beta1.2,beta2 = beta2.2)$approx.nu
  
  
  rbeta0.2[w]=beta0.2
  rbeta1.2[w]=beta1.2
  rbeta2.2[w]=beta2.2
  ########################################################
  NAR3=NVA(Y.train,A.three)
  beta0.3=NAR3$Coefficients[1,1]
  beta1.3=NAR3$Coefficients[2,1]
  beta2.3=NAR3$Coefficients[3,1]
  
  w.3=Matrix.G(A.three,beta1=  beta1.3,beta2 = beta2.3)$w
  G.3=Matrix.G(A.three,beta1=  beta1.3,beta2 = beta2.3)$G
  
  NU.3[,w]=Influential.power(w.3,G.3)$nu
  APP.NU.3[,w]=Influential.power(w.3,G.3,beta1 = beta1.3,beta2 = beta2.3)$approx.nu
  
  
  rbeta0.3[w]=beta0.3
  rbeta1.3[w]=beta1.3
  rbeta2.3[w]=beta2.3
  
}

##################################### heap map

heatmap.2((A!=0)*1,p,main = 'True model',scale = "none",Rowv = F,Colv = F,trace='none',dendrogram="none",col=colorpanel(50,"darkgray","black","red"))

heatmap.2(final.VAR,p,p,main = 'Debiased Lasso VAR model',scale = "none",Rowv = F,Colv = F,trace='none',dendrogram="none",col=colorpanel(500,"darkgray","black","red"))
final.VAR2=log(final.VAR+1)

heatmap.2(final.VAR2,p,p,main = 'Debiased Lasso VAR model',scale = "none",Rowv = F,Colv = F,trace='none',dendrogram="none",col=colorpanel(500,"darkgray","black","red"))


heatmap.2(final.ADA,p,p,main = 'Adaptive Lasso VAR model',scale = "none",Rowv = F,Colv = F,trace='none',dendrogram="none",col=colorpanel(500,"darkgray","black","red"))
final.ADA2=log(final.ADA+1)

heatmap.2(final.ADA2,p,p,main = 'Adaptive Lasso VAR model',scale = "none",Rowv = F,Colv = F,trace='none',dendrogram="none",col=colorpanel(500,"darkgray","black","red"))

heatmap.2(final.three,p,p,main = 'OGA+HDIC+TRIM VAR',scale = "none",Rowv = F,Colv = F,trace='none',dendrogram="none",col=colorpanel(500,"darkgray","black","red"))
final.three2=log(final.three+1)

heatmap.2(final.three2,p,p,main = 'OGA+HDIC+TRIM VAR',scale = "none",Rowv = F,Colv = F,trace='none',dendrogram="none",col=colorpanel(500,"darkgray","black","red"))


#################################
a=A-A*diag(dim(A)[2])

mean(a[which(a!=0,arr.ind = T)])

mean(diag(A))

################################# 0.2675212  -0.3231112
mean(rbeta0.1)
mean(rbeta0.2)
mean(rbeta0.3)

mean(rbeta1.1)
mean(rbeta1.2)
mean(rbeta1.3)

mean(rbeta2.1)
mean(rbeta2.2)
mean(rbeta2.3)

mean(precition1)
mean(precition2)
mean(precition3)

mean(recall1)
mean(recall2)
mean(recall3)

nu.1=round(apply(NU.1,1,mean),4)
app.nu.1=round(apply(APP.NU.1,1,mean),4)

nu.2=round(apply(NU.2,1,mean),4)
app.nu.2=round(apply(APP.NU.2,1,mean),4)


nu.3=round(apply(NU.3,1,mean),4)
app.nu.3=round(apply(APP.NU.3,1,mean),4)


