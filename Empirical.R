source("NetworkVectorAutoregression.R")
source("DebiasLasso.R")
library(rugarch)
library(statnet)
library(openxlsx)
library(rqPen)
library(quantreg)
library(visNetwork)
source("OGA.R")
source("HDIC.R")
source("Trim.R")
source("threeS_loc.R")
source("threeS.R")
####################################################### Crisis
data=read.csv("garch_return.csv",header=TRUE, row.names=1)
p=71
n=36

X=as.matrix(data[((11):(11+n-1)),],nrow = n)
y=as.matrix(data[((12):(12+n-1)),],nrow = n)
Y=as.matrix(data[11:(12+n-1),])
########################################### debiased lasso
A_hat=matrix(0,ncol = p,nrow = p)
pval_hat=matrix(1,ncol = p,nrow = p)

M=M.proxy(X)
for(j in 1:p){
  
  th = SSLasso(X,y[,j],M,verbose = T,alpha = 0.05, lambda = sqrt(log(p)/n))
  
  rev=which(th$coef!=0)
  A_hat[j,]=th$unb
  pval_hat[j,]=th$pvals
}

adjacency=FDR(pval_hat,0.2)$Adjacency

####################################### adaptive lasso
A.VAR.ada=matrix(0,ncol = p,nrow = p)

for(j in 1:p){
  yj=y[,j]
  th2=Adalasso(X,yj,intercept = F,method = "lasso",lambda=1*seq(sqrt(2*log(p)/n),1.3*sqrt(2*log(p)/n),length.out = 100))$AdaLasso
  A.VAR.ada[j,]=th2
}

############################################# three stage
A.three=matrix(0,ncol = p,nrow = p)

for(j in 1:p){
  yj=y[,j]
  th3=threeS_loc(yj,X,2,5,3.01)$trim
  A.three[j,th3]=1
}

###################################### nar model
nva1=NVA(Y,adjacency)
nva2=NVA(Y,A.VAR.ada!=0)
nva3=NVA(Y,A.three)

beta0.1=nva1$Coefficients[1,1]
beta1.1=nva1$Coefficients[2,1]
beta2.1=nva1$Coefficients[3,1]

beta0.2=nva2$Coefficients[1,1]
beta1.2=nva2$Coefficients[2,1]
beta2.2=nva2$Coefficients[3,1]


beta0.3=nva3$Coefficients[1,1]
beta1.3=nva3$Coefficients[2,1]
beta2.3=nva3$Coefficients[3,1]
########################################### intervention


w.1=Matrix.G(adjacency,beta1 = beta1.1,beta2 = beta2.1)$w
G.1=Matrix.G(adjacency,beta1 = beta1.1,beta2 = beta2.1)$G

nu.1=Influential.power(w.1,G.1)$nu

###########################################
w.2=Matrix.G(A.VAR.ada!=0,beta1 = beta1.2,beta2 = beta2.2)$w
G.2=Matrix.G(A.VAR.ada!=0,beta1 = beta1.2,beta2 = beta2.2)$G

nu.2=Influential.power(w.2,G.2)$nu


#############################################
w.3=Matrix.G(A.three,beta1 = beta1.3,beta2 = beta2.3)$w
G.3=Matrix.G(A.three,beta1 = beta1.3,beta2 = beta2.3)$G

nu.3=Influential.power(w.3,G.3)$nu


############################################# Top 10
names(data)[order(nu.3,decreasing = T)][1:10]


sort(nu.3,decreasing = T)[1:10]



############################################ plot network DL

net=network(adjacency,matrix.type="adjacency")
edgelist=as.matrix(net,matrix.type = "edgelist")
edges=as.data.frame(edgelist);colnames(edges)=c("from","to")
edges$color=c(rep("blue",25),rep("red",17),rep("blue",22),rep("red",25),rep("blue",23))
nodes <- data.frame(id = 1:71,label = names(data),value=rep(1,71), size=26,color=c(rep("skyblue",24),"red",rep("skyblue",22),"red",rep("skyblue",23)))

visNetwork(nodes, edges,label=T)
############################################################## plot network ADA
net=network(A.VAR.ada!=0,matrix.type="adjacency")
edgelist=as.matrix(net,matrix.type = "edgelist")
edges=as.data.frame(edgelist);colnames(edges)=c("from","to")
edges$color=c(rep("blue",19),rep("red",8),rep("blue",18),rep("red",26),rep("blue",17))

nodes <- data.frame(id = 1:71,label = names(data),value=rep(1,71), size=26,color=c(rep("skyblue",24),"red",rep("skyblue",22),"red",rep("skyblue",23)))


visNetwork(nodes, edges,label=T)
################################################### plot network three
net=network(A.three,matrix.type="adjacency")
edgelist=as.matrix(net,matrix.type = "edgelist")
edges=as.data.frame(edgelist);colnames(edges)=c("from","to")
edges$color=c(rep("blue",13),rep("red",10),rep("blue",13),rep("red",26),rep("blue",10))

nodes <- data.frame(id = 1:71,label = names(data),value=rep(1,71), size=26,color=c(rep("skyblue",24),"red",rep("skyblue",22),"red",rep("skyblue",23)))

visNetwork(nodes, edges,label=T)

apply(A.three,2,sum)


################################################### quiet

data=read.csv("garch_return.csv",header=TRUE, row.names=1)
n=36
X=as.matrix(data[((nrow(data)-36):(nrow(data)-1)),])
y=as.matrix(data[((nrow(data)-35):nrow(data)),])
Y=as.matrix(data[((nrow(data)-36):(nrow(data))),])

n=dim(y)[1]
p=dim(y)[2]

A_hat=matrix(0,ncol = p,nrow = p)
pval_hat=matrix(1,ncol = p,nrow = p)

######################################################
M=M.proxy(X)
for(j in 1:p){
  
  th = SSLasso(X,y[,j],M,verbose = T,alpha = 0.05, lambda = sqrt(log(p)/n))
  
  rev=which(th$coef!=0)
  # OLS=summary(lm(y[,j]~X[,rev]-1))
  A_hat[j,]=th$unb
  pval_hat[j,]=th$pvals
  
}

adjacency=FDR(pval_hat,0.2)$Adjacency


sum(adjacency)

#######################################
A.VAR.ada=matrix(0,ncol = p,nrow = p)

for(j in 1:p){
  yj=y[,j]
  th2=Adalasso(X,yj,intercept = F,method = "lasso",lambda=0.65*seq(sqrt(2*log(p)/n),1.2*sqrt(2*log(p)/n),length.out = 100))$AdaLasso
  A.VAR.ada[j,]=th2
}

#############################################
A.three=matrix(0,ncol = p,nrow = p)

for(j in 1:p){
  yj=y[,j]
  th3=threeS_loc(yj,X,2,5,3.01)$trim
  A.three[j,th3]=1
}


######################################
net=network(adjacency,matrix.type="adjacency")
edgelist=as.matrix(net,matrix.type = "edgelist")
edges=as.data.frame(edgelist);colnames(edges)=c("from","to")
edges$color=rep("blue",nrow(edges))
edges$color[which(edges$to==31|edges$to==40|edges$to==57)]="red"
nodes <- data.frame(id = 1:71,label = names(data),value=rep(1,71), size=26)
nodes$color=rep("skyblue",nrow(nodes));nodes$color[c(31,40,57)]="red"

visNetwork(nodes, edges,label=T)
##############################################################
net=network(A.VAR.ada!=0,matrix.type="adjacency")
edgelist=as.matrix(net,matrix.type = "edgelist")
edges=as.data.frame(edgelist);colnames(edges)=c("from","to")
edges$color=rep("blue",nrow(edges))

edges$color[which(edges$to==31|edges$to==42|edges$to==57)]="red"
nodes <- data.frame(id = 1:71,label = names(data),value=rep(1,71), size=26)
nodes$color=rep("skyblue",nrow(nodes));nodes$color[c(31,42,57)]="red"

visNetwork(nodes, edges,label=T)
###################################################
net=network(A.three,matrix.type="adjacency")
edgelist=as.matrix(net,matrix.type = "edgelist")
edges=as.data.frame(edgelist);colnames(edges)=c("from","to")
edges$color=rep("blue",nrow(edges))

edges$color[which(edges$to==31|edges$to==42|edges$to==57)]="red"
nodes <- data.frame(id = 1:71,label = names(data),value=rep(1,71), size=26)
nodes$color=rep("skyblue",nrow(nodes));nodes$color[c(31,42,57)]="red"
visNetwork(nodes, edges,label=T)


######################################
nva1=NVA(Y,adjacency)
nva2=NVA(Y,A.VAR.ada!=0)
nva3=NVA(Y,A.three)

beta0.1=nva1$Coefficients[1,1]
beta1.1=nva1$Coefficients[2,1]
beta2.1=nva1$Coefficients[3,1]

beta0.2=nva2$Coefficients[1,1]
beta1.2=nva2$Coefficients[2,1]
beta2.2=nva2$Coefficients[3,1]


beta0.3=nva3$Coefficients[1,1]
beta1.3=nva3$Coefficients[2,1]
beta2.3=nva3$Coefficients[3,1]
###########################################


w.1=Matrix.G(adjacency,beta1 = beta1.1,beta2 = beta2.1)$w
G.1=Matrix.G(adjacency,beta1 = beta1.1,beta2 = beta2.1)$G

nu.1=Influential.power(w.1,G.1)$nu



#######################################
w.2=Matrix.G(A.VAR.ada!=0,beta1 = beta1.2,beta2 = beta2.2)$w
G.2=Matrix.G(A.VAR.ada!=0,beta1 = beta1.2,beta2 = beta2.2)$G

nu.2=Influential.power(w.2,G.2)$nu


#############################################
w.3=Matrix.G(A.three,beta1 = beta1.3,beta2 = beta2.3)$w
G.3=Matrix.G(A.three,beta1 = beta1.3,beta2 = beta2.3)$G

nu.3=Influential.power(w.3,G.3)$nu


#############################################
names(data)[order(nu.3,decreasing = T)][1:10]

sort(nu.1,decreasing = T)[1:10]



######################################### Quantile Regression
rawdata=read.csv("monthly.return.csv",header=TRUE, row.names=1)
data=as.data.frame(apply(rawdata,2,function(x)scale(x,scale = F)))
rownames(data)=rownames(rawdata)
n=36
X=as.matrix(data[((nrow(data)-36):(nrow(data)-1)),])
y=as.matrix(data[((nrow(data)-35):nrow(data)),])
Y=as.matrix(data[((nrow(data)-36):(nrow(data))),])

n=dim(y)[1]
p=dim(y)[2]
dim(Y)

tau=0.95
A_hat=matrix(0,ncol = p,nrow = p)
pval_hat=matrix(1,ncol = p,nrow = p)
lambda=LAMBDA(x = X,tau = tau)
# c=Backtest(X = X,y = y[,j],lambda = lambda,tau=tau,from = 0.1,to = 2,length = 50)
for(j in 1:p){
  c=Backtest(X = X,y = y[,j],lambda = lambda,tau=tau,from = 0.1,to = 2,length = 50)
  lasso.QR=QR_lasso(x = X,y = y[,j],tau = tau,lambda = c*lambda,intercept = F)$coefficients
  if (sum(lasso.QR^2)==0){
    next
  }
  relevant=as.numeric(which(lasso.QR!=0))
  postlasso=rq(y[,j]~X[,relevant]-1,tau = tau)
  pv=summary.rq(postlasso,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
  A_hat[j,relevant]=postlasso$coefficients
  pval_hat[j,relevant]=pv
  
}
Adjacency=FDR(pval_hat,0.2)$Adjacency

Qnva=QR_NVA(Y,Adjacency,tau=tau);Qnva
beta1=Qnva$Coefficients[2,1]
beta2=Qnva$Coefficients[3,1]

Matrix.G(Adjacency,beta1 = beta1,beta2 = beta2)

w=Matrix.G(Adjacency,beta1 = beta1,beta2 = beta2)$w
G=Matrix.G(Adjacency,beta1 = beta1,beta2 = beta2)$G

nu=Influential.power(w,G)$nu
# nu.app=Influential.power(w,G)$approx.nu
plot(nu.app,nu)


order(nu,decreasing = T)
# order(nu.app,decreasing = T)


sort(nu,decreasing = T)[1:10]
# sort(nu.app,decreasing = T)[1:10]

names(data)[order(nu,decreasing = T)][1:10]
# names(data)[order(nu.app,decreasing = T)][1:10]


a=matrix(sort(nu,decreasing = T)[1:10],1)

xtable(a,digits = 4)
names(data)[order(nu,decreasing = T)][1:10]




apply(Adjacency,2,sum)
order(apply(Adjacency,2,sum),decreasing = T)

#################################################

tau=0.9
A_hat=matrix(0,ncol = p,nrow = p)
pval_hat=matrix(1,ncol = p,nrow = p)
lambda=LAMBDA(x = X,tau = tau)
# c=Backtest(X = X,y = y[,j],lambda = lambda,tau=tau,from = 0.1,to = 2,length = 50)
for(j in 1:p){
  c=Backtest(X = X,y = y[,j],lambda = lambda,tau=tau,from = 0.1,to = 2,length = 50)
  lasso.QR=QR_lasso(x = X,y = y[,j],tau = tau,lambda = c*lambda,intercept = F)$coefficients
  if (sum(lasso.QR^2)==0){
    next
  }
  relevant=as.numeric(which(lasso.QR!=0))
  postlasso=rq(y[,j]~X[,relevant]-1,tau = tau)
  pv=summary.rq(postlasso,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
  A_hat[j,relevant]=postlasso$coefficients
  pval_hat[j,relevant]=pv
  
}
Adjacency=FDR(pval_hat,0.2)$Adjacency
apply(Adjacency,2,sum)
order(apply(Adjacency,2,sum),decreasing = T)



W=list()
W_pv=list()
W_adj=list()

q=seq(0.05,0.95,by = 0.02)
A_hat=matrix(0,ncol = p,nrow = p)
pval_hat=matrix(1,ncol = p,nrow = p)


for(t in 1:length(q)){
  lambda=LAMBDA(x = X,tau = q[t])
  A_hat=matrix(0,ncol = p,nrow = p)
  pval_hat=matrix(1,ncol = p,nrow = p)
  
  for(j in 1:p){
    c=Backtest(X = X,y = y[,j],lambda = lambda,tau=q[t],from = 0.1,to = 2,length = 50)
    lasso.QR=QR_lasso(x = X,y = y[,j],tau = q[t],lambda = 0.8755102*lambda,intercept = F)$coefficients
    if (sum(lasso.QR^2)==0){
      next
    }
    relevant=as.numeric(which(lasso.QR!=0))
    postlasso=rq(y[,j]~X[,relevant]-1,tau = q[t])
    pv=summary.rq(postlasso,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
    A_hat[j,relevant]=postlasso$coefficients
    pval_hat[j,relevant]=pv
    
  }
  Adjacency=FDR(pval_hat,0.2)$Adjacency
  W[[t]]=A_hat
  W_pv[[t]]=pval_hat
  W_adj[[t]]=Adjacency
}
################################################################################
