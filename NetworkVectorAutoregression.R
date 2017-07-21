library(Matrix)
library(quantreg)
########################################################################
NVA=function(Y,A,Z=NULL,intercept=T,subset=NULL){
  Y=t(Y)
  if (is.null(subset)){
    if (is.null(Z)){
      n=dim(Y)[1]# # of firms
      t=dim(Y)[2]#time
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,n*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1,ncol=1)
      X[,2]=matrix(Yt_1,ncol = 1)
      Yt =  matrix(Y[,(2:t)],ncol = 1)
      
      
      if(intercept==T){
        beta=matrix(0,3,1)
        pval=matrix(0,3,1)
        conf=matrix(NA,3,2)
        
        OLS=lm(Yt~X)$coefficients
        beta[,1]=OLS
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          pval[c(1,3),1]=summary(lm(Yt~X))$coefficients[,4]
          conf[c(1,3),]=confint(lm(Yt~X))[c(1,3),]
        }else{
          pval[,1]=summary(lm(Yt~X))$coefficients[,4]
          conf=confint(lm(Yt~X))
        }
        
        rownames(beta)=c("beta0","beta1","beta2")
        rownames(pval)=c("beta0","beta1","beta2")
        rownames(conf)=c("beta0","beta1","beta2")
        colnames(conf)=c("2.5 %","97.5 %")
        # beta=matrix(0,2,1)
        # pval=matrix(0,2,1)
        # 
        # OLS=lm(Yt~X[,2])$coefficients
        # beta[,1]=OLS
        # pval[,1]=summary(lm(Yt~X[,2]))$coefficients[,4]
        # rownames(beta)=c("beta0","beta2")
        # rownames(pval)=c("beta0","beta2")
        
      } else{
        beta=matrix(0,2,1)
        pval=matrix(0,2,1)
        conf=matrix(NA,2,2)
        
        
        
        OLS=lm(Yt~X-1)$coefficients
        beta[,1]=OLS
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          pval[2,1]=summary(lm(Yt~X-1))$coefficients[,4]
          conf[2,]=confint(lm(Yt~X-1))[2,]
        }else{
          pval[,1]=summary(lm(Yt~X-1))$coefficients[,4]
          conf=confint(lm(Yt~X-1))
        }
        
        rownames(beta)=c("beta1","beta2")
        rownames(pval)=c("beta1","beta2")
        rownames(conf)=c("beta1","beta2")
        colnames(conf)=c("2.5 %","97.5 %")
      }
      
    }else{
      n=dim(Y)[1] #  # of firms
      t=dim(Y)[2] # time
      p=dim(Z)[2]
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,n*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1,ncol=1)
      X[,2]=matrix(Yt_1,ncol = 1)
      z=Z
      for(i in 1:(t-2)){z=rbind(z,Z)}
      X=cbind(X,z)
      Yt =  matrix(Y[,2:t],ncol = 1)
      if(intercept==T){
        beta=matrix(0,3+p,1)
        pval=matrix(0,3+p,1)
        OLS=lm(Yt~X)$coefficients
        
        beta[,1]=OLS
        
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          pval[c(1,3:(3+p)),1]=summary(lm(Yt~X))$coefficients[,4]
        }else{
          pval[,1]=summary(lm(Yt~X))$coefficients[,4]
        }
        
        
        rownames(beta)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
      } else{
        beta=matrix(0,2+p,1)
        pval=matrix(0,2+p,1)
        OLS=lm(Yt~X-1)$coefficients
        beta[,1]=OLS
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          pval[2:(2+p),1]=summary(lm(Yt~X-1))$coefficients[,4]
        }else{
          pval[,1]=summary(lm(Yt~X-1))$coefficients[,4]
        }
        
        rownames(beta)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
      }
    }
  }else{
    if (is.null(Z)){
      n=dim(Y)[1]# # of firms
      t=dim(Y)[2]#time
      s=length(subset)
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,s*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1[subset,],ncol=1)
      X[,2]=matrix(Yt_1[subset,],ncol = 1)
      Yt =  matrix(Y[,2:t][subset,],ncol = 1)
      
      
      if(intercept==T){
        beta=matrix(0,3,1)
        pval=matrix(0,3,1)
        conf=matrix(NA,3,2)
        
        OLS=lm(Yt~X)$coefficients
        beta[,1]=OLS
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          pval[c(1,3),1]=summary(lm(Yt~X))$coefficients[,4]
          conf[c(1,3),]=confint(lm(Yt~X))[c(1,3),]
        }else{
          pval[,1]=summary(lm(Yt~X))$coefficients[,4]
          conf=confint(lm(Yt~X))
        }
        
        rownames(beta)=c("beta0","beta1","beta2")
        rownames(pval)=c("beta0","beta1","beta2")
        rownames(conf)=c("beta0","beta1","beta2")
        colnames(conf)=c("2.5 %","97.5 %")
        # beta=matrix(0,2,1)
        # pval=matrix(0,2,1)
        # 
        # OLS=lm(Yt~X[,2])$coefficients
        # beta[,1]=OLS
        # pval[,1]=summary(lm(Yt~X[,2]))$coefficients[,4]
        # rownames(beta)=c("beta0","beta2")
        # rownames(pval)=c("beta0","beta2")
        
      } else{
        beta=matrix(0,2,1)
        pval=matrix(0,2,1)
        conf=matrix(NA,2,2)
        
        
        
        OLS=lm(Yt~X-1)$coefficients
        beta[,1]=OLS
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          pval[2,1]=summary(lm(Yt~X-1))$coefficients[,4]
          conf[2,]=confint(lm(Yt~X-1))[2,]
        }else{
          pval[,1]=summary(lm(Yt~X-1))$coefficients[,4]
          conf=confint(lm(Yt~X-1))
        }
        
        rownames(beta)=c("beta1","beta2")
        rownames(pval)=c("beta1","beta2")
        rownames(conf)=c("beta1","beta2")
        colnames(conf)=c("2.5 %","97.5 %")
      }
      
    }else{
      n=dim(Y)[1] #  # of firms
      t=dim(Y)[2] # time
      p=dim(Z)[2]
      s=length(subset)
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,s*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1[subset,],ncol=1)
      X[,2]=matrix(Yt_1[subset,],ncol = 1)
      z=Z[subset,]
      for(i in 1:(t-2)){z=rbind(z,Z[subset,])}
      X=cbind(X,z)
      Yt =  matrix(Y[,2:t][subset,],ncol = 1)
      if(intercept==T){
        beta=matrix(0,3+p,1)
        pval=matrix(0,3+p,1)
        OLS=lm(Yt~X)$coefficients
        
        beta[,1]=OLS
        
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          pval[c(1,3:(3+p)),1]=summary(lm(Yt~X))$coefficients[,4]
        }else{
          pval[,1]=summary(lm(Yt~X))$coefficients[,4]
        }
        
        
        rownames(beta)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
      } else{
        beta=matrix(0,2+p,1)
        pval=matrix(0,2+p,1)
        OLS=lm(Yt~X-1)$coefficients
        beta[,1]=OLS
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          pval[2:(2+p),1]=summary(lm(Yt~X-1))$coefficients[,4]
        }else{
          pval[,1]=summary(lm(Yt~X-1))$coefficients[,4]
        }
        
        rownames(beta)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
      }
    }
    
  }
  beta=round(beta,4)
  pval=round(pval,4)
  returnlist <- list("Coefficients"= beta,
                     "Pvals" = pval,
                     "Conf" = conf
  )
  return(returnlist)
  
}
########################################################################
INVA=function(Y,A,intercept=T){
  t=dim(Y)[1]
  n=dim(Y)[2]
  a=A-A*diag(n)
  w=a
  for(i in 1:n){
    if(sum(w[i,]!=0)){
      w[i,]=w[i,]/sum(w[i,])
    }
  }
  yt_1=Y[1:t-1,]
  yt=Y[2:t,]
  w.yt_1=yt_1%*%t(w)
  
  if(intercept==T){
    
    beta=matrix(0,3,n)
    pval=matrix(0,3,n)
    for(i in 1:n){
      OLS=lm(yt[,i]~w.yt_1[,i]+yt_1[,i])
      beta[,i]=OLS$coefficients
      if(is.na(beta[2,i])){
        pval[2,i]=NA
        pval[c(1,3),i]=summary(OLS)$coefficients[,4]
      }else{
        pval[,i]=summary(OLS)$coefficients[,4]
      }}
    rownames(beta)=c("beta0","beta1","beta2")
    rownames(pval)=c("beta0","beta1","beta2")
  }else{
    beta=matrix(0,2,n)
    pval=matrix(0,2,n)
    for(i in 1:n){
      OLS=lm(yt[,i]~w.yt_1[,i]+yt_1[,i]-1)
      beta[,i]=OLS$coefficients
      if(is.na(beta[1,i])){
        pval[1,i]=NA
        pval[2,i]=summary(OLS)$coefficients[,4]
      }else{
        pval[,i]=summary(OLS)$coefficients[,4]
      }}
    
    rownames(beta)=c("beta1","beta2")
    rownames(pval)=c("beta1","beta2")
  }
  beta=round(beta,4)
  pval=round(pval,4)
  returnlist <- list("Coefficients"= beta,
                     "Pvals" = pval
  )
  return(returnlist)
}
########################################
Pagerank=function(A){
  A=t(A)
  n=dim(A)[1]
  S=apply(A,2,sum)
  P=A
  m=0.15
  for(i in 1:n){
    if(sum(P[,i])!=0){
      P[,i]=P[,i]/S[i]
    }else{P[,i]=1/n}
  }
  # if(rankMatrix(P)!=n){
    M=(1-m)*P+m*matrix(1/n,n,n)
    P=M
  # }
  value=eigen(P)$values
  vector=eigen(P)$vectors
  if(sum((Im(eigen(P)$values))^2)!=0){
    value=Re(value[-which(Im(eigen(P)$values)!=0)])
    vector=vector[,-which(Im(eigen(P)$values)!=0)]
  }
  vector=(abs(Re(eigen(P)$vectors[,which(value==max(value))])))/sum(abs(Re(eigen(P)$vectors[,which(value==max(value))])))
  return(matrix(vector*100,nrow = n))
}
################################################################################
QR_NVA=function(Y,A,Z=NULL,intercept=T,subset=NULL,tau){
  Y=t(Y)
  if (is.null(subset)){
    if (is.null(Z)){
      n=dim(Y)[1]# # of firms
      t=dim(Y)[2]#time
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,n*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1,ncol=1)
      X[,2]=matrix(Yt_1,ncol = 1)
      Yt =  matrix(Y[,2:t],ncol = 1)
      
      
      if(intercept==T){
        beta=matrix(0,3,1)
        pval=matrix(0,3,1)
        conf=matrix(NA,3,2)
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2],tau=tau)
          beta[c(1,3),1]=QR$coefficients
          beta[2,1]=NA
          
        }else{
          QR=rq(Yt~X,tau=tau)
          beta[,1]=QR$coefficients
        }
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[c(1,3),1]=summary[,4]
          conf[c(1,3),1]=(summary[,1]-(qt(0.975,dim(Yt)[1]-3)*summary[,2]))[c(1,3)]
          conf[c(1,3),2](summary[,1]+(qt(0.975,dim(Yt)[1]-3)*summary[,2]))[c(1,3)]
        }else{
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[,1]=summary[,4]
          conf[,1]=(summary[,1]-(qt(0.975,dim(Yt)[1]-3)*summary[,2]))
          conf[,2]=(summary[,1]+(qt(0.975,dim(Yt)[1]-3)*summary[,2]))
        }
        
        rownames(beta)=c("beta0","beta1","beta2")
        rownames(pval)=c("beta0","beta1","beta2")
        rownames(conf)=c("beta0","beta1","beta2")
        colnames(conf)=c("2.5 %","97.5 %")
        # beta=matrix(0,2,1)
        # pval=matrix(0,2,1)
        # 
        # OLS=lm(Yt~X[,2])$coefficients
        # beta[,1]=OLS
        # pval[,1]=summary(lm(Yt~X[,2]))$coefficients[,4]
        # rownames(beta)=c("beta0","beta2")
        # rownames(pval)=c("beta0","beta2")
        
      } else{
        beta=matrix(0,2,1)
        pval=matrix(0,2,1)
        conf=matrix(NA,2,2)
        
        
        
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2]-1,tau=tau)
          beta[2,1]=QR$coefficients
          beta[1,1]=NA
        }else{
          QR=rq(Yt~X-1,tau=tau)
          beta[,1]=QR$coefficients
        }
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[2,1]=summary[,4]
          conf[2,1]=(summary[,2]-(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
          conf[2,2]=(summary[,2]+(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
        }else{
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[,1]=summary[,4]
          conf[,1]=(summary[,2]-(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
          conf[,2]=(summary[,2]+(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
        }
        
        rownames(beta)=c("beta1","beta2")
        rownames(pval)=c("beta1","beta2")
        rownames(conf)=c("beta1","beta2")
        colnames(conf)=c("2.5 %","97.5 %")
      }
      
    }else{
      n=dim(Y)[1] #  # of firms
      t=dim(Y)[2] # time
      p=dim(Z)[2]
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,n*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1,ncol=1)
      X[,2]=matrix(Yt_1,ncol = 1)
      z=Z
      for(i in 1:(t-2)){z=rbind(z,Z)}
      X=cbind(X,z)
      Yt =  matrix(Y[,2:t],ncol = 1)
      if(intercept==T){
        beta=matrix(0,3+p,1)
        pval=matrix(0,3+p,1)
        
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2],tau=tau)
          beta[c(1,3:(3+p)),1]=QR$coefficients
          beta[2,1]=NA
        }else{
          QR=rq(Yt~X,tau=tau)
          beta[,1]=QR$coefficients
        }
        
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          pval[c(1,3:(3+p)),1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }else{
          pval[,1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }
        
        
        rownames(beta)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
      } else{
        beta=matrix(0,2+p,1)
        pval=matrix(0,2+p,1)
        
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2]-1,tau=tau)
          beta[2:(2+p),1]=QR$coefficients
          beta[1,1]=NA
        }else{
          QR=rq(Yt~X-1,tau=tau)
          beta[,1]=QR$coefficients
        }
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          pval[2:(2+p),1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }else{
          pval[,1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }
        
        rownames(beta)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
      }
    }
  }else{
    if (is.null(Z)){
      n=dim(Y)[1]# # of firms
      t=dim(Y)[2]#time
      s=length(subset)
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,s*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1[subset,],ncol=1)
      X[,2]=matrix(Yt_1[subset,],ncol = 1)
      Yt =  matrix(Y[,2:t][subset,],ncol = 1)
      
      
      if(intercept==T){
        beta=matrix(0,3,1)
        pval=matrix(0,3,1)
        conf=matrix(NA,3,2)
        
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2],tau=tau)
          beta[c(1,3),1]=QR$coefficients
          beta[2,1]=NA
        }else{
          QR=rq(Yt~X,tau=tau)
          beta[,1]=QR$coefficients
        }
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[c(1,3),1]=summary[,4]
          conf[c(1,3),1]=(summary[,1]-(qt(0.975,dim(Yt)[1]-3)*summary[,2]))
          conf[c(1,3),2]=(summary[,1]+(qt(0.975,dim(Yt)[1]-3)*summary[,2]))
        }else{
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[,1]=summary[,4]
          conf[,1]=(summary[,1]-(qt(0.975,dim(Yt)[1]-3)*summary[,2]))
          conf[,2]=(summary[,1]+(qt(0.975,dim(Yt)[1]-3)*summary[,2]))
        }
        
        rownames(beta)=c("beta0","beta1","beta2")
        rownames(pval)=c("beta0","beta1","beta2")
        rownames(conf)=c("beta0","beta1","beta2")
        colnames(conf)=c("2.5 %","97.5 %")
        # beta=matrix(0,2,1)
        # pval=matrix(0,2,1)
        # 
        # OLS=lm(Yt~X[,2])$coefficients
        # beta[,1]=OLS
        # pval[,1]=summary(lm(Yt~X[,2]))$coefficients[,4]
        # rownames(beta)=c("beta0","beta2")
        # rownames(pval)=c("beta0","beta2")
        
      } else{
        beta=matrix(0,2,1)
        pval=matrix(0,2,1)
        conf=matrix(NA,2,2)
        
        
        
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2]-1,tau=tau)
          beta[2,1]=QR$coefficients
          beta[1,1]=NA
        }else{
          QR=rq(Yt~X-1,tau=tau)
          beta[,1]=QR$coefficients
        }
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[2,1]=summary[,4]
          conf[2,1]=(summary[,2]-(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
          conf[2,2]=(summary[,2]+(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
        }else{
          summary=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients
          pval[,1]=summary[,4]
          conf[,1]=(summary[,2]-(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
          conf[,2]=(summary[,2]+(qt(0.975,dim(Yt)[1]-2)*summary[,2]))
        }}
      
    }else{
      n=dim(Y)[1] #  # of firms
      t=dim(Y)[2] # time
      p=dim(Z)[2]
      s=length(subset)
      a=A-A*diag(dim(A)[2])
      w=a
      for(i in 1:n){
        if(sum(w[i,]!=0)){
          w[i,]=w[i,]/sum(w[i,])
        }
      }
      X=matrix(0,s*(t-1),2)
      Yt_1=Y[,1:(t-1)]
      w.Yt_1=w%*%Yt_1
      X[,1]=matrix(w.Yt_1[subset,],ncol=1)
      X[,2]=matrix(Yt_1[subset,],ncol = 1)
      z=Z[subset,]
      for(i in 1:(t-2)){z=rbind(z,Z[subset,])}
      X=cbind(X,z)
      Yt =  matrix(Y[,2:t][subset,],ncol = 1)
      if(intercept==T){
        beta=matrix(0,3+p,1)
        pval=matrix(0,3+p,1)
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2],tau=tau)
          beta[c(1,3:(3+p)),1]=QR$coefficients
          beta[2,1]=NA
        }else{
          QR=rq(Yt~X,tau=tau)
          beta[,1]=QR$coefficients
        }
        
        if(is.na(beta[2,1])){
          pval[2,1]=NA
          pval[c(1,3:(3+p)),1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }else{
          pval[,1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }
        
        
        rownames(beta)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta0","beta1","beta2",paste("gamma",c(1:p),sep=""))
      } else{
        beta=matrix(0,2+p,1)
        pval=matrix(0,2+p,1)
        if(length(unique(X[,1]))==1){
          QR=rq(Yt~X[,2]-1,tau=tau)
          beta[2:(2+p),1]=QR$coefficients
          beta[1,1]=NA
        }else{
          QR=rq(Yt~X-1,tau=tau)
          beta[,1]=QR$coefficients
        }
        if(is.na(beta[1,1])){
          pval[1,1]=NA
          pval[2:(2+p),1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }else{
          pval[,1]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
        }
        
        rownames(beta)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
        rownames(pval)=c("beta1","beta2",paste("gamma",c(1:p),sep=""))
      }
    }
  }
  beta=round(beta,4)
  pval=round(pval,4)
  returnlist <- list("Coefficients"= beta,
                     "Pvals" = pval,
                     "Conf" = conf
  )
  return(returnlist)
  
}
##############################################################################
QR_INVA=function(Y,A,intercept=T,tau){
  t=dim(Y)[1]
  n=dim(Y)[2]
  a=A-A*diag(n)
  w=a
  for(i in 1:n){
    if(sum(w[i,]!=0)){
      w[i,]=w[i,]/sum(w[i,])
    }
  }
  yt_1=Y[1:t-1,]
  yt=Y[2:t,]
  w.yt_1=yt_1%*%t(w)
  
  if(intercept==T){
    
    beta=matrix(0,3,n)
    pval=matrix(0,3,n)
    for(i in 1:n){
      X=matrix(c(w.yt_1[,i],yt_1[,i]),ncol = 2)
      if(length(unique(w.yt_1[,i]))==1){
        QR=rq(yt[,i]~X[,2],tau=tau)
        beta[c(1,3),i]=QR$coefficients
        beta[2,i]=NA
      }else{
        QR=rq(yt[,i]~X,tau=tau)
        beta[,i]=QR$coefficients
      }
      if(is.na(beta[2,i])){
        pval[2,i]=NA
        pval[c(1,3),i]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
      }else{
        pval[,i]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
      }}
    rownames(beta)=c("beta0","beta1","beta2")
    rownames(pval)=c("beta0","beta1","beta2")
  }else{
    beta=matrix(0,2,n)
    pval=matrix(0,2,n)
    for(i in 1:n){
      X=matrix(c(w.yt_1[,i],yt_1[,i]),ncol = 2)
      if(length(unique(w.yt_1[,i]))==1){
        QR=rq(yt[,i]~X[,2]-1,tau=tau)
        beta[2,i]=QR$coefficients
        beta[1,i]=NA
      }else{
        QR=rq(yt[,i]~X-1,tau=tau)
        beta[,i]=QR$coefficients
      }
      if(is.na(beta[1,i])){
        pval[1,i]=NA
        pval[2,i]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
      }else{
        pval[,i]=summary.rq(QR,se = "boot", bsmethod= "wild",R=1000)$coefficients[,4]
      }}
    
    rownames(beta)=c("beta1","beta2")
    rownames(pval)=c("beta1","beta2")
  }
  beta=round(beta,4)
  pval=round(pval,4)
  returnlist <- list("Coefficients"= beta,
                     "Pvals" = pval
  )
  return(returnlist)
}
##################################################################################
FDR=function(A,alpha){
  p=dim(A)[2]
  Adjacency=matrix(0,p,p)
  for(j in 1:p){
    yj=A[j,]
    Adjacency[j,order(yj)]=(sort(yj)<=((alpha-alpha/p)/(p-1)*(1:p)))
  }
  returnlist=list("Adjacency"=Adjacency)
  return(returnlist)
}
##################################################################################
# CV_LassoVAR=function(X,y,kfold,range){
#   source("lasso_inference.r")
#   p=ncol(X)
#   n=nrow(X)
#   s=sample(nrow(X))
#   X<-X[s,]
#   y=matrix(y,ncol=1)[s,]
#   y=as.matrix(y,ncol=1)
#   folds <- cut(seq(1,nrow(X)),breaks=kfold,labels=FALSE)
#   result=matrix(c(range,rep(NA,length(range))),ncol = 2);colnames(result)=c("Lambdas","CV")
#   
#   for(j in 1:length(range)){
#     cv.error=NULL
#     
#     for(i in 1:kfold){
#       
#       testIndexes <- which(folds==i,arr.ind=TRUE)
#       X_testData <- X[testIndexes, ]
#       y_testData <- y[testIndexes, ]
#       X_trainData <- X[-testIndexes, ]
#       y_trainData <- y[-testIndexes, ]
#       
#       beta = SSLasso(X_trainData,y_trainData,verbose = T,alpha = 0.05, lambda = range[j])$unb
#       pred=X_testData%*%beta
#       MSE=mean((y_testData-as.numeric(pred))^2)
#       cv.error=c(cv.error,MSE)
#     }
#     
#     cv=mean(cv.error)
#     result[j,2]=cv
#   }
#   return(result)
# }
# 
# 
# 
# 
# 
# 
# 
Matrix.G=function(A,beta1=NULL,beta2=NULL){
  n=dim(A)[1]
  a=A-A*diag(dim(A)[2])
  w=a
  for(i in 1:n){
    if(sum(w[i,]!=0)){
      w[i,]=w[i,]/sum(w[i,])
    }
  }
  if(is.null(beta1)&is.null(beta2)){
    
    return(list("w"=w))
  }
  
  G=beta1*w + beta2*diag(1,n)
  return(list("w"=w,"G"=G))
}

###################################################################################################
Influential.power=function(w,G,beta1=NULL,beta2=NULL){
  n=dim(G)[1]
  nu=solve((diag(1,n))-(t(G)))%*%matrix(1,n,1)
  
  if((is.null(beta1)&is.null(beta2))==F){
  approx.nu=((1-beta2)^-1)*matrix(1,n,1)+((1-beta2)^-2)*beta1*(t(w)%*%matrix(1,n,1))
  returnlist=list("nu"=nu,"approx.nu"=approx.nu)
  return(returnlist)
  }
  
  returnlist=list("nu"=nu)
  return(returnlist)
}
#########################################################################################
QR_lasso=function (x, y, tau = 0.5, lambda = NULL, intercept = TRUE, 
                   coef.cutoff = 1e-08, method = "br") 
{
  library(rqPen)
  if (is.null(dim(x))) {
    stop("x needs to be a matrix with more than 1 column")
  }
  
  p <- dim(x)[2]
  n <- dim(x)[1]
  xx=apply(x,2,function(y)scale(y,scale = F))
  sigma=diag(apply(xx,2,function(y)sqrt(sum(y^2)/n)))
  
  if (n != length(y)) {
    stop("length of y and rows of x do not match")
  }
  if (is.null(lambda) == TRUE | (length(lambda) != 1 & length(lambda) != 
                                 dim(x)[2])) {
    stop(paste("input of lambda must be of length 1 or", 
               dim(x)[2]))
  }
  if (sum(lambda < 0) > 0) {
    stop(paste("lambda must be positive and we have a lambda of ", 
               lambda, sep = ""))
  }
  
  
  
  lambda=lambda*sqrt(tau*(1-tau))/n
  lambda <- lambda * n
  
  if (length(lambda) == 1) {
    pen_x <- rbind(lambda*sigma, -lambda*sigma)
  }
  else {
    pen_x <- rbind(lambda*sigma, -lambda*sigma)
    pen_x <- pen_x[rowSums(pen_x == 0) != dim(pen_x)[2], 
                   ]
  }
  aug_n <- dim(pen_x)[1]
  aug_x <- rbind(x, pen_x)
  if (intercept) {
    aug_x <- cbind(c(rep(1, n), rep(0, aug_n)), aug_x)
  }
  aug_y <- c(y, rep(0, aug_n))
  
    model <- rq(aug_y ~ aug_x + 0, tau = tau)
  
  
  p_star <- p + intercept
  coefs <- coefficients(model)[1:p_star]
  return_val <- NULL
  return_val$coefficients <- coefs
  if (is.null(colnames(x))) {
    x_names <- paste("x", 1:p, sep = "")
  }
  else {
    x_names <- colnames(x)
  }
  if (intercept) {
    x_names <- c("intercept", x_names)
  }
  attributes(return_val$coefficients)$names <- x_names
  return_val$coefficients[abs(return_val$coefficients) < coef.cutoff] <- 0
  return_val$PenRho <- model$rho
  return_val$residuals <- model$residuals[1:n]
  
    return_val$rho <- sum(sapply(return_val$residuals, check, 
                                 tau))
  
  
  return_val$tau <- tau
  return_val$n <- n
  return_val$intercept <- intercept
  class(return_val) <- c("rq.pen", "rqLASSO")
  return_val
}
###############################################################################
LAMBDA=function(x,tau,B=500,c=1,alpha=0.9){
  x2=apply(x,2,function(y)scale(y,scale = F))
  n=nrow(x)
  p=ncol(x)
  tau=tau
  B=B
  lambda=0
  while(round(lambda,5)==0){
  Lambda=numeric(B)
  for(i in 1:B){
    k=numeric(p)
    set.seed(Sys.time())
    U=runif(n)
    for(j in 1:p){
      k[j]=abs(sum((x2[,j]*(tau-(U<=tau)*1))/(sqrt(sum(x2[,j]^2)/n)*sqrt(tau*(1-tau)))))/n
    }
    Lambda[i]=max(k)*n
  }
  lambda=c*as.numeric(quantile(Lambda,1-alpha))
  }
  
  return(lambda)
}

######################################################
Backtest=function(X,y,lambda,tau,from=0.01,to=2,length=100){
  n=nrow(X)
  library(caTools)
  interval=seq(from,to,length.out = length)
  C=NULL
  for(k in interval){
    lasso.QR=QR_lasso(x = X,y = y,tau = tau,lambda = k*lambda,intercept = F)$coefficients
    C=c(C,((sum((y<c(X%*%lasso.QR)))/n-tau)^2))
  }
  
  which.min(runmean(C,5))
  c=interval[which.min(runmean(C,5))]
  return(c)
}
###################################################################
Table=function(pre,true){
  table=table(pre!=0,true!=0)
  table[,"TRUE"]=as.numeric(table[,"TRUE"])
  table[,"FALSE"]=as.numeric(table[,"FALSE"])
  Precision=table["TRUE","TRUE"]/(table["TRUE","FALSE"]+table["TRUE","TRUE"])
  Recall=table["TRUE","TRUE"]/(table["FALSE","TRUE"]+table["TRUE","TRUE"])
  # MCC=(table["TRUE","TRUE"]*table["FALSE","FALSE"]-table["TRUE","FALSE"]*table["FALSE","TRUE"])/
  # sqrt(sum(table["TRUE",])*sum(table[,"TRUE"])*sum(table[,"FALSE"])*sum(table["FALSE",]))
  # ERR=sum((pre!=0)[,1:(dim(pre)[2]-p) ])/sum(pre!=0)
  A=c(Precision,Recall)
  names(A)=c("Precision","Recall")
  return(A)
}
####################################################################
MSE=function(A,X.test,y.test){
  iMSE=numeric(nrow(X.test))
  MSE=numeric(nrow(X.test))
  for(i in 1:nrow(X.test)){
    iMSE[i] = mean((y.test[i,]-(A %*% X.test[i,]))^2)
  }
  for(j in 1:nrow(X.test)){
    MSE[j]=mean(iMSE[1:j])
    
  }
  return(MSE)
}
######################################################################
Adalasso=function(X.train,y.train,method,lambda=seq(0.1*sqrt(2*log(p)/n),10*sqrt(2*log(p)/n),length.out = 100),intercept=intercept){
  n=nrow(X.train)
  p=ncol(X.train)
  
  if(method=="lasso"){
  lasso <- glmnet(X.train, y.train,alpha=1,lambda = lambda,intercept = intercept)
  
  hqc=numeric(length(lambda))
  
  for(i in 1:length(lambda)){
    betas=lasso$beta[,i]
    k=sum(betas!=0)
    hqc[i]=n*log(mean((y.train-(X.train%*%betas))^2))+2*k*log(log(n))
  }
  HQC.lasso.index=which.min(hqc)
  HQC.lasso=min(hqc)
  
  postlasso=glmnet(X.train, y.train,alpha=1,lambda = lambda[HQC.lasso.index],intercept = intercept)$beta
  
  w3 <- 1/abs(matrix(postlasso))
  
  w3[w3[,1] == Inf] <- 999999999 ## Replacing values estimated as Infinite for 999999999
  
  }else if(method=="LSE"){
    
    if(intercept==T){
    LSE=lm(y.train~X.train)$coefficients
    }else{LSE=lm(y.train~X.train-1)$coefficients}
    w3 <- 1/abs(matrix(LSE))
  }
  ## Adaptive Lasso
  
  ad.lasso <- glmnet(X.train, y.train,  alpha=1, standardize=TRUE,  penalty.factor=w3,lambda = lambda,intercept = intercept)
  
  hqc2=numeric(length(lambda))
  
  for(i in 1:length(lambda)){
    betas2=ad.lasso$beta[,i]
    k2=sum(betas2!=0)
    hqc2[i]=n*log(mean((y.train-(X.train%*%betas2))^2))+2*k2*log(log(n))
  }
  HQC.adlasso.index=which.min(hqc2)
  HQC.adlasso=min(hqc2)
  adaLasso=glmnet(X.train, y.train,  alpha=1, standardize=TRUE,  penalty.factor=w3,lambda = lambda[HQC.adlasso.index],intercept = intercept)$beta
  AdaLasso=as.numeric(adaLasso)
  return(list("AdaLasso"=AdaLasso,"adlasso.index"=HQC.adlasso.index))
}
