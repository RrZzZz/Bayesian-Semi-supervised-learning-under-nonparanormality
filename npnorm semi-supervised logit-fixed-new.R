# install.packages('statmod')
# install.packages('quadprog')
# install.packages('matrixcalc')
# install.packages('mvtnorm')
# install.packages('splines')
# install.packages('tmg')
# install.packages('tmvmixnorm')

setwd('//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/research3/sim-difftrans')

library(mvtnorm)
library(splines)
library(tmg)
library(tmvmixnorm)
library(statmod)
library(quadprog)
library(matrixcalc)
library(Matrix)

set.seed(999)
seed=sample(1:1000, 30)
fp=fn=array(dim=c(3,2,3,30))
pvalues=c(5,10,15)
n.total=c(50, 100)
n.label=c(3, 5, 10)


for(p.num in 1:3){
  set.seed(999)
  p=pvalues[p.num]
  cov.gen <- function(p){
    A1 <- matrix(runif(p^2)*2-1, ncol=p) 
    return(t(A1) %*% A1)
  }
  
  mean1.true=runif(p)*4
  mean2.true=runif(p)*4
  sigma1.true=cov.gen(p)
  sigma2.true=cov.gen(p)
  
  for(total.num in 1:2){
    for(n.num in 1:3){
      for(sim.num in 1:30){
        sim=seed[sim.num]
        set.seed(sim)
        niter=500
        #################generate data##########################
        n1=n2=n.total[total.num]
        n=n1+n2
        
        x.1=rmvnorm(n1, mean = mean1.true, sigma = sigma1.true)
        x.2=rmvnorm(n2, mean = mean2.true, sigma = sigma2.true)
        for(i in 1:p){
          x.1[,i]=1/(1+exp(-(x.1[,i]-0.5*(mean1.true[i]+mean2.true[i]))/((0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i]))))))
          x.2[,i]=1/(1+exp(-(x.2[,i]-0.5*(mean1.true[i]+mean2.true[i]))/((0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i]))))))
        }
        
        x=rbind(x.1,x.2)
  
        L.true=c(rep(0,n1),rep(1,n2))
        L.org=rep(2, n)
        L.org[sample(n1, size=n.label[n.num])]=0
        L.org[n1+sample(n2, size=n.label[n.num])]=1
        
        
        midpoint=rep(0, 15)
        
        for(J in 8:15){
          knotvec=c(rep(0,3), seq(0, 1, by=1/J), rep(1,3))
          A.bas=splineDesign(knots=knotvec, x = c(1/4, 1/2, 3/4))
          A=rbind(A.bas[2,], A.bas[3,]-A.bas[1,])
          
          J1=dim(A)[2]
          mu.org=qnorm((1:J1-0.375)/(J1+0.25))
          mu.prior=mu.org+t(A)%*%solve(A%*%t(A))%*%(matrix(c(0,1),2,1)-A%*%mu.org)
          Sigma.prior=5*(diag(rep(1, J1))-t(A)%*%solve(A%*%t(A))%*%A)
          
          # solve A*theta=(0 1)^T
          index1=which.max(A[1,])
          index2=which.max(A[2,])
          if(index1==index2){index2=index2+1}
          
          a.index1=A[1,index1]
          a.index2=A[1,index2]
          b.index1=A[2,index1]
          b.index2=A[2,index2]
          
          W=W.org=A[, -c(index1, index2)]
          W[2,]=(W.org[2,]-W.org[1,]*b.index1/a.index1)/-(b.index2-a.index2*b.index1/a.index1)
          W[1,]=-W.org[1,]/a.index1-W[2,]*a.index2/a.index1
          q=c(1/(b.index2-a.index2*b.index1/a.index1)*(-a.index2/a.index1), 1/(b.index2-a.index2*b.index1/a.index1))
          
          # check the calculation
          A[1,c(index1, index2)]%*%W+c(A[1,-c(index1, index2)])
          A[1,c(index1, index2)]%*%q
          A[2,c(index1, index2)]%*%W+c(A[2,-c(index1, index2)])
          A[2,c(index1, index2)]%*%q
          
          mu.prior.r=mu.prior[-c(index1, index2)]
          Sigma.prior.r=Sigma.prior[-c(index1, index2), -c(index1, index2)]
          Sigma.prior.r.inv=chol2inv(chol(Sigma.prior.r))
          
          newindex=function(index){
            count=(index<index1)+(index<index2)
            if(count==2){return(index)}
            if(count==1){return(index-1)}
            if(count==0){return(index-2)}
          }
          
          F.org=cbind(diag(rep(-1, J1-1)), rep(0, J1-1))+cbind(rep(0, J1-1), diag(rep(1, J1-1)))
          F=matrix(, 1,J1-2)
          g=rep(0, J1)
          for(i in 2:J1){
            new=rep(0, J1-2)
            if(i==index1){
              new=W[1,]
              g[i]=q[1]
            }else if(i==index2){
              new=W[2,]
              g[i]=q[2]
            }else{
              new[newindex(i)]=1
            }
            if(i-1==index1){
              new=new-W[1,]
              g[i]=g[i]-q[1]
            }else if(i-1==index2){
              new=new-W[2,]
              g[i]=g[i]-q[2]
            }else{
              new[newindex(i-1)]=new[newindex(i-1)]-1
            }
            F=rbind(F, new)
          }
          F=F[-1,]
          g=g[-1]
          
          theta.iter=array(dim=c(niter, J1-2, p))
          B=array(dim=c(n, p, J1-2)) 
          B.star=array(dim=c(n, p, 2))
          for(d in 1:p){
            middle=splineDesign(knots=knotvec, x = x[,d])
            B[,d,] = middle[,-c(index1, index2)]
            B.star[,d,] = middle[,c(index1, index2)]
          }
          
          # initial for theta
          out=gauss.quad(20,kind="hermite",alpha=0,beta=0)
          bspline.matrix.cdf=splineDesign(knots = knotvec, x=pnorm(out$nodes))
          fn.LHS=bspline.matrix.cdf*out$nodes
          fn.LHS.weighted=out$weights*fn.LHS
          LHS.gq=colSums(fn.LHS.weighted)
          
          fn.RHS.weighted=matrix(0, J1, J1)
          RHS.gq=matrix(0, J1, J1)
          for(pts in 1:20){
            for(j in 1:J1){
              for(k in 1:J1){
                fn.RHS.weighted[j,k]=out$weights[pts]*bspline.matrix.cdf[pts, j]*bspline.matrix.cdf[pts, k]
              }
            }
            RHS.gq=RHS.gq+fn.RHS.weighted
          }
          theta.initial=solve.QP(nearPD(RHS.gq%*%RHS.gq)$'mat', as.numeric(LHS.gq%*%RHS.gq), Amat=t(rbind(A, F.org)), bvec=c(0, 1, rep(10^(-4), J1-1)), meq=2)$'solution'
          theta=matrix(rep(theta.initial[-c(index1, index2)], p), ,p)
          
          Y=matrix(,n,p)
          for(d in 1:p){
            Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
          }
          
          # initial values
          if(sum(L.org==0)==1){
            mu1=Y[L.org==0,]
          }else{
            mu1=colMeans(Y[L.org==0,]) 
          }
          
          if(sum(L.org==1)==1){
            mu2=Y[L.org==1,]
          }else{
            mu2=colMeans(Y[L.org==1,]) 
          }
          
          
          # Assign initial labels according to distance
          L=L.org
          for(i in 1:n){
            if(L.org[i]==2){
              d1=sum((Y[i,]-mu1)^2)
              d2=sum((Y[i,]-mu2)^2)
              if(d1<d2){
                L[i]=0
              }else{
                L[i]=1
              }
            }
          }
          n1=sum(L==0)
          n2=sum(L==1)
          print(sum(L!=L.true))
          
          ################check the original L#################
          
          Sigma1=cov(Y[L==0,])
          Sigma2=cov(Y[L==1,])
          Sigma1inv=chol2inv(chol(Sigma1))
          is.positive.definite(Sigma1inv)
          Sigma2inv=chol2inv(chol(Sigma2))
          is.positive.definite(Sigma2inv)
          
          mu1.iter=matrix(,niter,p)
          Sigma1.iter=array(dim=c(niter, p,p))
          Sigma1inv.iter=array(dim=c(niter, p,p))
          mu2.iter=matrix(,niter,p)
          Sigma2.iter=array(dim=c(niter, p,p))
          Sigma2inv.iter=array(dim=c(niter, p,p))
          L.iter=matrix(,niter,n)
          
          mean.prior=mu.prior.r%*%Sigma.prior.r.inv
          mid=array(dim=c(n, p, J1-2)) 
          midsq=array(dim=c(n, p, J1-2, J1-2))
          mid2=matrix(,n, p)
          
          for(d in 1:p){
            mid[,d,]=B[,d,]+B.star[,d,]%*%W
            mid2[,d]=B.star[,d,]%*%matrix(q)
            for(i in 1:n){
              midsq[i,d,,]=(mid[i,d,])%*%t(mid[i,d,])
            }
          }
          
          
          for(iter in 1:niter){
            for(d in 1:p){
              Cov=matrix(0, J1-2, J1-2)
              mean1=0
              mean2=0
              var1=matrix(0,J1-2,J1-2)
              var2=matrix(0,J1-2,J1-2)
              for(i in 1:n){
                if(L[i]==0){
                  mean1=mean1+mid[i,d,]*((mid2[i,d]-mu1[d])*Sigma1inv[d, d]+sum(Sigma1inv[-d, d]*(Y[i,-d]-mu1[-d])))
                  var1=var1+midsq[i,d,,]
                }
                if(L[i]==1){
                  mean2=mean2+mid[i,d,]*((mid2[i,d]-mu2[d])*Sigma2inv[d, d]+sum(Sigma2inv[-d, d]*(Y[i,-d]-mu2[-d])))
                  var2=var2+midsq[i,d,,]
                }
              }
              Cov=Sigma1inv[d, d]*var1+Sigma2inv[d, d]*var2+Sigma.prior.r.inv
              mean.final=mean.prior-mean1-mean2
              Cov=chol2inv(chol(Cov))
              if(!is.positive.definite(Cov)){
                Cov=nearPD(Cov)$'mat'
              }
              theta[,d]=rtmvn(n=1, Mean=as.numeric(Cov%*%matrix(mean.final)), Sigma = Cov, D=F, lower=-g, upper=rep(Inf, length(g)), int = theta.initial[-c(index1, index2)])
              Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
              
            }
            L.iter[iter,]=L
            theta.iter[iter,,]=theta
            Sigma1inv.iter[iter,,]=Sigma1inv=rWishart(1, n1, chol2inv(chol(cov(Y[L==0,])*(n1-1))))[,,1]
            Sigma1.iter[iter,,]=  Sigma1=chol2inv(chol(Sigma1inv))
            mu1.iter[iter,]=  mu1=rmvnorm(1, colMeans(Y[L==0,]), Sigma1/n1)
            Sigma2inv.iter[iter,,]=Sigma2inv=rWishart(1, n2, chol2inv(chol(cov(Y[L==1,])*(n2-1))))[,,1]
            Sigma2.iter[iter,,]=  Sigma2=chol2inv(chol(Sigma2inv))
            mu2.iter[iter,]=  mu2=rmvnorm(1, colMeans(Y[L==1,]), Sigma2/n2)
            
            for(i in 1:n){
              if(L.org[i]==2){
                if(dmvnorm(Y[i,], mu1, Sigma1)>dmvnorm(Y[i,], mu2, Sigma2)){
                  L[i]=0
                }else{
                  L[i]=1
                }
              }
            }
            n1=sum(L==0)
            n2=sum(L==1)
            if(iter%%100==0){
              print(iter)
              print(sum(L!=L.true))
            }
          }
          theta.final=apply(theta.iter[(niter/10):niter,,], 2:3, mean)
          mu1.final=apply(mu1.iter[(niter/10):niter,], 2, mean)
          Sigma1.final=apply(Sigma1.iter[(niter/10):niter,,], 2:3, mean)
          mu2.final=apply(mu2.iter[(niter/10):niter,], 2, mean)
          Sigma2.final=apply(Sigma2.iter[(niter/10):niter,,], 2:3, mean)
          
          for(d in 1:p){
            Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta.final[,d]+B.star[,d,]%*%matrix(q)
          }
          
          p1=dmvnorm(Y, mu1.final, Sigma1.final)
          p2=dmvnorm(Y, mu2.final, Sigma2.final)
          r=p1/p2
          midpoint[J]=sum(r[which(r>1/3)]<3)
        }
        
        ############################################
        niter=10000
        J=7+which.min(midpoint[8:15])
        knotvec=c(rep(0,3), seq(0, 1, by=1/J), rep(1,3))
        A.bas=splineDesign(knots=knotvec, x = c(1/4, 1/2, 3/4))
        A=rbind(A.bas[2,], A.bas[3,]-A.bas[1,])
        
        J1=dim(A)[2]
        mu.org=qnorm((1:J1-0.375)/(J1+0.25))
        mu.prior=mu.org+t(A)%*%solve(A%*%t(A))%*%(matrix(c(0,1),2,1)-A%*%mu.org)
        Sigma.prior=5*(diag(rep(1, J1))-t(A)%*%solve(A%*%t(A))%*%A)
        
        # solve A*theta=(0 1)^T
        index1=which.max(A[1,])
        index2=which.max(A[2,])
        if(index1==index2){index2=index2+1}
        
        a.index1=A[1,index1]
        a.index2=A[1,index2]
        b.index1=A[2,index1]
        b.index2=A[2,index2]
        
        W=W.org=A[, -c(index1, index2)]
        W[2,]=(W.org[2,]-W.org[1,]*b.index1/a.index1)/-(b.index2-a.index2*b.index1/a.index1)
        W[1,]=-W.org[1,]/a.index1-W[2,]*a.index2/a.index1
        q=c(1/(b.index2-a.index2*b.index1/a.index1)*(-a.index2/a.index1), 1/(b.index2-a.index2*b.index1/a.index1))
        
        # check the calculation
        A[1,c(index1, index2)]%*%W+c(A[1,-c(index1, index2)])
        A[1,c(index1, index2)]%*%q
        A[2,c(index1, index2)]%*%W+c(A[2,-c(index1, index2)])
        A[2,c(index1, index2)]%*%q
        
        mu.prior.r=mu.prior[-c(index1, index2)]
        Sigma.prior.r=Sigma.prior[-c(index1, index2), -c(index1, index2)]
        Sigma.prior.r.inv=chol2inv(chol(Sigma.prior.r))
        
        newindex=function(index){
          count=(index<index1)+(index<index2)
          if(count==2){return(index)}
          if(count==1){return(index-1)}
          if(count==0){return(index-2)}
        }
        
        F.org=cbind(diag(rep(-1, J1-1)), rep(0, J1-1))+cbind(rep(0, J1-1), diag(rep(1, J1-1)))
        F=matrix(, 1,J1-2)
        g=rep(0, J1)
        for(i in 2:J1){
          new=rep(0, J1-2)
          if(i==index1){
            new=W[1,]
            g[i]=q[1]
          }else if(i==index2){
            new=W[2,]
            g[i]=q[2]
          }else{
            new[newindex(i)]=1
          }
          if(i-1==index1){
            new=new-W[1,]
            g[i]=g[i]-q[1]
          }else if(i-1==index2){
            new=new-W[2,]
            g[i]=g[i]-q[2]
          }else{
            new[newindex(i-1)]=new[newindex(i-1)]-1
          }
          F=rbind(F, new)
        }
        F=F[-1,]
        g=g[-1]
        
        theta.iter=array(dim=c(niter, J1-2, p))
        B=array(dim=c(n, p, J1-2)) 
        B.star=array(dim=c(n, p, 2))
        for(d in 1:p){
          middle=splineDesign(knots=knotvec, x = x[,d])
          B[,d,] = middle[,-c(index1, index2)]
          B.star[,d,] = middle[,c(index1, index2)]
        }
        
        # initial for theta
        out=gauss.quad(20,kind="hermite",alpha=0,beta=0)
        bspline.matrix.cdf=splineDesign(knots = knotvec, x=pnorm(out$nodes))
        fn.LHS=bspline.matrix.cdf*out$nodes
        fn.LHS.weighted=out$weights*fn.LHS
        LHS.gq=colSums(fn.LHS.weighted)
        
        fn.RHS.weighted=matrix(0, J1, J1)
        RHS.gq=matrix(0, J1, J1)
        for(pts in 1:20){
          for(j in 1:J1){
            for(k in 1:J1){
              fn.RHS.weighted[j,k]=out$weights[pts]*bspline.matrix.cdf[pts, j]*bspline.matrix.cdf[pts, k]
            }
          }
          RHS.gq=RHS.gq+fn.RHS.weighted
        }
        theta.initial=solve.QP(nearPD(RHS.gq%*%RHS.gq)$'mat', as.numeric(LHS.gq%*%RHS.gq), Amat=t(rbind(A, F.org)), bvec=c(0, 1, rep(10^(-4), J1-1)), meq=2)$'solution'
        theta=matrix(rep(theta.initial[-c(index1, index2)], p), ,p)
        
        Y=matrix(,n,p)
        for(d in 1:p){
          Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
        }
        
        # initial values
        if(sum(L.org==0)==1){
          mu1=Y[L.org==0,]
        }else{
          mu1=colMeans(Y[L.org==0,]) 
        }
        
        if(sum(L.org==1)==1){
          mu2=Y[L.org==1,]
        }else{
          mu2=colMeans(Y[L.org==1,]) 
        }
        
        
        # Assign initial labels according to distance
        L=L.org
        for(i in 1:n){
          if(L.org[i]==2){
            d1=sum((Y[i,]-mu1)^2)
            d2=sum((Y[i,]-mu2)^2)
            if(d1<d2){
              L[i]=0
            }else{
              L[i]=1
            }
          }
        }
        n1=sum(L==0)
        n2=sum(L==1)
        sum(L!=L.true)
        
        ################check the original L#################
        
        Sigma1=cov(Y[L==0,])
        Sigma2=cov(Y[L==1,])
        Sigma1inv=chol2inv(chol(Sigma1))
        is.positive.definite(Sigma1inv)
        Sigma2inv=chol2inv(chol(Sigma2))
        is.positive.definite(Sigma2inv)
        
        mu1.iter=matrix(,niter,p)
        Sigma1.iter=array(dim=c(niter, p,p))
        Sigma1inv.iter=array(dim=c(niter, p,p))
        mu2.iter=matrix(,niter,p)
        Sigma2.iter=array(dim=c(niter, p,p))
        Sigma2inv.iter=array(dim=c(niter, p,p))
        L.iter=matrix(,niter,n)
        
        mean.prior=mu.prior.r%*%Sigma.prior.r.inv
        mid=array(dim=c(n, p, J1-2)) 
        midsq=array(dim=c(n, p, J1-2, J1-2))
        mid2=matrix(,n, p)
        
        for(d in 1:p){
          mid[,d,]=B[,d,]+B.star[,d,]%*%W
          mid2[,d]=B.star[,d,]%*%matrix(q)
          for(i in 1:n){
            midsq[i,d,,]=(mid[i,d,])%*%t(mid[i,d,])
          }
        }
        
        
        for(iter in 1:niter){
          for(d in 1:p){
            Cov=matrix(0, J1-2, J1-2)
            mean1=0
            mean2=0
            var1=matrix(0,J1-2,J1-2)
            var2=matrix(0,J1-2,J1-2)
            for(i in 1:n){
              if(L[i]==0){
                mean1=mean1+mid[i,d,]*((mid2[i,d]-mu1[d])*Sigma1inv[d, d]+sum(Sigma1inv[-d, d]*(Y[i,-d]-mu1[-d])))
                var1=var1+midsq[i,d,,]
              }
              if(L[i]==1){
                mean2=mean2+mid[i,d,]*((mid2[i,d]-mu2[d])*Sigma2inv[d, d]+sum(Sigma2inv[-d, d]*(Y[i,-d]-mu2[-d])))
                var2=var2+midsq[i,d,,]
              }
            }
            Cov=Sigma1inv[d, d]*var1+Sigma2inv[d, d]*var2+Sigma.prior.r.inv
            mean.final=mean.prior-mean1-mean2
            Cov=chol2inv(chol(Cov))
            if(!is.positive.definite(Cov)){
              Cov=nearPD(Cov)$'mat'
            }
            theta[,d]=rtmvn(n=1, Mean=as.numeric(Cov%*%matrix(mean.final)), Sigma = Cov, D=F, lower=-g, upper=rep(Inf, length(g)), int = theta.initial[-c(index1, index2)])
            Y[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta[,d]+B.star[,d,]%*%matrix(q)
            
          }
          L.iter[iter,]=L
          theta.iter[iter,,]=theta
          Sigma1inv.iter[iter,,]=Sigma1inv=rWishart(1, n1, chol2inv(chol(cov(Y[L==0,])*(n1-1))))[,,1]
          Sigma1.iter[iter,,]=  Sigma1=chol2inv(chol(Sigma1inv))
          mu1.iter[iter,]=  mu1=rmvnorm(1, colMeans(Y[L==0,]), Sigma1/n1)
          Sigma2inv.iter[iter,,]=Sigma2inv=rWishart(1, n2, chol2inv(chol(cov(Y[L==1,])*(n2-1))))[,,1]
          Sigma2.iter[iter,,]=  Sigma2=chol2inv(chol(Sigma2inv))
          mu2.iter[iter,]=  mu2=rmvnorm(1, colMeans(Y[L==1,]), Sigma2/n2)
          
          for(i in 1:n){
            if(L.org[i]==2){
              if(dmvnorm(Y[i,], mu1, Sigma1)>dmvnorm(Y[i,], mu2, Sigma2)){
                L[i]=0
              }else{
                L[i]=1
              }
            }
          }
          n1=sum(L==0)
          n2=sum(L==1)
          if(iter%%100==0){
            print(iter)
            print(sum(L!=L.true))
          }
        }
        theta.final=apply(theta.iter[(niter/10):niter,,], 2:3, mean)
        mu1.final=apply(mu1.iter[(niter/10):niter,], 2, mean)
        Sigma1.final=apply(Sigma1.iter[(niter/10):niter,,], 2:3, mean)
        mu2.final=apply(mu2.iter[(niter/10):niter,], 2, mean)
        Sigma2.final=apply(Sigma2.iter[(niter/10):niter,,], 2:3, mean)
        
        
        
        
        
        ############testing###########
        n.train=5000
        train.1=rmvnorm(n.train, mean = mean1.true, sigma = sigma1.true)
        train.2=rmvnorm(n.train, mean = mean2.true, sigma = sigma2.true)
        
        train=rbind(train.1, train.2)
        for(i in 1:p){
          train[,i]=1/(1+exp(-(train[,i]-0.5*(mean1.true[i]+mean2.true[i]))/(0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i])))))
        }
        B=array(dim=c(nrow(train), p, J1-2)) 
        B.star=array(dim=c(nrow(train), p, 2))
        Y.train=matrix(,2*n.train,p)
        for(d in 1:p){
          middle=splineDesign(knots=knotvec, x = train[,d])
          B[,d,] = middle[,-c(index1, index2)]
          B.star[,d,] = middle[,c(index1, index2)]
          Y.train[,d]=(B[,d,]+B.star[,d,]%*%W)%*%theta.final[,d]+B.star[,d,]%*%matrix(q)
        }
        p1=dmvnorm(Y.train, mu1.final, Sigma1.final)
        p2=dmvnorm(Y.train, mu2.final, Sigma2.final)
        r=p1<p2
        fp[p.num, total.num, n.num, sim.num]=sum(r[1:n.train])/n.train
        fn[p.num, total.num, n.num, sim.num]=1-sum(r[n.train+1:n.train])/n.train
        
        save.image(file='result-mean4-diff-p1.RData')
      }

    }
  }
}
