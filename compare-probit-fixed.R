library(RSSL)
library(mvtnorm)

setwd('//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/research3/sim-new1')

set.seed(999)
seed=sample(1:1000, 30)
pvalues=c(5,10,15)
n.total=c(50, 100)
n.label=c(3, 5, 10)

class=list()
fp=fn=list()
for(i in 1:22){
  fp[[i]]=array(dim=c(3,2,3,30))
  fn[[i]]=array(dim=c(3,2,3,30))
}


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
        
        x=rbind(x.1,x.2)
        for(i in 1:p){
          x[,i]=pnorm(x[,i], 0.5*(mean1.true[i]+mean2.true[i]), 0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i])))
        }
        
        L.true=c(rep(0,n1),rep(1,n2))
        L.org=rep(NA, n)
        L.org[sample(n1, size=n.label[n.num])]=0
        L.org[n1+sample(n2, size=n.label[n.num])]=1
      
      L.org=as.factor(L.org)
      df=data.frame(x, L.org)
      
      n.train=5000
      train.1=rmvnorm(n.train, mean = mean1.true, sigma = sigma1.true)
      train.2=rmvnorm(n.train, mean = mean2.true, sigma = sigma2.true)
      
      train=rbind(train.1, train.2)
      # train=exp(train)/(1+exp(train))
      for(i in 1:p){
        train[,i]=pnorm(train[,i], 0.5*(mean1.true[i]+mean2.true[i]), 0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i])))
      }
      L.train=as.factor(c(rep(0, n.train), rep(1, n.train)))
      df.train=data.frame(train, L.train)
      
      class[[1]] <- EntropyRegularizedLogisticRegression(L.org~.,df, lambda=0.01,lambda_entropy = 100)
      class[[2]] <- EMLeastSquaresClassifier(L.org~.,df, max_iter = 5000)
      # class[[3]] <- EMLinearDiscriminantClassifier(L.org~., df, max_iter = 5000)
      class[[4]] <- ICLeastSquaresClassifier(L.org~., df)
      # class[[5]] <- ICLinearDiscriminantClassifier(L.org~., df)
      class[[6]] <- LaplacianSVM(L.org~., df)
      class[[7]] <- LeastSquaresClassifier(L.org~., df)
      # class[[8]] <- LinearDiscriminantClassifier(L.org~., df)
      class[[9]] <- LinearSVM(L.org~., df)
      class[[10]] <- LogisticLossClassifier(L.org~., df)
      class[[11]] <- LogisticRegression(L.org~., df)
      class[[12]] <- MCLinearDiscriminantClassifier(L.org~., df)
      class[[13]] <- MCNearestMeanClassifier(L.org~., df)
      # class[[14]] <- MCPLDA(L.org~., df)
      class[[15]] <- MajorityClassClassifier(L.org~., df)
      class[[16]] <- NearestMeanClassifier(L.org~., df)
      # class[[17]] <- QuadraticDiscriminantClassifier(L.org~., df)
      class[[18]] <- SVM(L.org~., df)
      class[[19]] <- USMLeastSquaresClassifier(L.org~., df)
      class[[20]] <- WellSVM(L.org~., df)
      class[[21]] <- svmlin(L.org~., df)
      class[[22]] <- EMNearestMeanClassifier(L.org~.,df)
      # class[[23]] <- KernelLeastSquaresClassifier(L.org~., df)
      # class[[24]] <- LaplacianKernelLeastSquaresClassifier(L.org~., df)
      # class[[25]] <- LinearTSVM(L.org~., df)
      # class[[26]] <- S4VM(L.org~., df)
      # class[[27]] <- SelfLearning(L.org~., df)
      # class[[28]] <- TSVM(L.org~., df)
      irange=c(1:2,4,6:7,9:11,13,15:16, 18:22)
      #       
      for(i in irange){
        fp[[i]][p.num, total.num, n.num, sim.num]=sum(predict(class[[i]],df.train)[1:n.train]==1)/n.train
        fn[[i]][p.num, total.num, n.num, sim.num]=sum(predict(class[[i]],df.train)[n.train+1:n.train]==0)/n.train
      }

     } 
    }
  }
}
#17,

save.image('compare-probit.RData')
# 
for(i in irange){
  print(i)
  print(apply(fp[[i]],1:3, mean))
  print(apply(fn[[i]],1:3, mean))

}
# fp[[5]][3,,]
# 
# for(i in 1:22){
#   print(i)
#   print(sum(predict(class[[i]],df.train)[1:n.train]==1)/n.train)
#   print(sum(predict(class[[i]],df.train)[n.train+1:n.train]==0)/n.train)
#   
# }
