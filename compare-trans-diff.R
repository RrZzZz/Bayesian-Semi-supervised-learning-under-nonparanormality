library(RSSL)
library(mvtnorm)

setwd('//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/research3/sim-difftrans')

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

for(p.num in 1){
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
          x.2[,i]=pnorm(x.2[,i], 0.5*(mean1.true[i]+mean2.true[i]), 0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i])) )
        }
        
        x=rbind(x.1,x.2)
        
        L.true=c(rep(0,n1),rep(1,n2))
        L.org=rep(NA, n)
        L.org[sample(n1, size=n.label[n.num])]=0
        L.org[n1+sample(n2, size=n.label[n.num])]=1
      
      L.org=as.factor(L.org)
      df=data.frame(x, L.org)
      
      n.train=5000
      train.1=rmvnorm(n.train, mean = mean1.true, sigma = sigma1.true)
      train.2=rmvnorm(n.train, mean = mean2.true, sigma = sigma2.true)
      for(i in 1:p){
        train.1[,i]=1/(1+exp(-(train.1[,i]-0.5*(mean1.true[i]+mean2.true[i]))/((0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i]))))))
        train.2[,i]=pnorm(train.2[,i], 0.5*(mean1.true[i]+mean2.true[i]), 0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i])) )
      }
      
      train=rbind(train.1, train.2)
      # for(i in 1:p){
      #   train[,i]=pnorm(train[,i], 0.5*(mean1.true[i]+mean2.true[i]), 0.5*(sqrt(sigma1.true[i,i])+sqrt(sigma2.true[i,i])))
      # }
      L.train=as.factor(c(rep(0, n.train), rep(1, n.train)))
      df.train=data.frame(train, L.train)
      
      # class[[1]] <- EntropyRegularizedLogisticRegression(L.org~.,df, lambda=0.01,lambda_entropy = 100)
      # class[[2]] <- EMLeastSquaresClassifier(L.org~.,df, max_iter = 5000)
      # class[[3]] <- EMLinearDiscriminantClassifier(L.org~., df, max_iter = 5000)
      class[[4]] <- ICLeastSquaresClassifier(L.org~., df)
      # class[[5]] <- ICLinearDiscriminantClassifier(L.org~., df)
      class[[6]] <- LaplacianSVM(L.org~., df)
      # class[[7]] <- LeastSquaresClassifier(L.org~., df)
      # class[[8]] <- LinearDiscriminantClassifier(L.org~., df)
      # class[[9]] <- LinearSVM(L.org~., df)
      # class[[10]] <- LogisticLossClassifier(L.org~., df)
      # class[[11]] <- LogisticRegression(L.org~., df)
      # class[[12]] <- MCLinearDiscriminantClassifier(L.org~., df)
      # class[[13]] <- MCNearestMeanClassifier(L.org~., df)
      # # class[[14]] <- MCPLDA(L.org~., df)
      # class[[15]] <- MajorityClassClassifier(L.org~., df)
      # class[[16]] <- NearestMeanClassifier(L.org~., df)
      # # class[[17]] <- QuadraticDiscriminantClassifier(L.org~., df)
      # class[[18]] <- SVM(L.org~., df)
      # class[[19]] <- USMLeastSquaresClassifier(L.org~., df)
      class[[20]] <- WellSVM(L.org~., df)
      class[[21]] <- svmlin(L.org~., df)
      class[[22]] <- EMNearestMeanClassifier(L.org~.,df)
      # class[[23]] <- KernelLeastSquaresClassifier(L.org~., df)
      # class[[24]] <- LaplacianKernelLeastSquaresClassifier(L.org~., df)
      # class[[25]] <- LinearTSVM(L.org~., df)
      # class[[26]] <- S4VM(L.org~., df)
      # class[[27]] <- SelfLearning(L.org~., df)
      # class[[28]] <- TSVM(L.org~., df)
      irange=c(4,6,20:22)
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

save.image('compare-logit-new.RData')
# 
for(i in c(4,6,20:22)){
  print(i)
  # print('fp')
  # print(apply(fp[[i]][1,,,] ,1:2, mean))
  # print('fn')
  # print(apply(fn[[i]][1,,,] ,1:2, mean))
  print('error')
  print((apply(fn[[i]][1,,,] ,1:2, mean)+apply(fp[[i]][1,,,] ,1:2, mean))/2*100)
  
}

i=6
print('error')
print((apply(fn[[i]][1,,,] ,1:2, function(x){mean(x, na.rm = T)})+apply(fp[[i]][1,,,] ,1:2, function(x){mean(x, na.rm = T)}))/2*100)
# fp[[5]][3,,]
# 
# for(i in 1:22){
#   print(i)
#   print(sum(predict(class[[i]],df.train)[1:n.train]==1)/n.train)
#   print(sum(predict(class[[i]],df.train)[n.train+1:n.train]==0)/n.train)
#   
# }
