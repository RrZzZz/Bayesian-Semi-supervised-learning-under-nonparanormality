library(RSSL)
library(mvtnorm)

iono=read.table('//wolftech.ad.ncsu.edu/cos/stat/Redirect/rzhu4/Downloads/ionosphere.data', sep=',')

summary(iono)

fp=list()
fn=list()
error=list()
mcc=list()
for(i in 1:22){
  fp[[i]]=rep(0, 10)
  fn[[i]]=rep(0, 10)
  error[[i]]=rep(0, 10)
  mcc[[i]]=rep(0, 10)
}

set.seed(999)
seed=sample(1000, 10)
class=list()

for(idx in 1:10){
  set.seed(seed[idx])
  train=sample(dim(iono)[1], 0.7*dim(iono)[1])
  x=iono[train,3:34]
  test=iono[-train,3:34]
  
  L.true=as.numeric(iono$V35[train])-1
  
  n=length(L.true)
  
  L.org=L.true
  L.org[sample(n, size=n*0.85)]=NA
  sum(L.org==0, na.rm = T)
  sum(L.org==1, na.rm = T)
  
  sum(L.org==0, na.rm=T)
  sum(L.org==1, na.rm=T)
  L.org=as.factor(L.org)
  df=data.frame(x, L.org)
  
  L.test=as.numeric(iono$V35[-train])-1
  df.test=data.frame(test, L.test)
  
  
  param=list()
  p=dim(x)[2]
  midpoint=rep(0, 15)
  
  class[[1]] <- EntropyRegularizedLogisticRegression(L.org~.,df, lambda=0.01,lambda_entropy = 100)
  class[[2]] <- EMLeastSquaresClassifier(L.org~.,df, max_iter = 5000)
  # class[[3]] <- EMLinearDiscriminantClassifier(L.org~., df, max_iter = 5000)
  class[[4]] <- ICLeastSquaresClassifier(L.org~., df)#################
  # class[[5]] <- ICLinearDiscriminantClassifier(L.org~., df)
  class[[6]] <- LaplacianSVM(L.org~., df) ########################
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
  class[[20]] <- WellSVM(L.org~., df) ###################
  class[[21]] <- svmlin(L.org~., df)
  class[[22]] <- EMNearestMeanClassifier(L.org~.,df)#############
  # class[[23]] <- KernelLeastSquaresClassifier(L.org~., df)
  # class[[24]] <- LaplacianKernelLeastSquaresClassifier(L.org~., df)
  # class[[25]] <- LinearTSVM(L.org~., df)
  # class[[26]] <- S4VM(L.org~., df)
  # class[[27]] <- SelfLearning(L.org~., df)
  # class[[28]] <- TSVM(L.org~., df)
  irange=c(4,6,20:22)
  #       
  for(i in irange){
    fp[[i]][idx]=sum(predict(class[[i]],df.test)[which(L.test==0)]==1)/sum(L.test==0)
    print(fp[[i]][idx])
    fn[[i]][idx]=sum(predict(class[[i]],df.test)[which(L.test==1)]==0)/sum(L.test==1)
    print(fn[[i]][idx])
    count=table(L.test)
    fpc=fp[[i]][idx]*count[1]
    fnc=fn[[i]][idx]*count[2]
    tpc=count[2]-fnc
    tnc=count[1]-fpc
    mcc[[i]][idx]=(tpc*tnc-fpc*fnc)/sqrt((tpc+fpc)*(tpc+fnc)*(tnc+fpc)*(tnc+fnc))
    error[[i]][idx]=(fpc+fnc)/(count[1]+count[2])
  }

}

for(i in c(4,6,20:22)){
  print(i)
  print(mean(fp[[i]]))
  print(mean(fn[[i]]))
  print(mean(error[[i]]))
  print(mean(mcc[[i]]))
}
