library(blockmodels)
library(blockcluster)
library(mclust)
library(FactoMineR)

set.seed(72)

#################
#####Binary######
#################

data(binarydata)  
head(binarydata)
dim(binarydata)
heatmap(binarydata)
x=as.data.frame(binarydata)
for(i in 1:ncol(x)){
  x[,i]=as.factor(x[,i])
}
y=MCA(x)

###Estimation best model for blockclusters
strategy=coclusterStrategy()
summary(strategy)
vect=c()
keep=-10000000000
k.max=4
start_time <- Sys.time()
for(i in 2:k.max){
  for(j in 2:k.max){
    out<-coclusterBinary(binarydata,nbcocluster=c(i,j))
    print(cat("i vaut: ", i," j vaut: ",j,""))
    vect=cbind(vect,slot(out,"ICLvalue"))
    print(slot(out,"ICLvalue"))
    if(slot(out,"ICLvalue") > keep){
      keep = slot(out,"ICLvalue")
      keeprow=i
      keepcol=j
    }
  }
}
end_time <- Sys.time();end_time-start_time
plot(1:length(vect), vect,
     type="b", pch = 19, frame = FALSE,
     xaxt = "n",
     xlab="Number of clusters K(row,col)",
     ylab="Value of ICL")
axis(1, at=1:length(vect), labels=c("2,2","2,3","2,4","3,2","3,3","3,4","4,2","4,3","4,4"))

plot(cat("best params: row =",keeprow," col =",keepcol))

###blockcluster binary
out<-coclusterBinary(binarydata,nbcocluster=c(keeprow,keepcol))
summary(out)
dev.off()
plot(out)
plot(out, type = "distribution")
slotNames(out)
slot(out,"rowclass")
slot(out,"colclass")
slot(out,"rowproportions")
slot(out,"classdispersion")

#blockmodel Binary
my_model <- BM_bernoulli("LBM",binarydata)
start_time <- Sys.time()
my_model$estimate()
end_time <- Sys.time();end_time-start_time
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(which.max(my_model$ICL))
my_model$memberships[[which.max(my_model$ICL)]]$Z1
my_model$memberships[[which.max(my_model$ICL)]]$Z2

###Comparaison
resrow=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z1)[1]){
  resrow=cbind(resrow,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z1[i,]))
}
table(slot(out,"rowclass"),resrow[1,])

rescol=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z2)[1]){
  rescol=cbind(rescol,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z2[i,]))
}
table(slot(out,"colclass"),rescol[1,])


#################
###Contingency###
#################
data(contingencydataunknown)
head(contingencydataunknown)
dim(contingencydataunknown)
heatmap(contingencydataunknown)
x=as.data.frame(contingencydataunknown)
for(i in 1:ncol(x)){
  x[,i]=as.factor(x[,i])
}
y=MCA(x)

###Estimation best model for blockclusters
strategy=coclusterStrategy()
strategy = coclusterStrategy( nbinititerations = 5, nbxem = 2, nbiterations_int = 2
                              , nbiterationsxem = 10, nbiterationsXEM = 100, epsilonXEM=1e-5)
summary(strategy)
vect=c()
keep=10000000000
k.max=4
start_time <- Sys.time()
for(i in 2:k.max){
  for(j in 2:k.max){
    out<-coclusterContingency(contingencydataunknown,nbcocluster=c(i,j))
    print(cat("i vaut: ", i," j vaut: ",j,""))
    vect=cbind(vect,slot(out,"ICLvalue"))
    print(slot(out,"ICLvalue"))
    if(slot(out,"ICLvalue") < keep){
      keep = slot(out,"ICLvalue")
      keeprow=i
      keepcol=j
    }
  }
}
end_time <- Sys.time();end_time-start_time
plot(1:length(vect), vect,
     type="b", pch = 19, frame = FALSE,
     xaxt = "n",
     xlab="Number of clusters K(row,col)",
     ylab="Value of ICL")
axis(1, at=1:length(vect), labels=c("2,2","2,3","2,4","3,2","3,3","3,4","4,2","4,3","4,4"))

plot(cat("best params: row =",keeprow," col =",keepcol))

###blockcluster contingency
out<-coclusterContingency(contingencydataunknown,nbcocluster=c(2,2),strategy=strategy)
out<-coclusterContingency(contingencydataunknown,nbcocluster=c(keeprow,keepcol),strategy=strategy)
summary(out)
dev.off()
plot(out)
plot(out, type = "distribution")
slotNames(out)
slot(out,"rowclass")
slot(out,"colclass")
slot(out,"rowproportions")
slot(out,"columnproportions")
slot(out,"classgamma")


my_model <- BM_poisson("LBM",contingencydataunknown,exploration_factor=10,explore_max=5)
my_model <- BM_poisson("LBM",contingencydataunknown)
start_time <- Sys.time()
my_model$estimate()
end_time <- Sys.time();end_time-start_time
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(which.max(my_model$ICL))
my_model$memberships[[which.max(my_model$ICL)]]$Z1
my_model$memberships[[which.max(my_model$ICL)]]$Z2

###Comparaison
resrow=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z1)[1]){
  resrow=cbind(resrow,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z1[i,]))
}
table(slot(out,"rowclass"),resrow[1,])

rescol=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z2)[1]){
  rescol=cbind(rescol,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z2[i,]))
}
table(slot(out,"colclass"),rescol[1,])

#################
###Continuous###
#################
data(gaussiandata)
res.pca=PCA(gaussiandata)
heatmap(gaussiandata)

summary(strategy)
vect=c()
keep=-10000000000
k.max=4
start_time <- Sys.time()
for(i in 2:k.max){
  for(j in 2:k.max){
    out<-coclusterContinuous(gaussiandata,nbcocluster=c(i,j))
    print(cat("i vaut: ", i," j vaut: ",j,""))
    vect=cbind(vect,slot(out,"ICLvalue"))
    print(abs(slot(out,"ICLvalue")))
    if(slot(out,"ICLvalue") > keep){
      print("oui")
      keep = slot(out,"ICLvalue")
      keeprow=i
      keepcol=j
    }
  }
}
end_time <- Sys.time();end_time-start_time
plot(1:length(vect), abs(vect),
     type="b", pch = 19, frame = FALSE,
     xaxt = "n",
     xlab="Number of clusters K(row,col)",
     ylab="Value of ICL")
axis(1, at=1:length(vect), labels=c("2,2","2,3","2,4","3,2","3,3","3,4","4,2","4,3","4,4"))

plot(cat("best params: row =",keeprow," col =",keepcol))

###blockcluster Continuous
out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3))
out<-coclusterContinuous(gaussiandata,nbcocluster=c(keeprow,keepcol))
summary(out)
dev.off()
plot(out)
plot(out, type = "distribution")
slotNames(out)
slot(out,"rowclass")
slot(out,"colclass")
slot(out,"rowproportions")


my_model <- BM_gaussian("LBM",gaussiandata)

start_time <- Sys.time()
my_model$estimate()
end_time <- Sys.time();end_time-start_time

end_time <- Sys.time();end_time-start_time
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(which.max(my_model$ICL))
my_model$memberships[[which.max(my_model$ICL)]]$Z1
my_model$memberships[[which.max(my_model$ICL)]]$Z2

resrow=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z1)[1]){
  resrow=cbind(resrow,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z1[i,]))
}
table(slot(out,"rowclass"),resrow[1,])

rescol=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z2)[1]){
  rescol=cbind(rescol,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z2[i,]))
}
table(slot(out,"colclass"),rescol[1,])

#################################
#############Donn??es R??elles#####
#################################
library("png")
image<- readPNG("/Users/jzk/Documents/M2/coclust/29.png");dim(image)
summary(strategy)
vect=c()
keep=10000000000
k.max=4
for(i in 2:k.max){
  for(j in 2:k.max){
    out<-coclusterContinuous(image,nbcocluster=c(i,j))
    print(cat("i vaut: ", i," j vaut: ",j,""))
    vect=cbind(vect,slot(out,"ICLvalue"))
    print(abs(slot(out,"ICLvalue")))
    if(slot(out,"ICLvalue") < keep){
      print("oui")
      keep = slot(out,"ICLvalue")
      keeprow=i
      keepcol=j
    }
  }
}
plot(1:length(vect), abs(vect),
     type="b", pch = 19, frame = FALSE,
     xaxt = "n",
     xlab="Number of clusters K(row,col)",
     ylab="Value of ICL")
axis(1, at=1:length(vect), labels=c("2,2","2,3","2,4","3,2","3,3","3,4","4,2","4,3","4,4"))

plot(cat("best params: row =",keeprow," col =",keepcol))
out<-coclusterContinuous(image,nbcocluster=c(keeprow,keepcol))
summary(out)
dev.off()
plot(out)
plot(out, type = "distribution")
slotNames(out)
slot(out,"rowclass")
slot(out,"colclass")
slot(out,"rowproportions")

my_model <- BM_gaussian("LBM",image, explore_max=6)
start_time <- Sys.time()
my_model$estimate()
end_time <- Sys.time();end_time-start_time
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(which.max(my_model$ICL))
my_model$memberships[[which.max(my_model$ICL)]]$Z1
my_model$memberships[[which.max(my_model$ICL)]]$Z2

resrow=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z1)[1]){
  resrow=cbind(resrow,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z1[i,]))
}
table(slot(out,"rowclass"),resrow[1,])

rescol=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z2)[1]){
  rescol=cbind(rescol,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z2[i,]))
}
table(slot(out,"colclass"),rescol[1,])

#################################
############Text Mining##########
#################################
library(R.matlab)
library(readr)
library(Matrix)
library(NMF)
library(tidytext)
library(tm)
library(slam)
library(dplyr)
library(SnowballC)
library(skmeans)
library(textir)
library(stm)
library(factoextra)
library(foreach)
library(doParallel)
library(fastICA)
library(wordcloud)
library(topicmodels)

brute <- read.table("/Users/jzk/Downloads/classicdocspreprocessed/docbyterm.txt",header=TRUE)
typeof(brute)
brute

brute$X1=brute$X7095
brute$X2=brute$X5896
brute$X3=brute$X247158
brute$X3[is.na(brute$X3)] <- 1
x1 <- as.factor(brute$X1)
levels(x1) <- 1:length(levels(x1))
x1 <- as.numeric(x1)

x2 <- as.factor(brute$X2)
levels(x2) <- 1:length(levels(x2))
x2 <- as.numeric(x2)

data_stm=simple_triplet_matrix(x1, x2, brute$X3, nrow = max(x1), ncol = max(x2))
simple_triplet_matrix_sparse=data_stm
data <-  sparseMatrix(i=simple_triplet_matrix_sparse$i, j=simple_triplet_matrix_sparse$j, x=simple_triplet_matrix_sparse$v,
                                              dims=c(simple_triplet_matrix_sparse$nrow, simple_triplet_matrix_sparse$ncol))
data=as.matrix(data)
max(data)
unique(data)
x=data
for(i in 1:ncol(x)){
  x[,i]=as.factor(x[,i])
}
y=MCA(x)
View(data)
out<-coclusterCategorical(head(data,nbcocluster=c(4,4)))
dev.new(width=25000, height=25000)
png(filename = "/Users/jzk/Documents/M2/coclust/Classic4cocl.png", width = 450, height = 450,
    pointsize = 12, bg = "white",  res = NA)
plot(out)
dev.off()


#blockmodel Binary
my_model <- BM_poisson("LBM",data,explore_min=6, explore_max=6)
my_model$estimate()
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(which.max(my_model$ICL))
my_model$memberships[[which.max(my_model$ICL)]]$Z1
my_model$memberships[[which.max(my_model$ICL)]]$Z2

########Donn??es r??elles Binary######
###################################
library(miic)
data("hematoData")
df=hematoData
View(df)
dim(df)
colnames(df)
rownames(df)

sparsity <-function(df){
  return (sum(df == 0)/(dim(df)[1]*dim(df)[2]))
}

sparsity(df)
df2=as.matrix(df)
binarydata=df2


###Estimation best model for blockclusters
strategy=coclusterStrategy()
summary(strategy)
vect=c()
keep=-10000000000
k.max=4
for(i in 2:k.max){
  for(j in 2:k.max){
    out<-coclusterCategorical(binarydata,nbcocluster=c(i,j))
    print(cat("i vaut: ", i," j vaut: ",j,""))
    vect=cbind(vect,slot(out,"ICLvalue"))
    print(slot(out,"ICLvalue"))
    if(slot(out,"ICLvalue") > keep){
      keep = slot(out,"ICLvalue")
      keeprow=i
      keepcol=j
    }
  }
}
plot(1:length(vect), vect,
     type="b", pch = 19, frame = FALSE,
     xaxt = "n",
     xlab="Number of clusters K(row,col)",
     ylab="Value of ICL")
axis(1, at=1:length(vect), labels=c("2,2","2,3","2,4","3,2","3,3","3,4","4,2","4,3","4,4"))

plot(cat("best params: row =",keeprow," col =",keepcol))

###blockcluster binary
out<-coclusterBinary(binarydata,nbcocluster=c(3,5))
summary(out)
dev.off()
plot(out)
plot(out, type = "distribution")
slotNames(out)
slot(out,"rowclass")
slot(out,"colclass")
slot(out,"rowproportions")
slot(out,"classdispersion")

#blockmodel Binary
my_model <- BM_bernoulli("LBM",binarydata,explore_max=8)
my_model$estimate()

which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(which.max(my_model$ICL))
my_model$memberships[[which.max(my_model$ICL)]]$Z1
my_model$memberships[[which.max(my_model$ICL)]]$Z2

###Comparaison
resrow=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z1)[1]){
  resrow=cbind(resrow,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z1[i,]))
}
table(slot(out,"rowclass"),resrow[1,])

rescol=c()
for(i in 1:dim(my_model$memberships[[which.max(my_model$ICL)]]$Z2)[1]){
  rescol=cbind(rescol,which.max(my_model$memberships[[which.max(my_model$ICL)]]$Z2[i,]))
}
table(slot(out,"colclass"),rescol[1,])


#################################
####Comparaison avec double GMM:#
#################################
start_time <- Sys.time()
clusterrow=Mclust(gaussiandata)
clustercol=Mclust(t(gaussiandata))
end_time <- Sys.time()
end_time - start_time



table(clusterrow$classification,resrow[1,])
table(clustercol$classification,rescol[1,])

x=kmeans(gaussiandata,2)
table(x$cluster,resrow[1,])

y=kmeans(t(gaussiandata),3)
table(y$cluster,rescol[1,])
