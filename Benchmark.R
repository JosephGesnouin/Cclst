library(blockmodels)
library(blockcluster)
## generation of one SBM network
npc <- 30 # nodes per class
Q <- 5 # classes
n <- npc * Q # nodes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
M<-matrix(rnorm(n*n,sd=10),n,n)+Z%*%Mu%*%t(Z) ## adjacency matrix

dim(M)
my_model <- BM_gaussian("SBM",M )
my_model$estimate()
which.max(my_model$ICL)

x=Mclust(binarydata)
x$classification

y=Mclust(t(binarydata))
y$classification

####Test cosmicCancer
my_model <- BM_gaussian("SBM",adj2)
my_model$estimate()
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(16)


keep=-10000000000
for(i in 2:4){
  for(j in 2:4){
    out<-coclusterBinary(binarydata,nbcocluster=c(i,j))
    #print("i vaut: " i + " j vaut: "+ j)
    print(slot(out,"ICLvalue"))
    if(slot(out,"ICLvalue") > keep){
      keep = abs(slot(out,"ICLvalue"))
      keeprow=i
      keepcol=j
    }
  }
}
keeprow
keepcol
slotNames(out)
slot(out,"rowclass")
out<-coclusterBinary(binarydata,nbcocluster=c(4,4))
summary(out)
dev.off()
plot(out)

###############  Binary
data(binarydata)  
head(binarydata)
dim(binarydata)
out<-coclusterBinary(binarydata,nbcocluster=c(2,3))
summary(out)
plot(out)
plot(out, type = "distribution")

my_model <- BM_bernoulli("LBM",binarydata)
my_model$estimate()
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]
my_model$plot_obs_pred(which.max(my_model$ICL))
my_model$memberships[[which.max(my_model$ICL)]]

#####Contingency
data(contingencydataunknown)
head(contingencydataunknown)
dim(contingencydataunknown)
strategy = coclusterStrategy( nbinititerations = 5, nbxem = 2, nbiterations_int = 2
                              , nbiterationsxem = 10, nbiterationsXEM = 100, epsilonXEM=1e-5)
out<-coclusterContingency( contingencydataunknown, nbcocluster=c(2,3), strategy = strategy)
summary(out)
dev.off()
plot(out)

my_model <- BM_poisson("LBM",contingencydataunknown)
my_model$estimate()
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]


#####Continuous
my_model <- BM_poisson("LBM",contingencydataunknown)
my_model$estimate()
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)][1]


####
data(gaussiandata)
out<-coclusterContinuous(gaussiandata,nbcocluster=c(2,3))
summary(out)
dev.off()
plot(out)

my_model <- BM_gaussian("LBM",gaussiandata)
dev.off()
my_model$estimate()
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)]


#####MNIST

my_model <- BM_gaussian("LBM",test)
dev.off()
my_model$estimate()
which.max(my_model$ICL)
my_model$memberships[which.max(my_model$ICL)]
dim(test)
out<-coclusterContinuous(test,nbcocluster=c(10,))
summary(out)
plot(out)
