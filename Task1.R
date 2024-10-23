#importing libraries

library(sna)
library(network)



#importing the data

setwd("C:/Users/Riccardo Ghetti/Downloads/Network modeling/assignment1/assignment1/Krackhardt-High-Tech_Multiplex_Social/Dataset")
list.files()
attributes=read.table("Krackhardt-High-Tech_nodes.txt",header=TRUE)
nodes=as.matrix(read.table("Krackhardt-High-Tech_multiplex.edges",header=FALSE))
colnames(nodes)=c("layerID","IDi","IDj","weight")
set.seed(0)
permutations=10000

#extract friendship as an adjacency matrix

friendship=nodes[191:292,2:4]
f.matrix=matrix(0,nrow=21,ncol=21)
colnames(f.matrix)=1:ncol(f.matrix)
for( k in 1:102){
  f.matrix[friendship[k,1],friendship[k,2]]=1
}

#extract advice as an adjacency matrix

advice=nodes[1:190,2:4]
ad.matrix=matrix(0,nrow=21,ncol=21)
colnames(ad.matrix)=1:ncol(ad.matrix)
for (k in 1:190){
  ad.matrix[advice[k,1],advice[k,2]]=1
}

#build the QAP test

g1=netlogit(ad.matrix,f.matrix, rep=permutations,nullhyp='qapy')
g1$names=c('intercept','friendship')
summary(g1)
#at the significance level of 0.05 we can see there is statistical 
#evidence that friendship is a good predictor of advice seeking

