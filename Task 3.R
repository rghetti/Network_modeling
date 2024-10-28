#importing libraries
install.packages("ergm")
library(sna)
library(network)
library('ergm')
library(igraph)
library(RColorBrewer)


#importing the data

setwd("./Dataset")
list.files()

attributes=read.table("Krackhardt-High-Tech_nodes.txt",header=TRUE)
nodes=as.matrix(read.table("Krackhardt-High-Tech_multiplex.edges",header=FALSE))
colnames(nodes)=c("layerID","IDi","IDj","weight")
set.seed(0)


#extract friendship as an adjacency matrix

friendship=nodes[191:292,2:4]
f.matrix=matrix(0,nrow=21,ncol=21)
colnames(f.matrix)=1:ncol(f.matrix)
for( k in 1:102){
  f.matrix[friendship[k,1],friendship[k,2]]=1
}

# Creating a network object 
netw <- network(f.matrix,directed=TRUE)

# Adding attributes
netw %v% "age" <- attributes$nodeAge
netw %v% "tenure" <- attributes$nodeTenure
netw %v% "level" <- attributes$nodeLevel
netw %v% "department" <- attributes$nodeDepartment

netw
# Visualize the network
plot(netw)

model0 <- ergm(netw ~ edges + nodecov('department'))
summary(model0)
coef(model0)
theta1=coef(model0)[1]
theta2=coef(model0)[2]

prob = exp(theta1 + theta2) / (1 + exp(theta1 + theta2))
prob_having_tie = exp(theta1) / (1 + exp(theta1))

#probability over tie belonging to the same department is 0.2 that is slightly more
# than probability of having a tie in general (0.188).

model1 <- ergm(netw ~ edges + nodecov('department') + mutual() +  gwesp(decay=0.3,fixed=TRUE) + idegree(6))
summary(model1)

#Under the significance of 0.05 we can see that popularity measured as number 
#of nodes with in degree 6 has no influence on the model,
#as opposed to transitivity and reciprocity which have a positive influence on the model.

#The algorithm experiences near-degeneracy when evaluated with ttriple statics, so we used the gwesp
#to resolve the problem with the decay of 0.3.

mcmc.diagnostics(model1)
gof_results = gof(model1)
par(mfrow=c(2,3),mar=c(5, 4, 4, 2))
plot(gof_results)

#The out degree distribution isn't accurately describing
#the real network. The in degree fits the network faithfully
#most of the time. The edge-wise shred partners does not accurately describe 
#the network for values under 5. The same behavior is seen in the distribution 
#of the minimum geodesic distance 

