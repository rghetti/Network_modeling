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
permutations=10000


friendship=nodes[191:292,2:4]
f.matrix=matrix(0,nrow=21,ncol=21)
colnames(f.matrix)=1:ncol(f.matrix)
for( k in 1:102){
  f.matrix[friendship[k,1],friendship[k,2]]=1
}

advice=nodes[1:190,2:4]
ad.matrix=matrix(0,nrow=21,ncol=21)
colnames(ad.matrix)=1:ncol(ad.matrix)
for (k in 1:190){
  ad.matrix[advice[k,1],advice[k,2]]=1
}

num_people = 21

tenure.matrix = matrix(0, nrow=num_people, ncol=num_people)
for(i in 1:num_people){
  for(j in 1:num_people){
    tenure.matrix[i, j] = attributes[["nodeTenure"]][i] 
  }
}


age_dif_matrix = matrix(0, nrow=num_people, ncol=num_people)
for(i in 1:num_people){
  for(j in 1:num_people){
    age_dif_matrix[i, j] = abs(attributes[["nodeAge"]][i]-attributes[["nodeAge"]][j])
  }
}


# Creating a network object 
netw <- network(f.matrix,directed=TRUE)

# Adding attributes
netw %v% "age" <- attributes$nodeAge
netw %v% "tenure" <- attributes$nodeTenure
netw %v% "level" <- attributes$nodeLevel
netw %v% "department" <- attributes$nodeDepartment
netw %e% "advice" <-ad.matrix
netw %e% "tenure" <- tenure.matrix

netw

m1=ergm(netw ~ edgecov(ad.matrix) + nodecov('department') + 
          edgecov(tenure.matrix) + edgecov(age_dif_matrix))
summary(m1)

#As in the MR-QAP test, there is evidence that advice-seeking and friendship correlate positively.
#The results about tenure are exactly the same, while homophily of department has a negative influence.
#Finally, in this model the age difference impact as a relevant negative predictor.

m2 = ergm(netw ~ edgecov(ad.matrix) + nodecov('department') +  edgecov(tenure.matrix)
          + edgecov(age_dif_matrix) + edges  + mutual() +  
            gwesp(decay=0.3,fixed=TRUE) + idegree(6))
summary(m2)

#After adding the structural terms, homophily of department and age difference become
#statistically  irrelevant. 
#Advice seeking still has a positive impact on the friendship, while tenure still is
#statistically irrelevant.
