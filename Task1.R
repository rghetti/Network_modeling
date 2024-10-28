#importing libraries

library(sna)
library(network)


#importing the data

setwd("./Dataset")
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


g1.5=netlogit(f.matrix,ad.matrix, rep=permutations,nullhyp='qapy')
g1.5$names=c('intercept', 'advice')
summary(g1.5)
#at the significance level of 0.05 we can see there is statistical 
#evidence that friendship is a good predictor of advice seeking


#extract same_department as a bool matrix (1 is departments are the same and 0 otherwise)

num_people = 21

same_department.matrix=matrix(0, nrow=num_people, ncol=num_people)
for (i in 1:num_people){
  for (j in 1:num_people){
    if(attributes[["nodeDepartment"]][i] == attributes[["nodeDepartment"]][j])
      same_department.matrix[i, j] = 1
  }
}

#make a matrix for testing that senior managers are less likely to nominate friends. The value at position (i, j) is tenure[i]

tenure.matrix = matrix(0, nrow=num_people, ncol=num_people)
for(i in 1:num_people){
  for(j in 1:num_people){
    tenure.matrix[i, j] = attributes[["nodeTenure"]][i] 
  }
}

#make a matrix for testing that a friendship nomination is more likely between a pair of managers of a similar age. The value at position (i, j) is |age[i] - age[j]|
    
age_dif_matrix = matrix(0, nrow=num_people, ncol=num_people)
for(i in 1:num_people){
  for(j in 1:num_people){
    age_dif_matrix[i, j] = abs(attributes[["nodeAge"]][i]-attributes[["nodeAge"]][j])
  }
}

#build the QAP test with additional features

zm <- list(ad.matrix, same_department.matrix, tenure.matrix, age_dif_matrix)

g2<-netlogit(f.matrix,zm, rep=permutations,nullhyp='qapspp')
g2$names<-c('intercept', 'advice', 'same_department', 'tenure', 'age_dif')
summary(g2)

#at the significance level of 0.05 we can see there is statistical 
#evidence that people more likely nominate their friends peoples from the same 
#department or people, who they are asking an advice.
#From the other hand, at the significance level of 0.05, tenure and 
#age difference do not correlate with friendship nomination


#extract reporting as adjacency matrix
rep = nodes[293:312,2:3]
reports.matrix=matrix(0,nrow=num_people,ncol=num_people)
colnames(reports.matrix)=1:ncol(reports.matrix)
for(k in 1:(312-293+1)){
  reports.matrix[rep[k,1],rep[k,2]]=1
}

#Test that if a person reports to another person as part of his function
#in organization, then the first person more likely nominate the second
#as a friend 

zm1 <- list(ad.matrix, same_department.matrix, tenure.matrix, age_dif_matrix, 
            reports.matrix)
g3<-netlogit(f.matrix,zm1, rep=permutations,nullhyp='qapspp')
g3$names<-c('intercept', 'advice', 'same_department', 'tenure', 'age_dif', 'report_matrix')
summary(g3)

#At the significance level of 0.05 we can see that there is not correlation 
#between friendship and reporting one person to be a part of others functions 
#in the organization
