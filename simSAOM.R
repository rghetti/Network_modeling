# Network Modeling - HS 2024
# C. Stadtfeld, A. Uzaheta and I. Smokovic
# Social Networks Lab
# Department of Humanities, Social and Political Sciences
# ETH Zurich
# 18 November 2024
#
# Assignment 2 - Task 3



# Task 3.1 ----------------------------------------------------------------
# The function "simulation" simulates the network evolution between 
# two time points. 
# Given the network at time t1, denoted by x1, the function simulates the 
# steps of the continuous-time Markov chain defined by a SAOM with outdegree,
# recip and dyadic covariate (W matrix) statistics.
# Unconditional simulation is used.
# The function returns the network at time t2.
# The structure of the algorithm is described in the file
# _Simulating from SAOM.pdf_ available in
# the Lecture notes and additional material section on Moodle.

#' Simulate the network evolution between two time points
#'
#' @param n number of actors in the network
#' @param x1 network at time t1
#' @param W matrix of dyadic covariate
#' @param lambda rate parameter
#' @param beta1 outdegree parameter
#' @param beta2 reciprocity parameter
#' @param beta3 dyadic covariate parameter
#'
#' @return network at time t2
#'
#' @examples
#' netT1 <- matrix(c(
#'   0, 1, 0, 0, 0,
#'   0, 0, 0, 1, 0,
#'   0, 0, 0, 1, 1,
#'   1, 0, 1, 0, 0,
#'   0, 1, 1, 0, 1
#'   ), 
#'   nrow = 5, ncol =  5, byrow = TRUE)
#' W <- matrix(0, nrow = 5, ncol = 5)
#' values <- runif(5 * 2, 0, 1)
#' W[upper.tri(W)] <- values
#' W <- W + t(W)
#' netT2 <- simulation(5, netT1, W, 4, -2, 0.5, 1)
#' 
#' 
library(RSiena)
library(sna)
library(parallel)
library(MASS)

comp_stat<-function(x1,W,i,mean_W,n){
  a=0
  b=0
  d=0
  for (k in 1:n){
    a=a+x1[i,k]  #outdegree
    b=b+x1[i,k]*x1[k,i]  #reciprocity
    d=d+x1[i,k]*(W[i,k]-mean_W)           #dyadic covariate centered effect
  }
  return(as.double(c(a,b,d)))
}

comp_prob<-function(b1,b2,b3,vect_stat,n){
  vect_prob=matrix(0,1,n)
  sum=0
  for (k in 1:n)  {   #compute the normalizer constant and save the weights
    p=exp(b1*vect_stat[1,k]+b2*vect_stat[2,k]+b3*vect_stat[3,k]) 
    sum=sum+p
    vect_prob[k]=p
    }
    
  for (k in 1:n)     #normalize the probabilities
    vect_prob[k]=vect_prob[k]/sum
  return(vect_prob)
}

simulation <- function(n, x1, W, lambda, beta1, beta2, beta3) {
  t <- 0 # time
  x <- x1
  mean_W=mean(matrix(W,1,n*n)) #compute the mean of the vectorized matrix
  while (t < 1) {
    dt <- rexp(1, n * lambda)
    # --- MISSING ---
    i=sample(1:n,size=1) #draw a random actor with uniform prob from the n available 
    vect_stat=matrix(0,3,n)  #initialize the container of all statistics
    
    for (k in 1:n){    #computing the statistics for all choices
      
      if(k==i){                                #if k==i we interpret it as the choice of leaving the network unchanged 
        vect_stat[1:3,k]=comp_stat(x1,W,k,mean_W,n)     #Break is necessary to stop the computation in this iteration
        break                                   
      }                                        
      
      x2=x                                     #copy the current network for computing statistics of its variation
      x2[i,k]=xor(x2[i,k],1)                   #invert the tie i->j in the hypothetical network
      vect_stat[1:3,k]=comp_stat(x2,W,k,mean_W,n)       #compute and store the statistics for the current choice
      
    }
    
    vect_prob=comp_prob(beta1,beta2,beta3,vect_stat,n)    #compute the multinomial probabilities for every choice
    j=sample(x=1:n,size=1,prob=vect_prob)                 #sample the actor according to the computed probabilities
    
    if (i!=j)    #invert the tie i->j in x else do nothing
      x[i,j]=xor(x[i,j],1)
    t=t+dt
   }
 
  return(x)
}


# Task 3.2 ----------------------------------------------------------------
# Consider the two adjacency matrices in the files net1.csv and net2.csv.
# Estimate the parameters of the SAOM with outdegree, reciprocity and
# dyadic covariate statistics using the function `siena07`.
# You can extract the estimated parameters from the `rate` and `theta`
#  components of the output object (e.g., res$rate and res$theta).

# ---MISSING---
setwd("C:/Users/Riccardo Ghetti/Downloads/Network modeling/Assignment2/Assignment2")
net1=as.matrix(read.csv("net1.csv",header=FALSE))
net2=as.matrix(read.csv("net2.csv",header=FALSE))
W=as.matrix(read.csv("W.csv",header=FALSE))

friendship=sienaDependent(array(c(net1,net2), dim=c(22,22,2)))

dyad=coDyadCovar(W)
mydata=sienaDataCreate(friendship,dyad)
mydata

myeff=getEffects(mydata)
myeff
effectsDocumentation(myeff)

myeff=includeEffects(myeff,X,interaction1 = "dyad")

myAlgorithm <- sienaAlgorithmCreate(
  projname = "friends_W",
  nsub = 4, n3 = 3000, seed = 1908
)

model0 <- siena07(
  myAlgorithm,
  data = mydata, effects = myeff,
  returnDeps = TRUE,
  useCluster = TRUE, nbrNodes = 4, batch = FALSE
)
rate=model0$rate
betas=model0$theta
rm(model0)
# Task 3.3 ----------------------------------------------------------------
# Conditioning on the first observation, generate 1,000 simulations of the 
# network evolution
# Compute the triad census counts for each simulated network.
# Save the results in an object, named `triadCensus`, in which rows are
# the index of a simulated network and columns are the type of triads.
# Column names should use the triad type name, e.g., "003", "012", "102", ... 

results=array(0,dim=c(1000,22,22))
beta1=betas[1]
beta2=betas[2]
beta3=betas[3]

for (i in 1:1000){
  results[i,1:22,1:22]=simulation(22,net1,W,rate,beta1,beta2,beta3)
}
triadCensus=triad.census(results)

# Task 3.4 ----------------------------------------------------------------
## i. standardized the simulated network stats. ----
##   Name the resulting object as triadCensusStd

triad_stat=matrix(0,16,2)

for (e in 1:16){
  triad_stat[e,1]=mean(triadCensus[1:1000,e]) #the mean of each column
  triad_stat[e,2]=sd(triadCensus[1:1000,e]) #its standard deviation
}

triadCensusStd=matrix(0,1000,16)

for(e in 1:16){
  average_e=triad_stat[e,1]
  sd_e=triad_stat[e,2]
  for (g in 1:1000){
    triadCensusStd[g,e]=(triadCensus[g,e]-average_e)/sd_e  #standardize the values for each network
  }
}




## ii. variance-covariance matrix and its generalized inverse.         ----

p_h=cov(triadCensusStd,triadCensusStd)  #compute the covariance matrix
p_inv=ginv(p_h)                         #compute its generalized inverse
## iii. standardized the observed values of the triad census counts    ----
##  in the second observation using values from i.

net2

Obs_stat=triad.census(net2)
for (e in 1:16){
  Obs_stat[e]=(Obs_stat[e]-triad_stat[e,1])/triad_stat[e,2] #normalize obs stat
}

colnames(triadCensusStd)<-colnames(triadCensus)

## iv. Monte-Carlo Mahalanobis distance computation                                ----
# Compute the Mahalanobis distance using the mhd function for 
# the auxiliar statistics of the simulated networks and the observed network.
# Remember to drop statistics with variance of 0 for the plot and
# Mahalanobis distance computation, report which statistics suffer this issue.

#' Compute the Mahalanobis distance
#'
#' @param auxStats numerical vector with the mean centered or standardized
#'   auxiliar statistics
#' @param invCov numerical matrix with the inverse of the variance-covariance
#'   matrix of the auxiliar statistics in the simulated networks
#'
#' @return numeric value with the Mahalanobis distance of auxiliar stats
#'
#' @examples
#' mhd(c(2, 4) - c(1.5, 2), solve(matrix(c(1, 0.8, 0.8, 1), ncol = 2)))
mhd <- function(auxStats, invCov) {
  t(auxStats) %*% invCov %*% auxStats
  }

mahal=matrix(0,1001,1)

for(e in 1:1000){
  mahal[e]=mhd(triadCensusStd[e,1:16], p_inv)
}

mahal[1001]=mhd(t(Obs_stat), p_inv)




## v. Monte-Carlo p-value computation                                ----
# Compute the proportion of simulated networks where the distance is 
# equal or greater than the distance in the observed network.

p_val=0
for (e in 1:1000){
  if(mahal[e]>=mahal[1001])
    p_val=p_val+1
}
p_val=p_val/1000
p_val

# violin plots ------------------------------------------------------------
# Fill out the missing part and run the code to obtain the violin plots

# install.packages(c("tidyverse", "ggplot2"))  # # run this line to install 
library(tidyverse)  # used: dplyr and tidyr
library(ggplot2)

# Given the array triadCensusStd, create a data frame from it in a long format, 
# do the same for the observed network statistics at time t2.
# Named the data frame "triadCensusDf" and "triadCensusObs".
# Drops statistics with variance of 0 for the plot.

triadCensusDf <- data.frame(triadCensusStd) |>
  select(where(~ var(.) > 0)) |>  # Drop statistics with zero variance
  pivot_longer(
    everything(),
    names_to = "triad", names_pattern = "^X(.+)$",
    values_to = "nnodes"
  )

# Compute the statistics of the observed network at time t2,
#  standardized using the stats from 2.4 literal i.
obs=data.frame(Obs_stat)
triadCensusObs <- obs |>
  pivot_longer(
    everything(),
    names_to = "triad", names_pattern = "^X(.+)$",
    values_to = "nnodes"
  ) |>
  filter(triad %in% unique(triadCensusDf$triad))

# The following code computes the 5% and the 95% quantiles
# of the triad counts by type
percTriad <- triadCensusDf |>
  group_by(triad) |>
  summarise(
    quant05 = quantile(nnodes, prob = 0.05),
    quant95 = quantile(nnodes, prob = 0.95)
  ) |>
  pivot_longer(
    starts_with("quant"),
    names_to = "quant", names_pattern = "quant(.+)",
    values_to = "nnodes"
  )

plot.new()
# The following code produces the violin plots
ggplot(triadCensusDf, aes(fct_inorder(triad), nnodes)) +
  geom_violin(trim = FALSE, scale = "width") +
  stat_summary(fun = mean, geom = "point", size = 2) +
  geom_boxplot(width = 0.2, fill = "gray", outlier.shape = 4) +
  geom_point(data = triadCensusObs, col = "red", size = 2) +
  geom_line(
    data = triadCensusObs, aes(group = 1), col = "red", linewidth = 0.5
  ) +
  geom_line(
    data = percTriad, mapping = aes(group = quant),
    col = "gray", linetype = "dashed"
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  xlab("triad type") +
  title("Goodness of Fit of Triad Census Counts")


#As it can be seen from the violin plots, the simulated networks for the most part
#don't have an accurate goodness of fit regarding the triad census. While 5 out of
#16 have their observed value within the 5-95% quantile interval, all of the other
#types show that the observed value in net2 is a very big outlier. 
#This bad goodness of fit is also exacerbated by the p_value of the Mahalanobis
#distance: 0; indeed it can be seen that the real value is 576 while the simulated
#ones rarely exceed 20. Therefore our simulations don't model in a satisfactory way
#the data.
