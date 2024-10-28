# Network Modeling - HS 2024
# C. Stadtfeld, A. Uzaheta and I. Smokovic
# Social Networks Lab
# Department of Humanities, Social and Political Sciences
# ETH Zurich
# 14 October 2024
#
# Assignment 1 - Task 2

library('ergm')
 #function calculating the given statics of an adjacency matrix:
 stat=function(net){
   z1=0
   z2=0
   z3=0
   indegree=0
   nvertices <- nrow(net) 
   
   for(k in 1:nvertices){
     for(q in 1:nvertices){
       z1=z1+net[k,q]         #calculating the num of edges
       if(k<q){
         z2=z2+net[k,q]*net[q,k]  #calculating reciprocity
       }
       if(net[q,k]==1){         #calculating indegree of node k
         indegree=indegree+1
       }
     }
     z3=z3+choose(indegree,2)        #calculating 2-istar and resetting for next node
     indegree=0
     
   }
   return(list(z1,z2,z3))
   
 }

# MHstep ------------------------------------------------------------------
#' Simulate the next step of a network in Markov chain using Metropolis-Hasting
#' 
#' The function `MHstep` simulates the the Metropolis-Hastings step that defines
#' the Markov chain whose stationary distribution is the ERGM with 
#' edge, mutual and nodematch statistics
#'
#' @param net an object of class `matrix`. Adjacency matrix of the network.
#' @param theta1 a numeric value. The value of the edge parameter of the ERGM.
#' @param theta2 a numeric value. The value of the mutual parameter of the ERGM.
#' @param theta3 a numeric value. The value of the istar(2) parameter of the ERGM.
#'
#' @return next state of the Markov Chain
#'
#' @examples
#' MHstep(
#'   matrix(c(0, 1, 0, 0, 0, 0, 1, 1, 0), nrow = 3, ncol = 3),
#'   -log(0.5), log(0.4), log(.8)
#' )
MHstep <- function(net, theta1, theta2, theta3){
  
  # Number of vertices in the network
  nvertices <- nrow(net) 
  
  # Choose randomly two vertices, prevent loops {i,i} with replace = FALSE
  tie <- sample(1:nvertices, 2, replace = FALSE) 
  i <- tie[1]
  j <- tie[2]
  
  # Compute the change statistics
  
  net2=net    #creating the matrix with different i->j
  
  if(net[i,j]==0)
    net2[i,j]=1
  else 
    net2[i,j]=0
  #calculating the statistics for the input
  stat1=stat(net)
  
  #computing the statistics for the changed graph
  stat2=stat(net2)
  
  
  # Compute the probability of the next state 
  # according to the Metropolis-Hastings algorithm
  
  # computing both exponential
  stat1=as.double(stat1)
  stat2=as.double(stat2)
  exp1=exp(theta1*stat1[1]+theta2*stat1[2]+theta3*stat1[3])
  exp2=exp(theta1*stat2[1]+theta2*stat2[2]+theta3*stat2[3])
  
  # computing transition probability 
  
  p=min(1,exp2/exp1)
  
  # Select the next state: 

  #sample with probability p the change of state
  outcome=sample(c(0,1),size=1,prob=c(1-p,p))
  if(outcome==1)
    net=net2
  
  # Return the next state of the chain
  return(net)
}

# Markov Chain simulation -------------------------------------------------
#' The function MarkovChain simulates the networks from the ERGM with 
#' edge, mutual and nodematch statistics
#'
#' @param net an object of class `matrix`. Adjacency matrix of the network.
#' @param theta1 a numeric value. The value of the edge parameter of the ERGM.
#' @param theta2 a numeric value. The value of the mutual parameter of the ERGM.
#' @param theta3 a numeric value. The value of the istar(2) parameter of the ERGM.
#' @param burnin an integer value.
#'   Number of steps to reach the stationary distribution.
#' @param thinning an integer value. Number of steps between simulated networks.
#' @param nNet an integer value. Number of simulated networks to return as output.
#'
#' @return a named list:
#'   - netSim: an `array` with the adjancency matrices of the simulated networks.
#'   - statSim: a `matrix` with the value of the statistic defining the ERGM.
#'
#' @examples
#' MarkovChain(
#'   matrix(c(0, 1, 0, 0, 0, 0, 1, 1, 0), nrow = 3, ncol = 3),
#'   -log(0.5), log(0.4), log(.8)
#' )
MarkovChain <- function(
    net,
    theta1, theta2, theta3,
    burnin = 10000, thinning = 1000, nNet = 1000){
  
  # Burnin phase: repeating the steps of the chain "burnin" times  
  nvertices <- nrow(net)
  burninStep <- 1 # counter for the number of burnin steps
  
  # Perform the burnin steps
  for(burninStep in 1:burnin){
    net=MHstep(net, theta1, theta2, theta3)
  }
  
  # After the burnin phase we draw the networks
  # The simulated networks and statistics are stored in the objects
  # netSim and statSim
  netSim <- array(0, dim = c(nvertices, nvertices, nNet))
  statSim <- matrix(0, nNet, 3)
  thinningSteps <- 0 # counter for the number of thinning steps
  netCounter <- 1 # counter for the number of simulated network
  
  while(netCounter<=nNet){
    while(thinningSteps<thinning){
      net=MHstep(net, theta1, theta2, theta3)  #performing 1000 transitions
      thinningSteps=thinningSteps+1
    }
    netSim[1:nvertices,1:nvertices,netCounter]=net #saving the current network and its statistics
    statSim[netCounter,1:3]=as.double(stat(net))
    netCounter=netCounter+1
    thinningSteps=0                             #resetting the counter
  }
  
  # Return the simulated networks and the statistics
  return(list(netSim = netSim, statSim = statSim))
}


t1=-2.76
t2=0.68
t3=0.05
ad.matrix=matrix(0, 21,21)

MC_simulation=MarkovChain(ad.matrix,t1,t2,t3)

guess1=-2.5
guess2=0.8
guess3=0.1
a.matrix=matrix(0,21,21)
MC2.0=MarkovChain(a.matrix,guess1,guess2,guess3)

alpha=-2
beta=0.9
gamma=0.2
b=matrix(0,21,21)
MC3.0=MarkovChain(b,alpha,beta,gamma)

t=-2.2
y=0.85
u=0.14
c=matrix(0,21,21)
MC4.0=MarkovChain(c,t,y,u)

q=-2.13
w=0.87
e=0.17
d=matrix(0,21,21)
MC5.0=MarkovChain(d,q,w,e)
avr.dens=mean(MC5.0$statSim[1:1000,1])
avr.rec=mean(MC5.0$statSim[1:1000,2])
avr.star=mean(MC5.0$statSim[1:1000,3])

z=-2.1225
x=0.8721
v=0.1721
f=matrix(0,21,21)
MC6.0=MarkovChain(f,z,x,v)
avr.dens=mean(MC6.0$statSim[1:1000,1])
avr.rec=mean(MC6.0$statSim[1:1000,2])
avr.star=mean(MC6.0$statSim[1:1000,3])
