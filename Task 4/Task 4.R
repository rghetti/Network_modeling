setwd("")

# Loading the packages.
library(RSiena)
library(sna)
library(parallel)

source("siena07ToConverge.R")

# Importing the adjacency matrices of the networks
net1 <- as.matrix(read.csv("Glasgow/f1.csv", header = FALSE))
net2 <- as.matrix(read.csv("Glasgow/f2.csv", header = FALSE))
net3 <- as.matrix(read.csv("Glasgow/f3.csv", header = FALSE))

num_pupils <- 129

# Reading the demographic and alcohol consumption characteristics of the students
attributes <- read.csv("Glasgow/demographic.csv", header = TRUE)
alcohol <- as.matrix(
  read.csv("Glasgow/alcohol.csv", header = TRUE)
)

# Reading the information on the distance between the houses of the pupils

log_distance <- as.matrix( read.csv("Glasgow/logdistance.csv", header = FALSE))

# Creation of a Siena network object
friendship <- sienaDependent(
  array(c(net1, net2, net3), dim = c(num_pupils, num_pupils, 3))
)


gender <- coCovar(attributes$gender)
age <- coCovar(attributes$age)
log_dist <- coDyadCovar(log_distance)

#The alcohol consumption is treated as a dependent variable.
alco_beh <- sienaDependent(alcohol, type = "behavior")

mydata <- sienaDataCreate(friendship, gender, age, log_dist, alco_beh)
mydata

## precondition of the analysis
# Data description
# Stability: Jaccard index
print01Report(mydata, modelname = "Glasgow")

#Jacard index is 0.304 for observations 1 and 2 and is 0.351 for observations 2 and 3. Both these numbers
#are greater than 0.3, so data contains enough information to investigate the evolution of the friendship network.


myeff <- getEffects(mydata)
myeff

effectsDocumentation

#We include  transitive triplets effect to represent the dynamics in local structure
myeff <- includeEffects(myeff, transTrip)

# To test the hypothesis H1, we include the in-degree popularity effect, which reflects tendencies to dispersion in in-degrees of the actors; or,
#tendencies for actors with high in-degrees to attract extra incoming ties ‘because’ of their high current in-degrees

myeff <- includeEffects(myeff, name = "friendship", inPop)

#To test the hypothesis H2, we include the covariate main effect of the logarithm of the distance between the
#houses of the pupils
myeff <- includeEffects(myeff, X, name = "friendship", interaction1 = "log_dist")

#To test the hypothesis H3, we include the indegree effect for behavioral evaluation function.
myeff <- includeEffects(myeff, indeg, name = "alco_beh", interaction1 = "friendship")

#To test the hypothesis H4, we include the average similarity effect, defined by the average of centered similarity scores 
#(of alcohol consumption) between pupil and the other pupils to whom he is tied.
myeff <- includeEffects(myeff, avSim, name = "alco_beh", interaction1 = 'friendship')
myeff


## Model estimation --------------------------------------------------------
# Specifying the parameter of the algorithm
myAlgorithm <- sienaAlgorithmCreate(
  projname = "friends_res",
  nsub = 4, n3 = 10000, seed = 239239
)

# Estimate the model
model1 <- siena07(
  myAlgorithm,
  data = mydata, effects = myeff,
  returnDeps = TRUE,
  useCluster = TRUE, nbrNodes = 4, batch = FALSE
)

model1


#Check the convergence of the model
t.conv <- apply(model1$sf, 2, mean) / apply(model1$sf, 2, sd)
overall <- sqrt(t(apply(model1$sf, 2, mean)) %*% solve(cov(model1$sf)) %*% apply(model1$sf, 2, mean))
#all numbers of t.conv are less than 0.1 by absolute value and the value of overall is less than 0.2. 
#Thus we can conclude that the model converges well.

#We double check the convergence using siena07ToConvergence function.

  model1 <- siena07(myAlgorithm,
  data = mydata, effects = myeff, returnDeps = TRUE, prevAns = model1,
  useCluster = TRUE, nbrNodes = 4
  )
  
siena07ToConvergence(myAlgorithm, dat = mydata, eff = myeff)

#The again see that all t-ratios are less than 0.1 by absolute value and the overall maximum convergence ratio
#is less than 0.2. Therefore, we received another evidence of that the model converges well. 

## Check the goodness of fit

# Indegree distribution
gofCoevId <- sienaGOF(
  model1,
  verbose = FALSE,
  varName = "friendship", IndegreeDistribution
)

# Outdegree distribution
gofCoevOd <- sienaGOF(
  model1,
  verbose = FALSE,
  varName = "friendship", OutdegreeDistribution
)

# Triad census
gofCoevTC <- sienaGOF(
  model1,
  verbose = FALSE,
  varName = "friendship", TriadCensus
)



# Behaviour distribution
gofCoevBeh <- sienaGOF(
  model1,
  verbose = FALSE,
  varName = "alco_beh", BehaviorDistribution,
)


plot(gofCoevId)
plot(gofCoevOd)
plot(gofCoevTC, center = TRUE, scale = TRUE)
plot(gofCoevBeh)

descriptives.sienaGOF(gofCoevBeh)

# If we consider the indegree distribution or behavior distribution as auxiliary statistics,
#the  model fits good (the red line lies in the 95% percent interval).
# But if we consider the outdegree distribution or triad census as auxiliary statistics,
#then the model fits not so good: sometimes red line is outside the 95% percent interval. 
#But overall, it happens quite rarely, so we may assume that the model fits well enough.


siena.table(model1, type = "html", sig = TRUE)

#1.The indegree popularity. This statistics is significant under the level of 0.001. 
#The sign of the estimated value is negative, which contradicts with the hypothesis H1.

#2.The centered covariate main effect of logarithm of the distance between the houses of the pupils.
#This statistics is significant under the level of 0.001. 
#The sign of the estimated value is negative, which supports the hypothesis H2.

#3.The the indegree effect for behavioral evaluation function.
#This statistics is not significant even under the level of 0.1.


#4. The average similarity effect.
#This statistics is significant under the level of 0.01.
#The sign of the estimated value is positive, which supports the hypothesis H4.


##Given the estimated model, do we have evidence for selection processes
#only, influece processes only, both selection and influence processes, or
#neither? Argue for your answer.

#The hypotheses H1 and H2 are examples of selection processes. H1 describes a scenario where 
#students are actively choosing their friends based on a specific attribute (popularity). 
#It highlights the selection mechanism where students are drawn to others with certain social traits.
#H2  also reflects a selection process, as the proximity (living nearby) serves as a criterion or facilitator
#for choosing friendships. This is about selecting friends based on geographic or situational factors.
#Thus we see evidence of selection processes.

#On the other hand, the hypotheses H3 and H4 are examples of influence processes. H3 is not statistically relevant, 
#so we concentrate on H4. It describes an influence process because it focuses on how students' behaviors (alcohol consumption)
#are shaped by their social relationships (friends' behavior). It reflects the adjustment or conformity to peer norms.
#Thus we see evidence of influence processes.


#Could you think of two other hypotheses concerning the dynamics of
#friendship and alcohol consumption dynamics that a researcher can test
#using SAOMs? State these hypotheses and how you would operationalize
#the corresponding effects in the evaluation function.


#The first hypothesis is that male students tend to consume more alcohol.
#To test this we can add the main covariate effect (effFrom) to the objective function.

#The second hypothesis is that students of the same gender tend to be friends more than students of the oppostite genders.
#To test this we can add covariate-related similarity (simX) effect.


 


