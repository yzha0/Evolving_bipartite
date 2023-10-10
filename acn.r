#title: "ACN4"
#author: "Eric Zhao and Carrie Diaz Eaton"
#date: "2023-10-09"
# This is an R version, translated from MATLAB

## Parameters set-up

# Clearing the workspace (equivalent to MATLAB's "clear all")
rm(list=ls())

# Setting up parameters 
endT <- 400  # stop time T - internally T is TRUE in R, need to use alternate
pX <- 0.01  # probability of speciation for X species (e.g., plant)
pY <- 0.01  # probability of speciation for Y species (e.g., animal)
maxNX <- 500  # maximum number of X species for preallocation
maxNY <- 1500  # maximum number of Y species for preallocation
v <- 1  # tolerance parameter for trait matching
theta <- 0.1  # assuming all theta values are the same in the model
thetax <- theta #phenotypic variance for x,
thetay <- theta #phenotypic variance for y,
xix <- 0.1  # standard deviation for drift (X species)
xiy <- 0.1  # standard deviation for drift (Y species)
xis <- 0.5  # standard deviation for trait value change during speciation
c <- 0.99  # proportion of connection probability based on historical association

# Initializations

# Starting with a single mutualism

# Trait values
# preallocate matrices
# karray prevents some funky subsetting by R down the road 
x<-karray(0, c(2, maxNX))
xall <- karray(0, c(endT, maxNX))
y<-karray(0, c(2, maxNY))
yall <- karray(0, c(endT, maxNY))

# Initial trait values
x[1,1] <- 0
y[1,1] <- 0

# Interaction matrix (1 if species i and j interact, 0 otherwise)
A <- karray(0, c(1, 1))
# We start with one species each that are connected to each other
A[1,1] <- 1

NsppX = 1; # total # of species of type X that have ever existed
iXextant = 1;  #the indicies of the extant species (out of the NsppX in x)

NsppY = 1; # total # of species of type X that have ever existed
iYextant = 1;  #the indicies of the extant species (out of the NsppX in x)

# List to keep track of extinct species
#total # of extant species 
NextantsppX <- length(iXextant)
NextantsppY <- length(iYextant)
#To keep track of phylogeny, we use an indexing pointer-like system  The
#first row holds the tree topology and the second row holds the metric of
#the tree (what generation speciation occured)
# For example, [1 3 2; 0 150 100] means
#       1     at start       
#      / \
#     /   \
#    1     2   at generation 100
#   / \     \
#  1   3     2 at generation 150

#phylogenyX<-matrix(0, nrow = 2, ncol = NsppX)
#phylogenyX<-rbind(c(1, 0), c(2, 0))
#phylogenyY<-matrix(0, nrow = 2, ncol = NsppY)
#phylogenyY<-rbind(c(1, 0), c(2, 0))


connectedness <- matrix(0, 1, endT)
connectX <- matrix(0, 1, endT)
connectY <- matrix(0, 1, endT)

# creating a list storing interaction matrix at each time stamp
ts_list<-list()

## Simulation Loop

for (n in 1:endT) {
  # Drift in traits for existing species
  xall[n,] <- x[1, ] 
  yall[n,] <- y[1, ] 
  
  #average trait value of the interactors with any
  #particular species
  Yextant <- karray(y[1, iYextant], c(1,NextantsppY)) #trait values of extant only
  AvgXconnection<-A %*% t(Yextant)/rowSums(A)
  
  Xextant <- karray(x[1, iXextant], c(NextantsppX,1))
  AvgYconnection<-t(A) %*% Xextant/rowSums(t(A))

  Aold <- A
  
  muI<-AvgXconnection
  muJ<-AvgYconnection
  
  w<-A
  for (i in 1:NextantsppX) {nI <- sum(A[i,])} #could prob use an apply
  for (j in 1:NextantsppY) {mJ <- sum(A[,j])}
  
  numerator1 <- (theta + v) / (mJ * sqrt(theta^2 + theta * v + v^2))
  numerator2 <- (theta + v) / (nI * sqrt(theta^2 + theta * v + v^2))
      
  expTerm1 <- -1 * (v * (x[1, iXextant[i]] - muJ[j]) * (x[1, iXextant[i]] + muJ[j] - 2 * muI[i]) + theta * (muI[i] - muJ[j]) * (muI[i] + muJ[j] - 2 * x[1, iXextant[i]]))
  expTerm2 <- -1 * (v * (y[1, iYextant[j]] - muI[i]) * (y[1, iYextant[j]] + muI[i] - 2 * muJ[j]) + theta * (muJ[j] - muI[i]) * (muI[i] + muJ[j] - 2 * y[1, iYextant[j]]))
      
  denominator <- theta^2 + theta * v + v^2
      
  #w[i,j] <- numerator1 * exp(expTerm1 / denominator)
  w[i,j] <- numerator2 * exp(expTerm2 / denominator)
   
  m<-0
  NoldsppX <- NsppX
  
  # only look at the extant species
  # i is the index of each of the extant of individuals/species
  # recall we need to put NaN where spec have gone extinct, not
  # delete them completely, so this limits us to only calling extant
  # species to speciate and evolve
  for (i in iXextant) {
    m <- m + 1  # counts which species we are on to use as an index for iXextant, A
    # where i acts as an index for x, i.e iXextant[m] = i
    
    # variation from one generation to the next
    nI <- sum(Aold[m,])
    
    thetaI <- thetay / nI  # CLT
    
    alpha <- (thetaI^2 + thetaI*v + v^2) / (thetaI^2 + thetaI*v + v^2 + thetax*v)
    gamma1 <- thetaI * thetax / (thetaI^2 + thetaI*v + v^2 + thetax*v)
    
    non_zero_indices <- which(Aold[m,] != 0)
    x[2, i] <- alpha * x[1, i] + (1 - alpha) * AvgXconnection[m] + (AvgXconnection[m] - mean(AvgYconnection[non_zero_indices])) * gamma1 + rnorm(1) * xix
    
    
    
    # speciation event - which occurs with probability pX
    if (runif(1) < pX) {
      NsppX <- NsppX + 1  # keeps track of how many spp there are now
      
      if (NsppX > maxNX) {
        cat("error, the number of species type X is greater than the number allocated\n")
        toobig <- toobig + 1
        break
      }
      
      # add new spp trait
      x[2, NsppX] <- x[1, i] + rnorm(1) * xis
      
      #add a new row for the new X species
      A <- rbind(A, A[m, ])  # inherits parental connections
      w <- rbind(w, w[m, ])  # and those connections inherit the same fitness value, at least initially
      
      
      # now we make sure that we know that this new species is a sister species to species i
      #boo <- which(phylogenyX[1, ] == i)
      #if((NsppX-1)>=(boo[1]+1)){
      #  for (j in seq((NsppX-1),(boo[1]+1),-1)) {
      #    phylogenyX <- cbind(phylogenyX,phylogenyX[, j])
      #  }
      #}
      #if((boo[1]+1)<=dim(phylogenyX)[2]){
      #  phylogenyX[, boo[1]+1] <- c(NsppX, n)
      #} else {phylogenyX<-cbind(phylogenyX, c(NsppX, n))}
    }
    
    
  }   
  
  m <- 0
  NoldsppY <- NsppY
  
  for (i in iYextant) {
    
    m <- m + 1; #counts which species we are on
    mJ<- sum(Aold[,m])

    thetaJ <- thetax / mJ  # CLT
    
    beta <- (thetaJ^2 + thetaJ*v + v^2) / (thetaJ^2 + thetaJ*v + v^2 + thetay*v)
    gamma2 <- thetaJ * thetax / (thetaJ^2 + thetaJ*v + v^2 + thetay*v)
    
    non_zero_indices <- which(Aold[,m] != 0)
    y[2, i] <- beta * y[1, i] + (1 - beta) * AvgYconnection[m] + (AvgYconnection[m] - mean(AvgXconnection[non_zero_indices])) * gamma2 + rnorm(1) * xix
    
    
    # speciation event - which occurs with probability pY
    if (runif(1) < pY) {
      NsppY <- NsppY + 1  # keeps track of how many spp there are now
      
      if (NsppY > maxNY) {
        cat("error, the number of species type Y is greater than the number allocated\n")
        toobig <- toobig + 1
        break
      }
      
      # add new spp trait
      y[2, NsppY] <- y[1, i] + rnorm(1) * xis
      
      #add a column for the new y spp
      A <- cbind(A, A[, m])  # inherits parental connections
      w <- cbind(w, w[, m])  # and those connections inherit the same fitness value, at least initially
      
      # now we make sure that we know that this new species is a sister species to species i
     # boo <- which(phylogenyY[1, ] == i)
      #if((NsppY-1)>=boo[1]){
      #  for (j in (NsppY-1):boo[1]) {
      #    phylogenyY <- cbind(phylogenyY,phylogenyY[, j])
      #  }
      #}
      #if((boo[1]+1)<=dim(phylogenyY)[2]){
      #  phylogenyY[, boo[1]+1] <- c(NsppY, n)
      #} else {phylogenyY<-cbind(phylogenyY, c(NsppY, n))}
    }
    
  }
  
  #interaction matrix - these are in the order in which spp appeared, but not including extinct species 
  count<-1
  im<-0
  jm<-0
  
  if((NoldsppX+1)<=NsppX){
    iXextantnew <- c(iXextant,(NoldsppX+1):NsppX)#give the lineages of X that are extant
    } else {
    iXextantnew<-iXextant
    }
  
  
  if((NoldsppY+1)<=NsppY){
    iYextantnew <- c(iYextant,(NoldsppY+1):NsppY) #give the linages of Y that are extant
  } else {
      iYextantnew<-iYextant
      }
  Anew <- karray(0, c(length(iXextantnew),length(iYextantnew)))
  
  for (i in iXextantnew) {
    im<-im+1
    for (j in iYextantnew) {
      jm<-jm+1
      matching <- c*A[im,jm]+ (1-c)*w[im,jm]
      
      if(matching>runif(1)){
        Anew[im, jm]<-1
      }
    }
    jm<-0
  }
  
  # Identify extant species with connections
  XextantA <- which(rowSums(Anew) != 0)
  Xextant <- iXextantnew[XextantA]
  NextantsppX <- length(iXextant)
  
  
  
  YextantA <- which(colSums(Anew) != 0)
  iYextant <- iYextantnew[YextantA]
  NextantsppY <- length(iYextant)
  
  # Check for extinction of species
  if (NsppX == 0 && NsppY == 0) {
    cat("extinction of both species\n")
    extinctions <- extinctions + 1
    break # replaces MATLAB's break in this context
    # note - next will skip an interation
    # we want break because it will stop the whole loop
  } else if (NsppX == 0) {
    cat("extinction of species X\n")
    extinctions <- extinctions + 1
    break
  } else if (NsppY == 0) {
    cat("extinction of species Y\n")
    extinctions <- extinctions + 1
    break
  }
  
  # Update interaction matrix
  A <- karray(Anew[XextantA, YextantA], c(NextantsppX, NextantsppY))
  
  # Update traits for next generation
  x[1, ] <- x[2, ]
  x[2, ] <- rep(NaN, maxNX)
  y[1, ] <- y[2, ]
  y[2, ] <- rep(NaN, maxNY)
  
  # Calculate connectedness
  connectedness[1, n] <- sum(sum(A)) / (NextantsppX * NextantsppY)
  # connectX[1, n] <- mean(rowSums(A)) / NextantsppY
  # connectY[1, n] <- mean(colSums(A)) / NextantsppX
}
  
  
## Visualization

# Plot the evolution of species traits over time for X and Y
par(mfrow=c(3,2))  # Set up a 3x2 plotting grid

plot(1:endT, xall[1:endT,1], type="p", cex=0.5, ylab="Mean Trait x", main="Evolution of Species Trait X")
plot(1:endT, yall[1:endT,1], type="p", cex=0.5, ylab="Mean Trait y", main="Evolution of Species Trait Y")


# Display the adjacency matrix of interactions(this is t=400)
x.at<-seq(0, NsppX+1)
y.at<-seq(0, NsppY+1)
image(t(A[1:NsppX, 1:NsppY]), axes=FALSE, main="Adjacency Matrix of Interactions", xlab="Pollinator species", ylab="Plant species")
axis(1, at=x.at)
axis(2, at=y.at)


# Display the distribution of connections per pollinator
hist(rowSums(A[1:NsppX, 1:NsppY]), main="Distribution of Connections per Pollinator", xlab="Number of connections", ylab="Frequency") 

# Display connectedness over time
plot(1:endT, connectedness, type="l", ylab="Connectedness", xlab="Iterations", main="Connectedness Over Time")
