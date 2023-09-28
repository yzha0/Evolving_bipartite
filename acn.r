#This is the R version of the ACN4.Rmd
# Original Code authored in Matlab by CDE
# The Rmd is an adaption to R by YZ
# Date - 9/28
# Clearing the workspace (equivalent to MATLAB's "clear all")
rm(list=ls())

# Setting up parameters (translated from MATLAB to R)
T <- 400  # stop time T
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
x<-matrix(0, nrow=2, ncol=maxNX)
xall <- matrix(0, nrow=T, ncol=maxNX)
y<-matrix(0, nrow = 2, ncol = maxNY)
yall <- matrix(0, nrow=T, ncol=maxNY)

# Initial trait values
x[1,1] <- 0
y[1,1] <- 0

# Interaction matrix (1 if species i and j interact, 0 otherwise)
A <- matrix(0, nrow=1, ncol=1)
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

#keep track of phylogeny using indexing pointer--The
#first row holds the tree topology and the second row holds the metric of
#the tree (what generation speciation occured)
phylogenyX<-matrix(0, nrow = 2, ncol = NsppX)
phylogenyX<-rbind(c(1, 0), c(2, 0))
phylogenyY<-matrix(0, nrow = 2, ncol = NsppY)
phylogenyY<-rbind(c(1, 0), c(2, 0))


connectedness <- matrix(0, 1, T)
connectX <- matrix(0, 1, T)
connectY <- matrix(0, 1, T)
```

## Simulation Loop
```{r}
# Main simulation loop
for (n in 1:T) {
  # Drift in traits for existing species
  xall[n,] <- x[1, ] 
  yall[n,] <- y[1, ] 
  
  #average trait value of the interactors with any
  #particular species
  if(is.matrix(A)){
    AvgXconnection<-A %*% as.matrix(y[1, iYextant])/as.matrix(rowSums(A))
    AvgYconnection<-t(A) %*% as.matrix(x[1, iXextant])/as.matrix(colSums(A))
  }
  
  else{
    AvgXconnection <- A * y[1, iYextant]
    AvgYconnection <- A * x[1, iXextant]
  }
  Aold <- A
  
  muI<-AvgXconnection
  muJ<-AvgYconnection
  
  w<-A
  for (i in 1:NextantsppX) {
    if(is.matrix(A)){
      nI <- sum(A[i,])
    }
    else{nI <-sum(A)}
    for (j in 1:NextantsppY) {
      if(is.matrix(A)){
        mJ <- sum(A[,j])
      }
      else{mJ<- sum(A)}
      
      numerator1 <- (theta + v) / (mJ * sqrt(theta^2 + theta * v + v^2))
      numerator2 <- (theta + v) / (nI * sqrt(theta^2 + theta * v + v^2))
      
      expTerm1 <- -1 * (v * (x[1, iXextant[i]] - muJ[j]) * (x[1, iXextant[i]] + muJ[j] - 2 * muI[i]) + theta * (muI[i] - muJ[j]) * (muI[i] + muJ[j] - 2 * x[1, iXextant[i]]))
      expTerm2 <- -1 * (v * (y[1, iYextant[j]] - muI[i]) * (y[1, iYextant[j]] + muI[i] - 2 * muJ[j]) + theta * (muJ[j] - muI[i]) * (muI[i] + muJ[j] - 2 * y[1, iYextant[j]]))
      
      denominator <- theta^2 + theta * v + v^2
      
      #w[i,j] <- numerator1 * exp(expTerm1 / denominator)
      if(is.matrix(w)){
        w[i,j] <- numerator2 * exp(expTerm2 / denominator)
      }
      else{w<-numerator2 * exp(expTerm2 / denominator)}
    }
  }
  
  
  
  m<-0
  NoldsppX <- NsppX
  
  for (i in iXextant) {  # only look at the extant species
    # i is the index of each of the extant of individuals/species
    # recall we need to put NaN where spec have gone extinct, not
    # delete them completely, so this limits us to only calling extant
    # species to speciate and evolve
    
    m <- m + 1  # counts which species we are on to use as an index for iXextant, A
    # where i acts as an index for x, i.e iXextant[m] = i
    
    # variation from one generation to the next
    if(is.matrix(Aold)){
      nI <- sum(Aold[m,])
    }
    else{nI<-sum(Aold)}
    thetaI <- thetay / nI  # CLT
    
    alpha <- (thetaI^2 + thetaI*v + v^2) / (thetaI^2 + thetaI*v + v^2 + thetax*v)
    gamma1 <- thetaI * thetax / (thetaI^2 + thetaI*v + v^2 + thetax*v)
    
    if(is.matrix(Aold)){
      non_zero_indices <- which(Aold[m,] != 0)
      x[2, i] <- alpha * x[1, i] + (1 - alpha) * AvgXconnection[m] + (AvgXconnection[m] - mean(AvgYconnection[non_zero_indices])) * gamma1 + rnorm(1) * xix
    }
    else{
      x[2, i] <- alpha * x[1, i] + (1 - alpha) * AvgXconnection[m] + (AvgXconnection[m] - mean(AvgYconnection)) * gamma1 + rnorm(1) * xix
    }
    
    
    
    
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
      if(is.matrix(A)){
        A <- rbind(A, A[m, ])  # inherits parental connections
        w <- rbind(w, w[m, ])  # and those connections inherit the same fitness value, at least initially
      }
      # else{
      #    A<-rbind(A, A[m])
      #    w<-rbind(w, w[m])
      #  }
      # now we make sure that we know that this new species is a sister species to species i
      boo <- which(phylogenyX[1, ] == i)
      if((NsppX-1)>=(boo[1]+1)){
        for (j in (NsppX-1):(boo[1]+1)) {
          phylogenyX <- cbind(phylogenyX,phylogenyX[, j])
        }
      }
      if((boo[1]+1)<=dim(phylogenyX)[2]){
        phylogenyX[, boo[1]+1] <- c(NsppX, n)
      }
      else{phylogenyX<-cbind(phylogenyX, c(NsppX, n))}
    }
    
    
  }   
  
  m <- 0
  NoldsppY <- NsppY
  
  for (i in iYextant) {
    
    m <- m + 1; #counts which species we are on
    if(is.matrix(Aold)){
      mJ<- sum(Aold[,m])
    }
    else{mJ<-sum(Aold)}
    thetaJ <- thetax / mJ  # CLT
    
    beta <- (thetaJ^2 + thetaJ*v + v^2) / (thetaJ^2 + thetaJ*v + v^2 + thetay*v)
    gamma2 <- thetaJ * thetax / (thetaJ^2 + thetaJ*v + v^2 + thetay*v)
    
    if(is.matrix(Aold)){
      non_zero_indices <- which(Aold[,m] != 0)
      y[2, i] <- beta * y[1, i] + (1 - beta) * AvgYconnection[m] + (AvgYconnection[m] - mean(AvgXconnection[non_zero_indices])) * gamma2 + rnorm(1) * xix
    }
    else{
      y[2, i] <- beta * y[1, i] + (1 - beta) * AvgYconnection[m] + (AvgYconnection[m] - mean(AvgXconnection)) * gamma2 + rnorm(1) * xix
    }
    
    
    
    # speciation event - which occurs with probability pY
    if (runif(1) < pY) {
      NsppY <- NsppY + 1  # keeps track of how many spp there are now
      
      if (NsppY > maxNY) {
        cat("error, the number of species type X is greater than the number allocated\n")
        toobig <- toobig + 1
        break
      }
      
      # add new spp trait
      y[2, NsppY] <- y[1, i] + rnorm(1) * xis
      if(is.matrix(A)){
        A <- cbind(A, A[, m])  # inherits parental connections
        w <- cbind(w, w[, m])  # and those connections inherit the same fitness value, at least initially
      }
      # else{
      #    A<-rbind(A, A[m])
      #    w<-rbind(w, w[m])
      #  }
      # now we make sure that we know that this new species is a sister species to species i
      boo <- which(phylogenyY[1, ] == i)
      if((NsppY-1)>=boo[1]){
        for (j in (NsppY-1):boo[1]) {
          phylogenyY <- cbind(phylogenyY,phylogenyY[, j])
        }
      }
      if((boo[1]+1)<=dim(phylogenyY)[2]){
        phylogenyY[, boo[1]+1] <- c(NsppY, n)
      }
      else{phylogenyY<-cbind(phylogenyY, c(NsppY, n))}
    }
    
  }
  
  #interaction matrix - these are in the order in which spp appeared, but not including extinct species 
  count<-1
  im<-0
  jm<-0
  
  if((NoldsppX+1)<=NsppX){
    iXextantnew <- c(iXextant,(NoldsppX+1):NsppX)#give the linages of X that are extant
  } 
  else{iXextantnew<-iXextant}
  if((NoldsppY+1)<=NsppY){
    iYextantnew <- c(iYextant,(NoldsppY+1):NsppY) #give the linages of Y that are extant
  }
  else{iYextantnew<-iYextant}
  Anew <- matrix(0, length(iXextantnew),length(iYextantnew))
  
  for (i in iXextantnew) {
    im<-im+1
    for (j in iYextantnew) {
      jm<-jm+1
      if(is.matrix(A)){
        matching <- c*A[im,jm]+ (1-c)*w[im,jm]
      }
      else{
        matching <- c*A+ (1-c)*w
      }
      if(matching>runif(1)){
        Anew[im, jm]<-1
      }
    }
    jm<-0
  }
  
  # Identify extant species with connections
  XextantA <- which(rowSums(Anew) != 0)
  iXextant <- iXextantnew[XextantA]
  NextantsppX <- length(iXextant)
  
  YextantA <- which(colSums(Anew) != 0)
  iYextant <- iYextantnew[YextantA]
  NextantsppY <- length(iYextant)
  
  # Check for extinction of species
  if (NsppX == 0 && NsppY == 0) {
    cat("extinction of both species\n")
    extinctions <- extinctions + 1
    next # replaces MATLAB's break in this context
  } else if (NsppX == 0) {
    cat("extinction of species X\n")
    extinctions <- extinctions + 1
    next
  } else if (NsppY == 0) {
    cat("extinction of species Y\n")
    extinctions <- extinctions + 1
    next
  }
  
  # Update interaction matrix
  A <- Anew
  
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
```


```{r Packages, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)

# Store string containing all required packages
my_packages <- c('bipartite', 'RColorBrewer', 'igraph')

# Store all installed packages
ya_loaded <- (.packages())

# Check whether required packages are already installed and grab only those that still need installation
need_load<-my_packages[!(my_packages %in% ya_loaded)]

# Load required packages
lapply(need_load, require, character.only = TRUE)
```


``` {r}
# assign plants' name to the plants
rname_vec<-rep(0, length(A[,1]))
for (i in 1:length(A[,1])) {
  rname_vec[i]<-paste("PL",i,sep = "")
}
row.names(A)<-rname_vec

# view the new matrix
View(A)

```

## Analyze the bipartite graphs

### First convert adjacency matrix to edgelist
```{r}
#
polli.g<- graph.incidence(A, weighted = TRUE)
polli_el <- get.edgelist(polli.g)
View(polli_el)
```



### Structuring the bipartite graphs
```{r}
# Set the arbitrary char as webID 
webID <- data.frame(matrix("pp", nrow = nrow(polli_el), ncol = 1))
polli_new <- cbind(polli_el, webID)
# assign column names
colnames(polli_new)<- c("plant", "pollinator", "webID")
View(polli_new)
```



### Plotting network interactions between two node types
```{r}
#Generate the network matrix in bipartite
polli_web<- frame2webs(polli_new, varnames = c("plant", "pollinator", "webID"), type.out = "list", emptylist = TRUE)
# create a color vector
cols1<-c( '#8214a0', '#005ac8', '#00a0fa', '#fa78fa', '#14d2dc', '#aa0a3c', '#fa7850', '#0ab45a', '#f0f032', '#a0fa82', '#fae6be') 
# plot the two-dimensional matrix to a bipartite graph
plotweb(polli_web$"pp", method='cca', labsize=1,col.interaction=cols1, bor.col.interaction=cols1)                    
```




### Visualizing the interaction matrix
This code will visualize the interaction/adjacency matrix. Look up the documentation, play with the command, and then comment within the code block to explain how to use visweb.

``` {r}
# set the color of levels of interactions
cols2<-c('white', '#006e82')
#you can also try other palettes
#cols2 <- brewer.pal(2, "Reds")
#visualize the interaction matrix
visweb(polli_web$'pp', type='nested', labsize=0.3, plotsize=28 ,square="def",def.col=cols2)
```



## Visualization
```{r}

# Plot the evolution of species traits over time for X and Y
par(mfrow=c(3,2))  # Set up a 3x2 plotting grid

plot(1:T, xall[1:T,1], type="p", cex=0.5, ylab="Mean Trait x", main="Evolution of Species Trait X")
plot(1:T, yall[1:T,1], type="p", cex=0.5, ylab="Mean Trait y", main="Evolution of Species Trait Y")


# Display the adjacency matrix of interactions
x.at<-seq(0, NsppX+1)
y.at<-seq(0, NsppY+1)
image(t(A[1:NsppX, 1:NsppY]), axes=FALSE, main="Adjacency Matrix of Interactions", xlab="Pollinator species", ylab="Plant species")
axis(1, at=x.at)
axis(2, at=y.at)


# Display the distribution of connections per pollinator
hist(rowSums(A[1:NsppX, 1:NsppY]), main="Distribution of Connections per Pollinator", xlab="Number of connections", ylab="Frequency") 

# Display connectedness over time
plot(1:T, connectedness, type="l", ylab="Connectedness", xlab="Iterations", main="Connectedness Over Time")
