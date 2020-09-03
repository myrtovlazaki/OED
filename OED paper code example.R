## Example script to assess the information content of experimental designs
## quantified by the determinant of the parameter variance-covariance matrix 

## The example given below reproduces the results in Figure 6B - exponential decay with high variance in a 2-compartment stochastic model 
##
## Written by: Myrto Vlazaki, David Price and Olivier Restif*
## Contact: mv382@cam.ac.uk
## Date: 2/92020

install.packages("SPEEDI.R")
library(SPEEDI.R)


# Function to calculate the KL divergence between predicted and experimental moments

KL.div.o2p <- 
function(N, obs,pred,subs=1:N){
  # Predicted moments
  mu0 <- pred[subs]
  cov0 <- v.2.m(pred[(N+1):length(pred)], N)[subs,subs]
  # Observed moments
  mu1<- obs[subs]
  cov1 <- v.2.m(obs[(N+1):length(pred)], N)[subs,subs]
  # KL divergence
  KL.div(mu0,cov0,mu1,cov1)
}



KL.div <- 
function(mu0,cov0,mu1,cov1){
  k <- length(mu0)
  # inv.C1 <- solve(cov1)
  # as.numeric(sum(diag(inv.C1 %*% cov0)) + (mu1-mu0) %*% inv.C1 %*% (mu1-mu0) - k + log(det(cov1)/det(cov0)))/2
  as.numeric(sum(diag(apply(cov0, 2, FUN = function(x){return(solve(cov1,x))}))) + (mu1-mu0) %*% solve(cov1, (mu1-mu0)) - k + log(det(cov1)/det(cov0)))/2
}


# Function to calculate divergence between predicted and experimental moments

div.measure  <- 
function(N,obs,pred,subs=1:N, div){
  switch(div,
         KL = KL.div.o2p(N, obs,pred,subs=subs),
         Hell = Hell.dist.o2p(N, obs,pred,subs=subs),
         Maha = Maha.dist.o2p(N, obs,pred,subs=subs),
         Chisq = Chisq.dist.o2p(obs,pred))
}


# Function to facilitate identification and replacement of parameters by name

replace.par 
function(x,rep)
{
  if(length(rep)==0) return(x)
  pos <- sapply(names(rep),function(n) which(names(x)==n))
  replace(x,pos,rep)
}


# Function that defines the structure of a network with N compartments

all.poss <- 
function(N){
  
  # All migrations allowed
  migrations <- matrix(1,N,N)
  
  # (except, of course, migrations from a compartment to itself)
  diag(migrations) <- 0
  
  # Birth + death allowed in all compartments
  kills <- rep(1,N)
  replications <- rep(1,N)
  
  return(list('migrations'=migrations, 'kills'=kills, 'replications'=replications))}




# Function that yields the parameter names relevant to a compartmental model with N compartments

network.params <- 
function(network.structure){
  
  capacitated <- !is.null(network.structure$capacities)
  
  # Split network structure into migrations, kills, replications
  migrations <- network.structure$migrations
  kills <- network.structure$kills
  replications <- network.structure$replications
  
  # Check if network is capacitated, act appropriately
  if (capacitated){
    capacities <- network.structure$capacities
  }
  
  else {capacities <- 0}
  
  # Get number of compartments
  N = length(kills)
  
  # Get total number of parameters
  total.parameters = sum(migrations) + sum(kills) + sum(replications) + sum(capacities)
  
  # Initialise vector of names
  out.names <- vector('character', total.parameters)
  
  # Start position tracker at 1
  current.parameter <- 1
  
  # Add permitted migrations
  # If migration permitted, add its name to the name vector
  # and increment the position tracker
  for (i in 1:N){
    for (j in 1:N){
      if (migrations[i,j] == 1){
        
        out.names[current.parameter] = paste('m', i, '.', j, sep='')
        
        current.parameter <- current.parameter + 1
        
      }
    }
  }
  
  
  # Add permitted kills
  # If kill permitted, add its name to the name vector
  # and increment the position tracker
  for (i in 1:N){
    if (kills[i] == 1){
      
      out.names[current.parameter] = paste('k', i, sep='')
      
      current.parameter = current.parameter + 1
      
    }
  }
  
  
  # Add permitted replications
  # If replication permitted, add its name to the name vector
  # and increment the position tracker
  for (i in 1:N){
    if (replications[i] == 1){
      
      out.names[current.parameter] = paste('r', i, sep='')
      
      current.parameter = current.parameter + 1
      
    }
  }
  
  # Add permitted capacitances
  # If capacitance permitted, add its name to the name vector
  # and increment the position tracker
  if (capacitated){
    for (i in 1:N){
      if (capacities[i] == 1){
        
        out.names[current.parameter] = paste('c', i, sep='')
        
        current.parameter = current.parameter + 1
        
      }
    }
  }
  
  
  # Output named vector of zeros with names equal to names of parameters
  out <- vector('numeric', total.parameters)
  
  names(out) <- out.names
  
  return(out)
  
}


# Function to calculate the moments for a network with all possible inter-compartmental migrations, and intra-organ replication and killing

moment.sol.N.general <- 
function(N,t,par,M0,met="Higham08.b"){
  
  # Retrieves r_i from par
  
  replication <- function(i){
    par.name <- paste('r', i, sep='')
    par.val <- par[par.name]
    
    if (is.na(par.val)) {
      return(0)
    }
    
    else {
      return(par.val)
    }
  }
  
  # Retrieves k_i from par
  
  kill <- function(i){
    par.name <- paste('k', i, sep='')
    par.val <- par[par.name]
    
    if (is.na(par.val)) {return(0)}
    else {return(par.val)}
  }
  
  # Retrieves m_i1,i2 from par
  
  migration <- function(i1, i2){
    par.name <- paste('m', i1, '.', i2, sep='')
    par.val <- par[par.name]
    
    if (is.na(par.val)) {return(0)}
    else {return(par.val)}
    
  }
  
  
  # Calculates total migration out of a compartment
  
  total_migration <- function(i){
    
    out <- 0
    
    for (j in 1:N){
      if (i == j) next
      out <- out + migration(i,j)
    }
    
    return(out)
    
  }
  
  # First moment: M.1(t) = (E[n1],E[n2],...,E[nN]) is solution of M.1'(t) = A * M.1
  A <- matrix (0, N, N)
  
  for (i in 1:N){
    A[i,i] = replication(i) -  kill(i) - total_migration(i)
  }
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      A[i,j] = migration(j,i)
      A[j,i] = migration(i,j)
    }
  }
  
  # Solution for first moment
  M.1 <- expm::expAtv(A,M0[1:N],t)$eAtv
  
  
  # Second moment: M.2(t) = (Var(n1),Var(n2),..,Var(nN),Cov(n1,n2),..,Cov(nN-1,nN)) is solution of M.2' = B*M.1 + C*M.2
  # B and C matrices
  B <- matrix(0, N*(N+1)/2, N)
  C <- matrix(0, N*(N+1)/2, N*(N+1)/2)
  
  # Calculates index for (n_i1.n_i2)
  
  index <- function(i1, i2){
    return(i1*N - choose(i1,2) + i2 - i1)
  }
  
  # Var(ni)
  for (i in 1:N){
    
    B[i,i] <- total_migration(i) + kill(i) + replication(i)
    
    C[i,i] <- 2*A[i,i]
  }
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      
      m_ij <- migration(i,j)
      m_ji <- migration(j,i)
      
      B[i,j] <- m_ji
      B[j,i] <- m_ij
      
      ij = index(i,j)
      
      C[i, ij] <- 2*m_ji
      C[j, ij] <- 2*m_ij
      
    }
  }
  
  # Cov(ni, nj)
  
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      
      ij <- index(i,j)
      
      m_ij <- migration(i,j)
      m_ji <- migration(j,i)
      
      B[ij, i] <- -m_ij
      B[ij, j] <- -m_ji
      C[ij, i] <- m_ij
      C[ij, j] <- m_ji
      
      C[ij, ij] <- (A[i,i] + A[j,j])
      
      for (k in 1:N){
        if (k != i && k != j){
          
          # mins and maxes used here as the index function requires i1 < i2
          ik = index(min(i,k),max(i,k))
          jk = index(min(j,k),max(j,k))
          
          C[ij, ik] <- migration(k,j)
          C[ij, jk] <- migration(k,i)
        }
      }
    }
  }
  
  # Matrix exponential solution to integral
  
  D <- matrix(0, N*(N+1)/2, N)
  K <- matrix(0, N*(N+3)/2, N*(N+3)/2)
  
  # As a block matrix, K is:
  
  # C B
  # 0 A
  
  K[1:(N*(N+1)/2),1:(N*(N+1)/2)] <- C[1:(N*(N+1)/2),1:(N*(N+1)/2)]
  K[1:(N*(N+1)/2),(N*(N+1)/2+1):(N*(N+3)/2)] <- B[1:(N*(N+1)/2),1:N]
  K[(N*(N+1)/2+1):(N*(N+3)/2),(N*(N+1)/2+1):(N*(N+3)/2)] <- A[1:N,1:N]
  
  ## We take the exponential of K*t
  
  M <- expm::expm(t*K,m=met)
  
  ## D represents the matrix solution of the second moment equation
  
  D[1:(N*(N+1)/2),1:N] <- M[1:(N*(N+1)/2),(N*(N+1)/2+1):(N*(N+3)/2)]
  
  ## Solution for second moment
  
  M.2 <- D %*% M0[1:N]
  M.2 <- as.numeric(expm::expm(t*C,m = met) %*% M0[(N+1):(N*(N+3)/2)] + M.2)
  
  # Solutions
  
  return(c(M.1,M.2))
  
}





# Modified mde.est from SPEEDI.R package to include all possible inter-compartmental migrations

general.mde.est.noper <- 
  function(N, t, obs.moments, network, divergence, par.start, 
           par.fixed=NULL, init.moments.fun=NULL, init.moments.vals=NULL, 
           subs=rep(list(1:N), times=length(t)), constraint.fun = function(x) {
             return(min(x) < 0)
           },optim.routine=minqa::bobyqa, 
           combine="sum", capacitated=FALSE, ...){
    
    # Check if there are more parameters to be estimated than moments from which to estimate them
    if (length(obs.moments) < length(par.start)){
      print('Warning: There are fewer observed moments than parameters to estimate. Estimates may be very inaccurate.')
    }
    
    # Check if initial conditions are to be estimated
    
    calc.i0 <- FALSE
    if (!is.null(init.moments.vals)){
      calc.i0 <- TRUE
      
      # If initial conditions are to be estimated, but no names included,
      # add names to initial values vector
      
      if (is.null(names(init.moments.vals))){
        names(init.moments.vals) <- paste("init.val.",1:length(init.moments.vals), sep="")
      }
    }
    
    # Initialise parameter vector
    all.par <- network.params(network)
    
    # Update fixed parameters to their value
    all.par <- replace.par(all.par, par.fixed)
    
    # Stop if par.start is not named at all
    if(is.null(names(par.start))){
      stop("Please name entries in par.start")
    }
    # Stop if names of par.start are incorrect
    if (!all(names(par.start) %in% names(all.par))){
      stop('Ensure that names of par.start are correct.')
    }
    
    # Make sure obs.moments input is a matrix
    if (length(t)==1){
      obs.moments <- matrix(obs.moments, nrow=1)
    }
    
    # Pre-define matrix to contain moments of predicted system
    mom.i <- matrix(NA, nrow=length(t), ncol= N * (N + 3)/2)
    div.i <- vector("numeric", length=length(t))
    
    
    # Initialise parameter vector (initial values + dynamic parameters)
    par.start <- c(init.moments.vals, par.start)
    
    # Calculate minimum divergence estimate
    mde.out <- do.call(optim.routine, list(par=par.start, fn = function(par.ext){
      
      names(par.ext) <- names(par.start)
      
      # If any parameters are guessed to be negative, return a prohibitively high value
      if(any(par.ext < 0)) return(1E100)
      
      if (calc.i0){
        # Get dynamic parameters (i0 only needed for initial conditions)
        par <- par.ext[-(1:length(init.moments.vals))]
        # Estimate proportion of inoculum
        i0 <- par.ext[1:length(init.moments.vals)]
      } else{
        par <- par.ext
        i0 <- NULL
      }
      # Update parameters being estimated
      par.i <- replace.par(all.par,par)
      
      # Calculate initial moments from initial conditions (e.g. Poisson from inoculum dose)
      init.mom.i <- init.moments.fun(i0)
      
      # Calculate moments + divergence measure for given parameters
      for (i in 1:length(t)){
        
        mom.i[i, ] <- do.call(moment.sol.N.general, list(N, t[i], 
                                                         par.i, init.mom.i))
        div.i[i] <- div.measure(N = N, obs = obs.moments[i, 
                                                         ], pred = mom.i[i, ], subs = subs[[i]], div = divergence)
      }
      # Combine divergence measures according to given function
      do.call(combine, as.list(div.i))
    }, ...))
    out.pars <- mde.out$par
    par.estimates <- c(out.pars[names(out.pars) %in% names(init.moments.vals)], 
                       replace.par(all.par, out.pars[names(out.pars) %in% names(all.par)]))
    obs.div <- mde.out$value
    return(c(par.estimates, obs.div = obs.div))
  }





######################################################################################
######################################################################################
######################################################################################

# Select 10 parameter sets
set.seed(1)
k <- runif(10,2,2.5)
set.seed(2)

r<- runif(10,2,2.5)

parameters.2a <- matrix(nrow=10, ncol=6)
for (i in 1:10)
{
  set.seed(i)
  parameters.2a[i,] <- c("m1.2"=runif(1,0,0.2), "m2.1"=runif(1,0,0.2),
                         
                         "k1"=k[i], "k2"=r[i], 
                         "r1"=k[i]-0.2, "r2"=r[i]-0.2) }



init <- c(1000,10000,100000)


# Produce simulated data from the model
s2a.small.bspar<- list()
initcondi <- list()
o <- list()

for (i in 1:10)
{
  transitios <-  list(
    c(B=-1, L=1), #m12
    c(B=1, L=-1), #m21
    
    c(B=-1, L=0), #k1
    c(B=0, L=-1), #k2
    
    c(B=1, L=0), #r1
    c(B=0, L=1)) #r2
  
  
  
  
  my.rates <- function(x, params, t) {
    return(c(params$m12*x["B"],
             params$m21*x["L"],
             
             params$k1*x["B"],
             params$k2*x["L"],
             
             params$r1*x["B"],
             params$r2*x["L"])
    )
  }
  params = list(
    m12=parameters.2a[i,1],
    m21=0,
    k1=parameters.2a[i,3],
    k2=parameters.2a[i,4],
    r1=parameters.2a[i,5],
    r2=parameters.2a[i,6])
  
  for (n in 1:3)
  {
    for (j in 1:200)
    {
      
      start.conditions <- round(c(rpois(1,init[n]), 1,1))
      names(start.conditions) <- c("B", "L", "S")
      
      o[[j]] <- ssa.adaptivetau(init.values=start.conditions, transitions=transitios,
                                rateFunc=my.rates, params=params, tf=24)
    }
    initcondi[[n]] <- o
  }
  s2a.small.bspar[[i]] <- initcondi
}

# For the sampling times (=times.of.interest3), determine the number of bacteria in the organs
xfg <- list()
times.of.interest3 <- c(1,2, 12,13, 23,24)
times.in.one.param3 <- matrix(ncol=4, nrow=200)
times.in.20.params3 <- list()
all.times.s2a.small <- list()
line <- vector(length=10000)

for (m in 1:length(times.of.interest3))
{
  for (i in 1:10)
  {
    for (n in 1:3)
    {
      for (j in 1:200)
      {
        line [j] <- 
          which(abs(s2a.small.bspar[[i]][[n]][[j]][,1]-times.of.interest3[m])==
                  min(abs(s2a.small.bspar[[i]][[n]][[j]][,1]-times.of.interest3[m])))
        times.in.one.param3[j,] <- s2a.small.bspar[[i]][[n]][[j]][line [j],]
      }
      xfg[[n]] <- times.in.one.param3
    }
    times.in.20.params3 [[i]] <-xfg
  }
  all.times.s2a.small[[m]] <- times.in.20.params3
}

sample.times.s2a.small <- vector(length=6)

sample.times.s2a.small[1] <- c(1)
sample.times.s2a.small[2] <- c(2)
sample.times.s2a.small[3] <- c(12)
sample.times.s2a.small[4] <- c(13)
sample.times.s2a.small[5] <- c(23)
sample.times.s2a.small[6] <- c(24)


indices3 <- vector(length=6)
indices3 <- sample.times.s2a.small


# Inference eon bootstrapped samples
times <-c(1,2, 12,13, 23,24)
par.list <- list()
par.init.list <- list()
par.init.time.s2a.small.list<- list()
y <- matrix(nrow=200,ncol=7)
x <- matrix(nrow=1, ncol=7)
start <- vector()
finish <- vector()
N <- 2
network <- all.poss(2)
for (i in 1:6)
{
  for (n in 1:3) 
  {
    for (j in 1:10)
    {
      for (b in 1:200)
      {
        
        tryCatch({ 
          
          
          t <- sample.times.s2a.small[i]
          
          
          obs.moments <- rbind(SPEEDI.R:::network.moments(
            cbind(sample(all.times.s2a.small[[i]][[j]][[n]][,2],500, replace=T),
                  sample(all.times.s2a.small[[i]][[j]][[n]][,3],500, replace=T))))
          
          divergence <- "KL"
          
          for (k in 1:1)
          {
            tryCatch({
              par.start <- c("m1.2"=parameters.2a[j,1],
                             
                             "k1"=parameters.2a[j,3],
                             "k2"=parameters.2a[j,4],
                             
                             "r1"=parameters.2a[j,5],
                             "r2"=parameters.2a[j,6])
              
              par.fixed <- c("m2.1"=0)
              init.moments.fun <- function(x){return(
                
                SPEEDI.R:::network.moments(cbind(rpois(500,init[n]),
                                                 rep(1,500),rep(1,500))))}
              subs <- list(c(1,2))
              x[k,] <- general.mde.est.noper(N=N, t=t, obs.moments=obs.moments, network=network, par.start=par.start,
                                             init.moments.fun = init.moments.fun,divergence=divergence, subs=subs,
                                             par.fixed=par.fixed,optim.routine = "optim")
            },
            error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
          }
          y[b,] <- x[which.min(x[,7]),]
        },
        error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
      par.list[[j]]  <- y
      print(y)
    }
    par.init.list [[n]] <- par.list
  }
  par.init.time.s2a.small.list[[i]] <- par.init.list
  
}

# Evaluate the determinant of the parameter variance-covariance matrix
det.init.rep.time.s2a.small.means <- matrix(ncol=6, nrow=3)
det <- vector(length=6)
det.par.rep.s2a.small <- matrix (nrow=10,ncol=3)
det.init.rep.time.s2a.small <- list()

for (u in 1:6)
{
  for (n in 1:3)
  {
    for (j in 1:10)
    {
      det[j] <- 
        det(cov(par.init.time.s2a.small.list[[u]][[n]][[j]][,c(1,3:6)]))
    }
    det.par.rep.s2a.small[,n] <- det
  }
  det.init.rep.time.s2a.small [[u]] <- 
    det.par.rep.s2a.small
  
  det.init.rep.time.s2a.small.means[,u] <- colMeans(det.par.rep.s2a.small)
}
