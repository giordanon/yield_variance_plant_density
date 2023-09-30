

# Process model for the mean
quad_plateau <- function(b1, AOPD,x){
  y  = ifelse(x > AOPD,
              b1*AOPD+(-b1/(2*AOPD))*AOPD^2,
              # else
              b1*x+(-b1/(2*AOPD))*x^2)
  
  return(y)
}

#Process Model for the variance
exp_plateau <- function(a1, maxVPD,minVPD,x){
  a2 = (-a1/(2*maxVPD))
  y = ifelse(x == 0, 0,
             ifelse(x > minVPD, 
                    exp(a1*minVPD+a2*minVPD^2),
                    # else
                    exp(a1*x+a2*x^2)))
  
  return(y)
}

# Likelihood function
likelihood <- function(parameters, 
                       x,
                       y
){
  b1 = parameters[1]
  AOPD = parameters[2]
  
  a1 = parameters[3]
  maxVPD = parameters[4]
  minVPD = parameters[5]
  
  pred = quad_plateau(b1, AOPD, x) # Process model to the mean
  
  predVar <- exp_plateau(a1,maxVPD,minVPD,x) # Process model for the variance
  
  singlelikelihoods = dnorm(y, 
                            mean = pred, 
                            sd = sqrt(predVar), 
                            log = T)
  sumll = sum(singlelikelihoods)
  return(sumll)
}

# Priors 
priors <- function(parameters){
  
  b1 = parameters[1]
  AOPD = parameters[2]
  
  a1 = parameters[3]
  maxVPD = parameters[4]
  minVPD = parameters[5]
  
  # Specify priors
  bprior = dunif(b1, min = 0.001, max = 0.2, log = T)
  AOPDprior = dunif(AOPD, min = 10, max = 400, log = T)
  
  aprior = dunif(a1, min = 0, max = 0.05, log = T)
  # Upper bound from Canadian data 192
  minVPDprior = dunif(minVPD, min = 10, max = 192, log = T)
  # Upper bound from Canadian data 124
  maxVPDprior = dunif(maxVPD, min = 10, max = 124, log = T)
  
  return(bprior+AOPDprior+aprior+minVPDprior+maxVPDprior)
}

# Likelihood * Prior
posterior <- function(parameters,x,y){
  
  return (likelihood(parameters,x,y) + priors(parameters))
}

# Proposal values
proposalFunction <- function(paramOld){
  # Expected value
  b1.try <- rnorm(1, mean = paramOld[1], sd = sigma.tune[1])
  AOPD.try <- rnorm(1, mean = paramOld[2], sd = sigma.tune[2])
  # Variance
  a1.try <- rnorm(1, mean = paramOld[3], sd = sigma.tune[3])
  maxVPD.try <- rnorm(1, mean = paramOld[4], sd = sigma.tune[4])
  minVPD.try <- rnorm(1, mean = paramOld[5], sd = sigma.tune[5])
  
  return(c(b1.try,AOPD.try, a1.try,maxVPD.try, minVPD.try))
}

# MCMC chain using metropolis hastings 
metropolisMCMC <- function(iterations,# number of iterations
                           burnIn,# burn in interval
                           thinning, # every how many samples?
                           startValue, # starting values for model parameters
                           x, 
                           y
                           ){
  # Chain for E(y)
  chain = array(dim = c(iterations+1,5))
  chain[1,] = startValue
  
  # Loop for E(y) and Var(y) together
  for (i in 1:iterations){
    
    proposal = proposalFunction(chain[i,])
    
    probability = exp(posterior(proposal,x,y) - posterior(chain[i,],x,y))
    if (runif(1) < probability){
      # Accept
      chain[i+1,] = proposal
    }else{
      # Reject
      chain[i+1,] = chain[i,]
    }
    
  }
  # Burn in interval
  burnSamples = iterations * burnIn
  chain = chain[-burnSamples,]
  # Thinning
  chain <- chain[seq(1,nrow(chain),thinning),]
  # Assign names to columns
  colnames(chain) = names(startValue)
  # MH Output
  return(chain)
}
