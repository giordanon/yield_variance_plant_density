

# Process model for the mean
quadPlateau <- function(b1, AOPD,x){
  y  = ifelse(x > AOPD,
              b1*AOPD+(-b1/(2*AOPD))*AOPD^2,
              # else
              b1*x+(-b1/(2*AOPD))*x^2)
  
  return(y)
}

#Process Model for the variance
expPlateau <- function(a1, maxVPD,minVPD,x){
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
  
  pred = quadPlateau(b1, AOPD, x) # Process model to the mean
  
  predVar <- expPlateau(a1,maxVPD,minVPD,x) # Process model for the variance
  
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
  # B1
  bprior = dunif(b1, min = 0.001, max = 0.2, log = T)
  #AOPD
  AOPDprior = dunif(AOPD, min = 10, max = 400, log = T)
  # A1
  aprior = dunif(a1, min = 0, max = 0.05, log = T)
  # minVPD
  # Upper bound from Canadian data 192
  minVPDprior = dunif(minVPD, min = 10, max = 192, log = T)
  # maxVPD
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
  b1.try <- rnorm(1, mean = paramOld[1], sd = sigmaTune[1])
  AOPD.try <- rnorm(1, mean = paramOld[2], sd = sigmaTune[2])
  # Variance
  a1.try <- rnorm(1, mean = paramOld[3], sd = sigmaTune[3])
  maxVPD.try <- rnorm(1, mean = paramOld[4], sd = sigmaTune[4])
  minVPD.try <- rnorm(1, mean = paramOld[5], sd = sigmaTune[5])
  
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

# Run Splits

runSplits <- function(data, 
                      groups, 
                      iterations,
                      burnIn,
                      thinning,
                      startValue,
                      sigmaTune,
                      environment
){
  
  
  iterations <- iterations
  burnIn <- burnIn
  thinning <- thinning
  startValue <- startValue
  sigmaTune <- sigmaTune
  
  
  
  filename = paste(groups, collapse = "_")
  # Parallel processing
  cores = ifelse(nrow(unique(data[groups]))>16,16,nrow(unique(data[groups])))
  cluster <- new_cluster(cores)
  cluster_library(cluster, c("dplyr", "purrr"))
  cluster_copy(cluster, environment)
  
  # Set seed
  set.seed(2023)
  
  data %>% 
    group_by_at(groups) %>% 
    nest() %>% 
    multidplyr::partition(cluster) %>% 
    mutate(mod = map(.x = data,
                     ~metropolisMCMC(iterations = iterations, 
                                     burnIn = burnIn, 
                                     thinning = thinning,
                                     startValue = startValue,
                                     x = .x$STAND,
                                     y = .x$YIELD)
    )
    ) %>%
    dplyr::collect() %>% 
    ungroup() %>% 
    saveRDS(paste0("../output/model/model.RData"))
  
  return(paste0("../output/model/model.RData"))
}



# Trace plots ans posteriors
tracePlots <- function(mod){
  
  for (i in 1:nrow(mod)){
    
    out <- mod$mod[[i]]
    
    title <- paste0(mod$YE[i],"_", mod$TRT[i], "_", mod$CLEAN[i])
    filename <- paste0("../output/plots/trace/", title, ".png")
    
    png(filename, width=800, height=350) 
    par(mfrow = c(2,5))
    hist(out[,1],nclass=30, main="b1" )
    hist(out[,2],nclass=30, main="AOPD")
    hist(out[,3],nclass=30, main="a1")
    hist(out[,4],nclass=30, main="maxVPD")
    hist(out[,5],nclass=30, main="minVPD")
    
    plot(out[,1], type = "l",  main = "b1")
    plot(out[,2], type = "l",  main = "AOPD")
    plot(out[,3], type = "l",  main = "a1")
    plot(out[,4], type = "l",  main = "maxVPD")
    plot(out[,5], type = "l",  main = "minVPD")
    title(title, outer = TRUE)
    # Close the PNG device
    dev.off()
    
  }
  
}

# Get model summaries

modelSum <- function(mod, step, xmax){
  
  xseq = seq(1,xmax,step)
  
  
  est_int = matrix(NA,nrow(mod), 2)
  pred_int = matrix(NA,nrow(mod), 1)
  
  
  out = matrix(NA, length(xseq), 9)
  
  system.time({
    for (x in 1:length(xseq)){
      for (i in 1:nrow(mod)){
        
        draw = mod[i,]
        
        a1 = draw[1]
        maxVPD = draw[2]
        minVPD = draw[3]
        b1 = draw[4]
        AOPD = draw[5]
        
        var_pred = expPlateau(a1, maxVPD, minVPD, x)
        mean_pred = quadPlateau(b1,AOPD, x)
        
        pi_draw = rnorm(1, mean_pred, var_pred)
        
        est_int[i,] = c(var_pred, mean_pred)
        pred_int[i,] = pi_draw
      }
      
      # Get summaries for variance
      var_q500 = mean(est_int[,1])
      var_q025 = quantile(est_int[,1], probs = 0.025)
      var_q975 = quantile(est_int[,1], probs = 0.975)
      
      # Get summaries for mean
      m_q500 = mean(est_int[,2])
      m_q025 = quantile(est_int[,2], probs = 0.025)
      m_q975 = quantile(est_int[,2], probs = 0.975)
      
      # Get summaries for the prediction interval
      
      pred_q500 = mean(pred_int)
      pred_q025 = quantile(pred_int, probs = 0.025)
      pred_q975 = quantile(pred_int, probs = 0.975)
      
      
      out[x,] = c(pred_q500, pred_q025,pred_q975, var_q500, var_q025, var_q975, m_q500,m_q025,m_q975)
    }
    
  })
  
  df <- as.data.frame(out) %>% 
    rename_all(~c("pred_q500", "pred_q025","pred_q975", 
                  "var_q500", "var_q025", "var_q975", 
                  "m_q500","m_q025","m_q975")) %>% 
    rowid_to_column("STAND") %>% 
    mutate(STAND = as.double(STAND))
  
  return(df)
}

getSummaries <- function(mod, groups, step){
  
  filename = paste(groups, collapse = "_")
  
  # Parallel processing
  cores = ifelse(nrow(unique(mod[groups]))>16,16,nrow(unique(mod[groups])))
  cluster <- new_cluster(cores)
  cluster_library(cluster, c("dplyr", "purrr", "tibble"))
  cluster_copy(cluster,c("modelSum", "expPlateau", "quadPlateau", "step","groups"))
  
  out <- mod %>% 
    group_by_at(groups) %>% 
    multidplyr::partition(cluster) %>% 
    mutate(mod_sum = pmap(list(..1 = data, ..2 = mod), 
                          ~modelSum(mod = ..2,
                                    step = step, 
                                    xmax = max(..1$STAND) )
    )
    ) %>% 
    dplyr::collect()
  
  
  outpath <- paste0("../output/model/summaries.RData")
  
  
  out %>% saveRDS(outpath)  
  return(out)
}


