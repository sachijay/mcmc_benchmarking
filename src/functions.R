#####
## A Normal proposal distribution: q(x*, xi)
d_proposal <- function(theta_s, 
                       theta_i, 
                       sd = 1){
  
  out <- dnorm(x = theta_s, 
               mean = theta_i,
               sd = sd)
  
  return(out)
}

r_proposal <- function(theta_i, 
                       sd = 1){
  
  out <- rnorm(1, 
               mean = theta_i,
               sd = sd)
  
  return(out)
}


#####
## Run the Metropolis-Hastings algorithm
run_mh <- function(data, sd = 1, 
                   hyper_params = NULL, 
                   n_iter = 10000,
                   ...){
  
  ## initialize vector to save the sampled values
  sampled_vals <- rep(NA, 
                      times = n_iter)
  
  ## initial parameter
  sampled_vals[1] <- 0.5
  
  
  for (j in 1:n_iter) {
    
    theta_i <- sampled_vals[j]
    
    theta_s <- r_proposal(theta_i, 
                          ...)
    
    r <- joint_distribution(data, theta_s, hyper_params) * d_proposal(theta_i, theta_s, ...) / 
      (joint_distribution(data, theta_i, hyper_params) * d_proposal(theta_s, theta_i, ...))
    
    accept <- rbinom(n = 1,
                     size = 1,
                     prob = min(1, r))
    
    sampled_vals[j+1] <- ifelse(accept, 
                                yes = theta_s,
                                no = theta_i)
  }
  
  return(sampled_vals)
}
