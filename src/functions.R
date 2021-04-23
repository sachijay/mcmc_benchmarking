## Run the main script
source(here::here("src", "main.R"))

#####
## A normal proposal distribution: q(x*, xi)
d_proposal <- function(x_s, 
                       x_i, 
                       sd = 1){
  
  out <- dnorm(x = x_s, 
               mean = x_i,
               sd = sd)
  
  return(out)
}

r_proposal <- function(theta_i, 
                       sd = 1){
  
  out <- rnorm(length(theta_i), 
               mean = theta_i,
               sd = sd)
  
  return(out)
}


#####
## Run the Metropolis-Hastings algorithm
run_mh <- function(data, 
                   n_samples,
                   init_vals,
                   hyper_params = NULL,
                   burn = TRUE,
                   burn_prop = 0.5, 
                   annealing = 1,
                   ...){
  
  n_unknown <- length(init_vals)
  
  ## Initialize vector to save the sampled values
  sampled_vals <- matrix(NA, 
                         nrow = n_samples+1,
                         ncol = n_unknown)
  
  ## Initial parameter
  sampled_vals[1,] <- init_vals
  
  
  for (j in 1:n_samples) {
    
    x_i <- sampled_vals[j,]
    
    x_s <- r_proposal(x_i, 
                      ...)
    
    ## add the proposals if using a non-symmetric proposal
    # r <- joint_distribution(data, x_s, hyper_params) * d_proposal(x_i, x_s, ...) / 
    #   (joint_distribution(data, x_i, hyper_params) * d_proposal(x_s, x_i, ...))
    r <- joint_distribution(data, x_s, hyper_params, annealing) / 
      joint_distribution(data, x_i, hyper_params, annealing)
    
    accept <- rbinom(n = 1,
                     size = 1,
                     prob = min(1, r))
    
    ## For debugging
    # print("x_i = ")
    # print(x_i)
    # 
    # print("x_s = ")
    # print(x_s)
    
    # print("pi_i = ")
    # print(joint_distribution(data, x_i, hyper_params, annealing))
    # 
    # print("pi_s = ")
    # print(joint_distribution(data, x_s, hyper_params, annealing))

    # print("r = ")
    # print(r)
    # 
    # print("accept = ")
    # print(accept)
    
    
    if(accept){
      sampled_vals[j+1, ] <- x_s
    } else {
      sampled_vals[j+1, ] <- x_i
    }
  }
  
  if(burn){
    ## Number burn-in
    n_burnin <- ceiling(burn_prop * n_samples)
    
    out <- sampled_vals[(n_burnin+1):n_samples, ]
  } else {
    out <- sampled_vals
  }
  
  return(out)
}


#####
## Run the parallel tempering algorithm
run_pt <- function(data,
                   n_samples,
                   n_chains,
                   init_vals,
                   hyper_params = NULL,
                   burn = TRUE,
                   burn_prop = 0.5,
                   sd = 1,
                   ...){
  
  annealing_params <- seq(0, 1, length.out = n_chains)
  
  ## Initialize vector to save the sampled values
  sampled_vals <- matrix(NA,
                         nrow = n_samples + 1,
                         ncol = n_chains)
  # sampled_vals <- rep(NA,
  #                     times = n_samples + 1)
  
  ## Initial parameter
  sampled_vals[1,] <- rep(init_vals, times = n_chains)
  
  
  for (j in 1:n_samples) {
    
    proposal <- rep(NA, times = n_chains)
    
    for (ann in 1:n_chains) {
      
      ann_param <- annealing_params[ann]
      
      proposal[ann] <- head(run_mh(data,
                                   n_samples = 1,
                                   hyper_params = hyper_params,
                                   init_vals = sampled_vals[j, ann],
                                   annealing = ann_param,
                                   sd = sd), 1)
    }
    
    swap_chain <- sample(1:(n_chains - 1), size = 1)
    
    x_i <- proposal[swap_chain]
    x_i1 <- proposal[swap_chain + 1]
    
    ann_i <- annealing_params[swap_chain]
    ann_i1 <- annealing_params[swap_chain + 1]
    
    r <- joint_distribution(data, x_i1, hyper_params, ann_i) * joint_distribution(data, x_i, hyper_params, ann_i1) / 
      joint_distribution(data, x_i, hyper_params, ann_i) * joint_distribution(data, x_i1, hyper_params, ann_i1)
    
    accept <- rbinom(n = 1,
                     size = 1,
                     prob = min(1, r))
    
    ## For debugging
    # print("x_i = ")
    # print(x_i)
    # 
    # print("x_s = ")
    # print(x_s)
    # 
    # print("r = ")
    # print(r)
    # 
    # print("accept = ")
    # print(accept)
    #
    # print("swap_chain = ")
    # print(swap_chain)
    
    
    if(accept){
      proposal[swap_chain] <- x_i1
      proposal[swap_chain + 1] <- x_i
    }
    
    sampled_vals[j+1, ] <- proposal
    
  }
  
  if(burn){
    ## Number burn-in
    n_burnin <- ceiling(burn_prop * n_samples)
    
    out <- sampled_vals[(n_burnin+1):n_samples, n_chains]
  } else {
    out <- sampled_vals[,n_chains]
  }
  
  return(out)
}
