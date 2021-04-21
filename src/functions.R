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
  
  out <- rnorm(1, 
               mean = theta_i,
               sd = sd)
  
  return(out)
}


#####
## Run the Metropolis-Hastings algorithm
run_mh <- function(data, 
                   n_samples,
                   hyper_params = NULL,
                   burn = TRUE,
                   burn_prop = 0.5,
                   ...){
  
  ## Initialize vector to save the sampled values
  sampled_vals <- rep(NA, 
                      times = n_samples)
  
  ## Initial parameter
  sampled_vals[1] <- 0.5
  
  
  for (j in 1:n_samples) {
    
    x_i <- sampled_vals[j]
    
    x_s <- r_proposal(x_i, 
                      ...)
    
    
    ## add the proposals if using a non-symmetric proposal
    # r <- joint_distribution(data, x_s, hyper_params) * d_proposal(x_i, x_s, ...) / 
    #   (joint_distribution(data, x_i, hyper_params) * d_proposal(x_s, x_i, ...))
    r <- joint_distribution(data, x_s, hyper_params) / 
      joint_distribution(data, x_i, hyper_params)
    
    
    accept <- rbinom(n = 1,
                     size = 1,
                     prob = min(1, r))
    
    sampled_vals[j+1] <- ifelse(accept, 
                                yes = x_s,
                                no = x_i)
  }
  
  if(burn){
    ## Number burn-in
    n_burnin <- ceiling(burn_prop * n_samples)
    
    out <- sampled_vals[(n_burnin+1):n_samples]
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
                   hyper_params = NULL,
                   burn = TRUE,
                   burn_prop = 0.5,
                   ...){
  
  n_data <- length(data)
  
  b <- 0:n_data/n_data
  
  return(b)
}

run_pt(data, 
       n_samples = n_samples,
       n_chains = 10)


#######################################################
# temper<-function(niter, 
#                  Bmin,
#                  swap.interval,
#                  rank,
#                  size){
#   
#   swap = 0
#   swaps.attempted = 0
#   swaps.accepted = 0
#   
#   # Higher ranks run the higher "temperatures" (~smaller fractional powers)
#   B = rep(0, size-1)
#   
#   for(r in 1:size-1){
#     temp = (r-1)/(size-2)
#     B[r] = Bmin^temp
#   }
#   
#   
#   # Create a list for proposal moves
#   prop = rep(0,2);
#   theta = matrix(0,niter,2)
#   
#   for(t in 2:niter){
#     
#     for(c in 1:length(prop))   
#       prop = theta[t-1, c] + rnorm(1,0,0.1);
#     
#     #Calculate Log-Density at proposed and current position
#     logdensity.current = logdensity(theta[t-1,])
#     logdensity.prop = logdensity(prop);
#     
#     #Calculate log acceptance probability
#     lalpha = B[rank]*(logdensity.prop-logdensity.current)
#     
#     if(log(runif(1))<lalpha){
#       #Accept proposed move
#       theta[t,] = prop;
#       logdensity.current = logdensity.prop;
#     }else{
#       #Otherwise do not move
#       theta[t,]=theta[t-1,];
#     }
#     
#     if(t%%swap.interval ==0){
#       for(evenodd in 0:1){
#         swap=0;
#         logdensity.partner=0;
#         if(rank%%2 == evenodd%%2){
#           rank.partner=rank + 1;
#           #ranks range from 1:size-1. Cannot have a partner rank == size
#           if(0<rank.partner && rank.partner<size){
#             #On first iteration, evens receive from above odd
#             #On second iteration, odds receive from above evens
#             logdensity.partner<-mpi.recv.Robj(rank.partner,rank.partner);
#             lalpha = (B[rank]-B[rank.partner])*(logdensity.partner-logdensity.current);
#             swaps.attempted=swaps.attempted+1;
#             if(log(runif(1))<lalpha){
#               swap=1;
#               swaps.accepted=swaps.accepted+1;
#             }
#             mpi.send.Robj(swap,dest=rank.partner,tag=rank)
#           }
#           if(swap==1){
#             thetaswap=theta[t,];
#             mpi.send.Robj(thetaswap,dest=rank.partner,tag=rank)
#             theta[t,]=mpi.recv.Robj(rank.partner,rank.partner)
#           }
#         }else{
#           rank.partner=rank-1;
#           #ranks range from 1:size-1. Cannot have a partner rank ==0
#           if(0<rank.partner && rank.partner<size){
#             #On first iteration, odds send to evens below
#             #On second iteration, evens sent to odds below
#             mpi.send.Robj(logdensity.current,dest=rank.partner,tag=rank);
#             swap=mpi.recv.Robj(rank.partner,rank.partner);
#           }
#           if(swap==1){
#             thetaswap=theta[t,];
#             theta[t,]=mpi.recv.Robj(rank.partner,rank.partner);
#             mpi.send.Robj(thetaswap,dest=rank.partner,tag=rank);
#           }
#         }
#       }
#     }
#   }
#   return(theta)
# }
