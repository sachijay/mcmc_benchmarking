source("functions.R")

## theta - latent variables
## data - observations

#####
## Beta prior
prior <- function(theta, hyper_params){
  
  beta_prior <- dbeta(theta,
                      shape1 = hyper_params[1],
                      shape2 = hyper_params[2])
  
  return(beta_prior)
}

## Bernoulli likelihood
likelihood <- function(data, theta){
  
  if(theta > 1 | theta < 0){
    bernoulli_likelihood <- 0
  } else {
    bernoulli_likelihood <- prod(dbinom(x = data,
                                        size = 1,
                                        prob = theta))
    
    # n_data <- length(data)
    # sum_data <- sum(data)
    # bernoulli_likelihood <- theta^sum_data * (1 - theta)^(n_data - sum_data)
  }
  
  return(bernoulli_likelihood)
}


## Joint distribution
joint_distribution <- function(data, theta, hyper_params){
  
  out <- prior(theta, hyper_params) * likelihood(data, theta)
  
  return(out)
}


#####
set.seed(2021)

data <- c(0, 0, 0)
alpha <- 1
beta <- 1

begin_time <- Sys.time()
mh_sample <- run_mh(data, 
                    n_samples = n_samples,
                    hyper_params = c(alpha, beta),
                    sd = 1)
end_time <- Sys.time()
print(end_time - begin_time)

plot(mh_sample,
     type = 'l')

hist(mh_sample, 
     probability = TRUE)

x <- seq(from = min(mh_sample), 
         to = max(mh_sample), 
         length.out = 1000)
y <- dbeta(x, 
           shape1 = alpha + sum(data), 
           shape2 = beta + length(data) - sum(data))
lines(x, y)
rm(x, y)
