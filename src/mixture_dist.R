#####################################################################
## This file runs simulations for the mixture distribution problem ##
#####################################################################


#####
## Run functions.R
source(here::here("src", "functions.R"))


#####
p <- 0.3
mu_1 <- 2
mu_2 <- 10
sd_1 <- 2
sd_2 <- 2
n <- 1000


#####
set.seed(520)

ind <- rbinom(n, 
              size = 1, 
              prob = p)
data <- ifelse(ind, 
               yes = rnorm(n, mean = mu_1, sd = sd_1),
               no = rnorm(n, mean = mu_2, sd = sd_2))

x <- seq(-5, 15, length.out = 1000)
y <- p*dnorm(x, mean = mu_1, sd = sd_1) + (1-p)*dnorm(x, mean = mu_2, sd = sd_2)

ggplot() + 
  geom_histogram(mapping = aes(x = data,
                               y = after_stat(density))) + 
  geom_line(mapping = aes(x = x,
                          y = y))


####
# library("rjags")
# jags_model <- "
#   model {
#     for (i in 1:n) {
#       y[i] ~ dnorm(mu[zeta[i]], tau[zeta[i]])
#       zeta[i] ~ dcat(pi[])
#     }
#     
#     for (i in 1:H) {
#       mu[i] ~ dnorm(0,1e-5)
#       tau[i] ~ dgamma(1,1)
#       sigma[i] <- 1/sqrt(tau[i])
#     }
#     
#     pi ~ ddirich(a)
# }"
# 
# 
# dat <- list(n = length(data), H = 2, y = data, a = rep(1,2))
# jm <- jags.model(textConnection(jags_model), data = dat, n.chains = 1)
# r <- coda.samples(jm, c('mu','sigma','pi'), 1e3)
# 
# str(r)
# hist(r[[1]][,"pi[1]"])

#####

## Joint distribution
## theta = (mu_1, mu_2, p)
joint_distribution <- function(data, theta, hyper_params){
  
  mu_1 <- theta[1]
  mu_2 <- theta[2]
  # p <- theta[3]
  p <- p

  
  out <- dnorm(mu_1, 
               mean = hyper_params["p_mu_1"], 
               sd = hyper_params["p_sd_1"]) * 
    dnorm(mu_2, 
          mean = hyper_params["p_mu_2"], 
          sd = hyper_params["p_sd_2"]) * 
    dbeta(p, 
          shape1 = hyper_params["p_alpha"], 
          shape2 = hyper_params["p_beta"]) *
    sum(p * dnorm(data, 
                  mean = mu_1, 
                  sd = hyper_params["sd_1"]) + 
       (1 - p) * dnorm(data, 
                       mean = mu_2, 
                       sd = hyper_params["sd_2"]))
  
  # print(paste("p = ", p, " out = ", out))
  # print(dbeta(p, shape1 = hyper_params["p_alpha"], shape2 = hyper_params["p_beta"]))
  
  if(out < 0){
    out <- 0
  }
  
  return(out)
}

hyper_params <- c(p_mu_1 = 1, 
                  p_mu_2 = 11,
                  p_sd_1 = 2,
                  p_sd_2 = 2,
                  sd_1 = sd_1,
                  sd_2 = sd_2,
                  p_alpha = 3,
                  p_beta = 2)

a <- run_mh(data, 
            init_vals = c(0, 0),
            n_samples = n_samples, 
            hyper_params = hyper_params,
            sd = 0.5)

mcerr <- mcmcse::mcse(a)

## Calculate the effective sample size
ess <- mcmcse::ess(a)

colMeans(a)

matplot(a, type = 'l')

hist(a[,1], probability = TRUE)
hist(a[,2], probability = TRUE)
hist(a[,3], probability = TRUE)
