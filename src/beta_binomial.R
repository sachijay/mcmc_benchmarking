## Run functions.R
source(here::here("src", "functions.R"))

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
## Run benchmark for a single set of hyper parameters
beta_bin_bench <- function(hyper_params, data, prior_len_out = 1000){
  
  begin_time <- Sys.time()
  mh_sample <- run_mh(data, 
                      n_samples = n_samples,
                      hyper_params = hyper_params,
                      sd = 1)
  end_time <- Sys.time()
  
  ## Calculate the run time for the MCMC algorithm
  run_time <- as.numeric(end_time - begin_time,
                         units = "secs")
  
  ## Calculate posterior mean and MC SE
  mcerr <- mcmcse::mcse(mh_sample)
  
  ## Calculate the effective sample size
  ess <- mcmcse::ess(mh_sample)
  
  ## Get the prior dist
  prior_x <- seq(0, 1, length.out = prior_len_out)
  prior_density <- prior(prior_x,
                         hyper_params = hyper_params)
  
  
  out <- list(param = hyper_params,
              run_time = run_time,
              ess = ess,
              ess_per_time = ess/run_time,
              prior_x = prior_x,
              prior = prior_density,
              pos_sample = mh_sample,
              pos_mean = mcerr$est,
              mcse = mcerr$se)
  
  return(out)
}


#####
set.seed(520)

## Data fixed to be all failures
data <- c(0, 0, 0)

## Define a range for the hyper parameters
## Prior parameters are set to alpha = beta
param_vec <- seq(1, 100, by = 2)

## Run the benchmark for all the parameters in the `param_vec`
beta_bin_bench_vals <- lapply(param_vec, function(param, data){
  
  out <- beta_bin_bench(hyper_params = c(alpha = param, 
                                         beta = param),
                        data = data)
  
  return(out)
  
}, data = data)


#####
## Extract ESS
ess <- sapply(beta_bin_bench_vals, function(x){
  return(x$ess)
})

## Extract running time
run_time <- sapply(beta_bin_bench_vals, function(x){
  return(x$run_time)
})

## Extract ESS per time
ess_per_time <- sapply(beta_bin_bench_vals, function(x){
  x$ess_per_time
})


## Plot the running time vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_time_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = run_time)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = "Run time (s)"
       # title = "Runtime vs parameters for Beta-Binomial"
  )
dev.off()

## Plot ESS vs parameter value
ggplot(mapping = aes(x = param_vec,
                     y = ess)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = "ESS"
       # title = "ESS vs parameters for Beta-Binomial"
  )

## Plot the ESS/time vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_esspt_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = log(ess_per_time))) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = latex2exp::TeX(r"($\log_{10}$ (ESS per second) )")
       # title = "ESS per sec vs parameters for Beta-Binomial"
  )
dev.off()


#####

## Plot end 2 posteriors and priors
s1 <- beta_bin_bench_vals[[2]]
x <- s1$prior_x
y_theoretical <- dbeta(x, 
                       shape1 = s1$param["alpha"] + sum(data), 
                       shape2 = s1$param["beta"] + length(data) - sum(data))

y_prior <- s1$prior


pdf(file = paste0(fig_dir,
                  "beta_bin_s1_dist.pdf"))
ggplot(mapping = aes(x = s1$pos_sample,
                     y = after_stat(density))) + 
  geom_histogram() +
  geom_line(mapping = aes(x = x,
                          y = y_theoretical,
                          colour = "darkblue"), 
            size = 1) + 
  geom_line(mapping = aes(x = x,
                          y = y_prior,
                          colour = "red")) + 
  scale_color_discrete(name = "Distribution", 
                       labels = c("Posterior", 
                                  "Prior")) +
  theme_minimal() + 
  labs(x = "Posterior sample",
       y = "Density", 
       title = latex2exp::TeX(sprintf(r"($\alpha = %d, \beta = %d$)", s1$param["alpha"], s1$param["beta"])))
dev.off()

s2 <- beta_bin_bench_vals[[length(beta_bin_bench_vals)]]
x <- s2$prior_x
y_theoretical <- dbeta(x, 
                       shape1 = s2$param["alpha"] + sum(data), 
                       shape2 = s2$param["beta"] + length(data) - sum(data))

y_prior <- s2$prior

pdf(file = paste0(fig_dir,
                  "beta_bin_s99_dist.pdf"))
ggplot(mapping = aes(x = s2$pos_sample,
                     y = after_stat(density))) + 
  geom_histogram() +
  geom_histogram() +
  geom_line(mapping = aes(x = x,
                          y = y_theoretical,
                          colour = "darkblue"), 
            size = 1) + 
  geom_line(mapping = aes(x = x,
                          y = y_prior,
                          colour = "red")) + 
  scale_color_discrete(name = "Distribution", 
                       labels = c("Posterior", 
                                  "Prior")) +
  theme_minimal() + 
  labs(x = "Posterior sample",
       y = "Density",
       title = latex2exp::TeX(sprintf(r"($\alpha = %d, \beta = %d$)", s2$param["alpha"], s2$param["beta"])))
dev.off()