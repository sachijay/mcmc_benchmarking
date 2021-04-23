##############################################################
## This file runs simulations for the Beta-Binomial problem ##
##############################################################

#####
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
joint_distribution <- function(data, theta, hyper_params, annealing = 1){
  
  out <- prior(theta, hyper_params) * likelihood(data, theta)^annealing
  
  return(out)
}

#####
## Run benchmark for a single set of hyper parameters
beta_bin_bench <- function(hyper_params, data, n_samples, prior_len_out = 1000){
  
  ## Run the MH algorithm
  begin_time <- Sys.time()
  mh_sample <- run_mh(data, 
                      n_samples = n_samples,
                      hyper_params = hyper_params,
                      init_vals = 0.5,
                      sd = 1)
  end_time <- Sys.time()
  
  ## Calculate the run time for the MH algorithm
  mh_run_time <- as.numeric(end_time - begin_time,
                            units = "secs")
  
  ## Calculate posterior mean and MC SE
  mh_mcerr <- mcmcse::mcse(mh_sample)
  
  
  ## Calculate the effective sample size
  mh_ess <- mcmcse::ess(mh_sample)
  
  begin_time <- Sys.time()
  pt_sample <- run_pt(data,
                      n_samples = n_samples,
                      n_chains = 5,
                      hyper_params = hyper_params,
                      init_vals = 0.5,
                      sd = 1)
  end_time <- Sys.time()
  
  ## Calculate the run time for the PT algorithm
  pt_run_time <- as.numeric(end_time - begin_time,
                            units = "secs")
  
  ## Calculate posterior mean and MC SE
  pt_mcerr <- mcmcse::mcse(pt_sample)
  
  ## Calculate the effective sample size
  pt_ess <- mcmcse::ess(pt_sample)
  
  ## Get the prior dist
  prior_x <- seq(0, 1, length.out = prior_len_out)
  prior_density <- prior(prior_x,
                         hyper_params = hyper_params)
  
  
  out <- list(param = hyper_params,
              mh_run_time = mh_run_time,
              mh_ess = mh_ess,
              mh_ess_per_time = mh_ess/mh_run_time,
              mh_sample = mh_sample,
              mh_pos_mean = mh_mcerr$est,
              mh_mcse = mh_mcerr$se,
              pt_run_time = pt_run_time,
              pt_ess = pt_ess,
              pt_ess_per_time = pt_ess/pt_run_time,
              pt_sample = pt_sample,
              pt_pos_mean = pt_mcerr$est,
              pt_mcse = pt_mcerr$se,
              prior_x = prior_x,
              prior = prior_density)
  
  return(out)
}

#####
set.seed(2021)

## Data fixed to be all failures
data <- c(0, 0, 0)

## Define a range for the hyper parameters
## Prior parameters are set to alpha = beta
param_vec <- seq(1, 100, by = 2)

## Run the benchmark for all the parameters in the `param_vec`
# beta_bin_bench_vals <- lapply(param_vec, function(param, data){
# 
#   out <- beta_bin_bench(hyper_params = c(alpha = param,
#                                          beta = param),
#                         n_samples = n_samples,
#                         data = data)
# 
#   return(out)
# 
# }, data = data)

## Save and load the benchmark for convenience
beta_bin_data_file <- paste0(data_dir,
                             "beta_bin_bench_vals.Rdata")
# save(beta_bin_bench_vals,
#      file = beta_bin_data_file)
load(file = beta_bin_data_file)
rm(beta_bin_data_file)


#####
## Define a range for the number of samples
n_samples_vec <- seq(500, n_samples, by = 500)

## Run the benchmark for all the parameters in the `n_samples_vec`
# beta_bin_bench_vals2 <- lapply(n_samples_vec, function(n_samples, data){
# 
#   param <- tail(param_vec, 1)
#   
#   out <- beta_bin_bench(hyper_params = c(param, param),
#                         n_samples = n_samples,
#                         data = data)
#   
#   return(out)
#   
# }, data = data)

## Save and load the benchmark for convenience
beta_bin_data_file2 <- paste0(data_dir,
                              "beta_bin_bench_vals2.Rdata")
# save(beta_bin_bench_vals2,
#      file = beta_bin_data_file2)
load(file = beta_bin_data_file2)
rm(beta_bin_data_file2)


##########################################
## Visualization of prior and posterior ##
##########################################

s1 <- beta_bin_bench_vals[[2]]
x <- s1$prior_x
y_theoretical <- dbeta(x, 
                       shape1 = s1$param["alpha"] + sum(data), 
                       shape2 = s1$param["beta"] + length(data) - sum(data))

y_prior <- s1$prior

pdf(file = paste0(fig_dir,
                  "beta_bin_mh_s1_dist.pdf"))
ggplot(mapping = aes(x = s1$mh_sample,
                     y = after_stat(density))) + 
  geom_histogram() +
  geom_line(mapping = aes(x = x,
                          y = y_theoretical,
                          colour = "darkblue"), 
            size = 1) + 
  geom_line(mapping = aes(x = x,
                          y = y_prior,
                          colour = "red"), 
            size = 1) + 
  scale_color_discrete(name = "Distribution", 
                       labels = c("Posterior", 
                                  "Prior")) +
  theme_minimal() + 
  labs(x = "Posterior sample",
       y = "Density", 
       title = latex2exp::TeX(sprintf(r"($\alpha = %d, \beta = %d$)", s1$param["alpha"], s1$param["beta"])))
dev.off()

rm(s1, x, y_theoretical, y_prior)

s2 <- beta_bin_bench_vals[[length(beta_bin_bench_vals)]]
x <- s2$prior_x
y_theoretical <- dbeta(x, 
                       shape1 = s2$param["alpha"] + sum(data), 
                       shape2 = s2$param["beta"] + length(data) - sum(data))

y_prior <- s2$prior

pdf(file = paste0(fig_dir,
                  "beta_bin_mh_s99_dist.pdf"))
ggplot(mapping = aes(x = s2$mh_sample,
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

rm(s2, x, y_theoretical, y_prior)


####################################
## Visualization of MH benchmarks ##
####################################

####
## vs. parameter value

#####
## Extract ESS
mh_ess <- sapply(beta_bin_bench_vals, function(x){
  return(x$mh_ess)
})

## Extract running time
mh_run_time <- sapply(beta_bin_bench_vals, function(x){
  return(x$mh_run_time)
})

## Extract ESS per time
mh_ess_per_time <- sapply(beta_bin_bench_vals, function(x){
  x$mh_ess_per_time
})


## Plot the running time vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_mh_time_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = mh_run_time)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = "Run time (s)"
       # title = "Runtime vs parameters for Beta-Binomial"
  )
dev.off()

## Plot ESS vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_mh_ess_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = mh_ess)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = "ESS"
       # title = "ESS vs parameters for Beta-Binomial"
  )
dev.off()

## Plot the ESS/time vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_mh_esspt_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = log10(mh_ess_per_time))) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = latex2exp::TeX(r"($\log_{10}$ (ESS per second) )")
       # title = "ESS per sec vs parameters for Beta-Binomial"
  )
dev.off()


####
## vs. number of samples

mh_ess <- sapply(beta_bin_bench_vals2, function(x){
  return(x$mh_ess)
})

## Extract running time
mh_run_time <- sapply(beta_bin_bench_vals2, function(x){
  return(x$mh_run_time)
})

## Extract ESS per time
mh_ess_per_time <- sapply(beta_bin_bench_vals2, function(x){
  x$mh_ess_per_time
})

## Plot the running time vs number of samples
pdf(file = paste0(fig_dir,
                  "beta_bin_mh_time_samples.pdf"))
ggplot(mapping = aes(x = log10(n_samples_vec),
                     y = mh_run_time)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\log_{10}$ (Samples) )"),
       y = "Run time (s)"
       # title = "Runtime vs number of samples for Beta-Binomial"
  )
dev.off()

## Plot ESS vs number of samples
pdf(file = paste0(fig_dir,
                  "beta_bin_mh_ess_samples.pdf"))
ggplot(mapping = aes(x = log10(n_samples_vec),
                     y = mh_ess)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\log_{10}$ (Samples) )"),
       y = "ESS"
       # title = "ESS vs number of samples for Beta-Binomial"
  )
dev.off()

## Plot the ESS/time vs number of samples
pdf(file = paste0(fig_dir,
                  "beta_bin_mh_esspt_samples.pdf"))
ggplot(mapping = aes(x = log10(n_samples_vec),
                     y = log10(mh_ess_per_time))) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\log_{10}$ (Samples) )"),
       y = latex2exp::TeX(r"($\log_{10}$ (ESS per second) )")
       # title = "ESS per sec vs number of samples for Beta-Binomial"
  )
dev.off()

rm(mh_ess, mh_run_time, mh_ess_per_time)


####################################
## Visualization of PT benchmarks ##
####################################

#####
## Extract ESS
pt_ess <- sapply(beta_bin_bench_vals, function(x){
  return(x$pt_ess)
})

## Extract running time
pt_run_time <- sapply(beta_bin_bench_vals, function(x){
  return(x$pt_run_time)
})

## Extract ESS per time
pt_ess_per_time <- sapply(beta_bin_bench_vals, function(x){
  x$pt_ess_per_time
})


## Plot the running time vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_pt_time_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = pt_run_time)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = "Run time (s)"
       # title = "Runtime vs parameters for Beta-Binomial"
  )
dev.off()

## Plot ESS vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_pt_ess_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = pt_ess)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = "ESS"
       # title = "ESS vs parameters for Beta-Binomial"
  )
dev.off()

## Plot the ESS/time vs parameter value
pdf(file = paste0(fig_dir,
                  "beta_bin_pt_esspt_params.pdf"))
ggplot(mapping = aes(x = param_vec,
                     y = log10(pt_ess_per_time))) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\alpha = \beta$)"),
       y = latex2exp::TeX(r"($\log_{10}$ (ESS per second) )")
       # title = "ESS per sec vs parameters for Beta-Binomial"
  )
dev.off()


####
## vs. number of samples

pt_ess <- sapply(beta_bin_bench_vals2, function(x){
  return(x$pt_ess)
})

## Extract running time
pt_run_time <- sapply(beta_bin_bench_vals2, function(x){
  return(x$pt_run_time)
})

## Extract ESS per time
pt_ess_per_time <- sapply(beta_bin_bench_vals2, function(x){
  x$pt_ess_per_time
})

## Plot the running time vs number of samples
pdf(file = paste0(fig_dir,
                  "beta_bin_pt_time_samples.pdf"))
ggplot(mapping = aes(x = log10(n_samples_vec),
                     y = pt_run_time)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\log_{10}$ (Samples) )"),
       y = "Run time (s)"
       # title = "Runtime vs number of samples for Beta-Binomial"
  )
dev.off()

## Plot ESS vs number of samples
pdf(file = paste0(fig_dir,
                  "beta_bin_pt_ess_samples.pdf"))
ggplot(mapping = aes(x = log10(n_samples_vec),
                     y = pt_ess)) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\log_{10}$ (Samples) )"),
       y = "ESS"
       # title = "ESS vs number of samples for Beta-Binomial"
  )
dev.off()

## Plot the ESS/time vs number of samples
pdf(file = paste0(fig_dir,
                  "beta_bin_pt_esspt_samples.pdf"))
ggplot(mapping = aes(x = log10(n_samples_vec),
                     y = log10(pt_ess_per_time))) + 
  geom_line() +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\log_{10}$ (Samples) )"),
       y = latex2exp::TeX(r"($\log_{10}$ (ESS per second) )")
       # title = "ESS per sec vs number of samples for Beta-Binomial"
  )
dev.off()

rm(pt_ess, pt_run_time, pt_ess_per_time)


#####
## Remove all variables from workspace
rm(list = ls())
