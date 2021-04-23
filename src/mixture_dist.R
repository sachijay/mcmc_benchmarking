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

## Joint distribution
## theta = (mu_1, mu_2, p)
joint_distribution <- function(data, theta, hyper_params, annealing = 1){
  
  mu_1 <- theta[1]
  mu_2 <- theta[2]
  p <- theta[3]
  
  prior_mu_1 <- dnorm(mu_1, 
                      mean = hyper_params["p_mu_1"], 
                      sd = hyper_params["p_sd_1"])
  
  prior_mu_2 <- dnorm(mu_2, 
                      mean = hyper_params["p_mu_2"], 
                      sd = hyper_params["p_sd_2"])
  
  prior_p <- dbeta(p, 
                   shape1 = hyper_params["p_alpha"],
                   shape2 = hyper_params["p_beta"])
  
  likelihood <- sum(p * dnorm(data, 
                              mean = mu_1, 
                              sd = hyper_params["sd_1"]) + 
                      (1 - p) * dnorm(data, 
                                      mean = mu_2, 
                                      sd = hyper_params["sd_2"]))
  
  ## For debugging
  # print("prior_mu_1 = ")
  # print(prior_mu_1)
  # 
  # print("prior_mu_2 = ")
  # print(prior_mu_2)
  # 
  # print("likelihood = ")
  # print(likelihood)
  
  out <- prior_mu_1 * prior_mu_2 * prior_p * likelihood^annealing
  
  # print(paste("p = ", p, " out = ", out))
  # print(dbeta(p, shape1 = hyper_params["p_alpha"], shape2 = hyper_params["p_beta"]))
  
  return(out)
}


############################
## Sample data generation ##
############################

set.seed(2021)

ind <- rbinom(n, 
              size = 1, 
              prob = p)
data <- ifelse(ind, 
               yes = rnorm(n, mean = mu_1, sd = sd_1),
               no = rnorm(n, mean = mu_2, sd = sd_2))

rm(ind)


#####################################
## Visualization of generated data ##
#####################################

x <- seq(-5, 15, length.out = 1000)
y <- p*dnorm(x, mean = mu_1, sd = sd_1) + (1-p)*dnorm(x, mean = mu_2, sd = sd_2)

pdf(file = paste0(fig_dir,
                  "mixture_sample_data.pdf"))
ggplot() + 
  geom_histogram(mapping = aes(x = data,
                               y = after_stat(density))) + 
  geom_line(mapping = aes(x = x,
                          y = y)) +
  theme_minimal() + 
  labs(x = "Simulation sample",
       y = "Density")
dev.off()

rm(x, y)


########################
## Run the simulation ##
########################

hyper_params <- c(p_mu_1 = 1, 
                  p_mu_2 = 11,
                  p_sd_1 = 2,
                  p_sd_2 = 2,
                  sd_1 = sd_1,
                  sd_2 = sd_2,
                  p_alpha = 3,
                  p_beta = 2)

mu_1_vec <- -5:2
mu_2_vec <- 10:3

# mixture_bench_vals <- lapply(1:length(mu_1_vec), function(j){
# 
#   m1 <- mu_1_vec[j]
#   m2 <- mu_2_vec[j]
# 
#   x <- seq(-5, 15, length.out = 1000)
#   p_mu_1 <- dnorm(x, mean = hyper_params["p_mu_1"], sd = hyper_params["p_sd_1"])
#   p_mu_2 <- dnorm(x, mean = hyper_params["p_mu_2"], sd = hyper_params["p_sd_2"])
# 
#   ind <- rbinom(n,
#                 size = 1,
#                 prob = p)
#   data <- ifelse(ind,
#                  yes = rnorm(n, mean = m1, sd = sd_1),
#                  no = rnorm(n, mean = m2, sd = sd_2))
# 
#   rm(ind)
# 
#   ## Run the MH algorithm
#   begin_time <- Sys.time()
#   mh_sample <- run_mh(data,
#                       init_vals = c(1, 1, 0.1),
#                       n_samples = n_samples,
#                       hyper_params = hyper_params,
#                       sd = 0.5)
#   end_time <- Sys.time()
# 
# 
#   ## Calculate the run time for the MH algorithm
#   mh_run_time <- as.numeric(end_time - begin_time,
#                             units = "secs")
# 
#   ## Calculate posterior mean and MC SE
#   mh_mcerr <- mcmcse::mcse(mh_sample)
# 
# 
#   ## Calculate the effective sample size
#   mh_ess <- mcmcse::multiESS(mh_sample)
# 
# 
#   out <- list(param = c(mu_1 = m1, mu_2 = m2),
#               p_x = x,
#               p_mu_1 = p_mu_1,
#               p_mu_2 = p_mu_2,
#               mh_run_time = mh_run_time,
#               mh_mcerr = mh_mcerr,
#               mh_ess = mh_ess,
#               mh_ess_per_time = mh_ess/mh_run_time,
#               mh_mcerr = mh_mcerr,
#               mh_sample = mh_sample)
# 
#   return(out)
# })

## Save and load the benchmark for convenience
mixture_data_file <- paste0(data_dir,
                             "mixture_bench_vals.Rdata")
# save(mixture_bench_vals,
#      file = mixture_data_file)
load(file = mixture_data_file)
rm(mixture_data_file)


##########################################
## Visualization of prior and posterior ##
##########################################

s1 <- mixture_bench_vals[[1]]

pdf(file = paste0(fig_dir,
                  "mixture_mh_s1_dist.pdf"))
ggplot() + 
  geom_histogram(mapping = aes(x = s1$mh_sample,
                               y = after_stat(density))) + 
  geom_line(mapping = aes(x = s1$p_x,
                          y = s1$p_mu_1,
                          colour = "darkblue"),
            size = 1) +
  geom_line(mapping = aes(x = s1$p_x,
                          y = s1$p_mu_2,
                          colour = "red"),
            size = 1) +
  scale_color_discrete(name = "Prior",
                       # labels = c(latex2exp::TeX(r"($\mu_1$)"),
                       #            latex2exp::TeX(r"($\mu_2$)")),
                       labels = c("Mean 1", "Mean 2")) +
  theme_minimal() +
  labs(x = "Posterior sample",
       y = "Density",
       title = latex2exp::TeX(sprintf(r"(Data means: $\mu_1 = %d, \mu_2 = %d$)", s1$param["mu_1"], s1$param["mu_2"])))
dev.off()

rm(s1)


s2 <- mixture_bench_vals[[length(mixture_bench_vals)]]

pdf(file = paste0(fig_dir,
                  "mixture_mh_s2_dist.pdf"))
ggplot() + 
  geom_histogram(mapping = aes(x = s2$mh_sample,
                               y = after_stat(density))) + 
  geom_line(mapping = aes(x = s2$p_x,
                          y = s2$p_mu_1,
                          colour = "darkblue"),
            size = 1) +
  geom_line(mapping = aes(x = s2$p_x,
                          y = s2$p_mu_2,
                          colour = "red"),
            size = 1) +
  scale_color_discrete(name = "Prior",
                       # labels = c(latex2exp::TeX(r"($\mu_1$)"),
                       #            latex2exp::TeX(r"($\mu_2$)")),
                       labels = c("Mean 1", "Mean 2")) +
  theme_minimal() +
  labs(x = "Posterior sample",
       y = "Density",
       title = latex2exp::TeX(sprintf(r"(Data means: $\mu_1 = %d, \mu_2 = %d$)", s2$param["mu_1"], s2$param["mu_2"])))
dev.off()

rm(s2)


####################################
## Visualization of MH benchmarks ##
####################################

####
## vs. parameter value

#####
## Extract ESS
mh_ess <- sapply(mixture_bench_vals, function(x){
  return(x$mh_ess)
})

## Extract running time
mh_run_time <- sapply(mixture_bench_vals, function(x){
  return(x$mh_run_time)
})

## Extract ESS per time
mh_ess_per_time <- sapply(mixture_bench_vals, function(x){
  x$mh_ess_per_time
})

mean_dist <- mu_2_vec - mu_1_vec


## Plot the running time vs parameter value
pdf(file = paste0(fig_dir,
                  "mixture_mh_time_diff.pdf"))
ggplot() + 
  geom_line(mapping = aes(x = mean_dist, 
                          y = mh_run_time)) +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\mu_2 - \mu_1$)"),
       y = "Run time (s)"
       # title = "Runtime vs difference in means for the mixture"
  )
dev.off()

## Plot ESS vs parameter value
pdf(file = paste0(fig_dir,
                  "mixture_mh_ess_diff.pdf"))
ggplot() + 
  geom_line(mapping = aes(x = mean_dist, 
                          y = mh_ess)) +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\mu_2 - \mu_1$)"),
       y = "ESS"
       # title = "Runtime vs difference in means for the mixture"
  )
dev.off()

## Plot the ESS/time vs parameter value
pdf(file = paste0(fig_dir,
                  "mixture_mh_esspt_diff.pdf"))
ggplot() + 
  geom_line(mapping = aes(x = mean_dist, 
                          y = log10(mh_ess_per_time))) +
  theme_minimal() + 
  labs(x = latex2exp::TeX(r"($\mu_2 - \mu_1$)"),
       y = latex2exp::TeX(r"($\log_{10}$ (ESS per second) )")
       # title = "Runtime vs difference in means for the mixture"
  )
dev.off()

rm(mh_ess, mh_run_time, mh_ess_per_time)


#####
## Remove all variables from workspace
rm(list = ls())


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




