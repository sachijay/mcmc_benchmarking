## Number of iterations used by default
n_iter <- 10000

## Burn-in proportion
burn_prop <- 0.5


# #######
# # specify the data, to be used in the likelihood function.
# my_data <- c(rep(0, 6), rep(1, 14))
# 
# # define the relative probability of the target distribution, 
# # as a function of vector theta. for our application, this
# # target distribution is the unnormalized posterior distribution
# target_rel_prob <- function(theta, data) {
#   target_rel_prob <- likelihood(theta, data) * prior_d(theta)
#   return(target_rel_prob)
# }
# 
# # specify the length of the trajectory, i.e., the number of jumps to try:
# traj_length <- 50000 # this is just an arbitrary large number
# 
# # initialize the vector that will store the results
# trajectory <- rep(0, traj_length)
# 
# # specify where to start the trajectory:
# trajectory[1] <- 0.01 # another arbitrary value
# 
# # specify the burn-in period
# burn_in <- ceiling(0.0 * traj_length) # arbitrary number, less than `traj_length`
# 
# # initialize accepted, rejected counters, just to monitor performance:
# n_accepted <- 0
# n_rejected <- 0
# 
# my_metropolis <- function(proposal_sd) {
#   
#   # now generate the random walk. the 't' index is time or trial in the walk.
#   # specify seed to reproduce same random walk
#   set.seed(47405)
#   
#   
#   ## I'm taking this section out and will replace it
#   
#   # # specify standard deviation of proposal distribution
#   # proposal_sd <- c(0.02, 0.2, 2.0)[2]
#   
#   ## end of the section I took out
#   
#   
#   for (t in 1:(traj_length - 1)) {
#     current_position <- trajectory[t]
#     # use the proposal distribution to generate a proposed jump
#     proposed_jump <- rnorm(1, mean = 0, sd = proposal_sd)
#     # compute the probability of accepting the proposed jump
#     prob_accept <- min(1,
#                        target_rel_prob(current_position + proposed_jump, my_data)
#                        / target_rel_prob(current_position, my_data))
#     # generate a random uniform value from the interval [0, 1] to
#     # decide whether or not to accept the proposed jump
#     if (runif(1) < prob_accept) {
#       # accept the proposed jump
#       trajectory[t + 1] <- current_position + proposed_jump
#       # increment the accepted counter, just to monitor performance
#       if (t > burn_in) {n_accepted <- n_accepted + 1}
#     } else {
#       # reject the proposed jump, stay at current position
#       trajectory[t + 1] <- current_position
#       # increment the rejected counter, just to monitor performance
#       if (t > burn_in) {n_rejected <- n_rejected + 1}
#     }
#   }
#   
#   # extract the post-burn_in portion of the trajectory
#   accepted_traj <- trajectory[(burn_in + 1) : length(trajectory)]
#   
#   tibble(accepted_traj = accepted_traj,
#          n_accepted    = n_accepted, 
#          n_rejected    = n_rejected)
#   # end of Metropolis algorithm
#   
# }