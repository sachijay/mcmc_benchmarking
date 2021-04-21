## Installs missing packages
source(here::here("src", "install_packages.R"))

## Load libraries
library(ggplot2)


## Number of samples obtained by default
n_samples <- 10000

## Burn-in proportion
burn_prop <- 0.5

## Figures directory
fig_dir <- paste0(here::here("figures"),
                  "/")
