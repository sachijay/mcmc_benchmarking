# Install missing packages from CRAN
cran_packages <- c("here",
                   "mcmcse",
                   "ggplot2",
                   "latex2exp")

if (length(missing_pkgs <- setdiff(cran_packages, rownames(installed.packages()))) > 0) {
  message("Installing missing package(s): ", 
          paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

rm(cran_packages, 
   missing_pkgs)