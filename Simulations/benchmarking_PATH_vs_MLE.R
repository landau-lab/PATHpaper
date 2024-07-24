# PATH vs. MLE simulations

source("PATHpaper_functions.R")

# Benchmarking PATH vs. MLE (varying cell states)
simPATHvsMLEvarn <- simulate_PATH_MLE_grid(N_set = 1000, 
					   n_set = seq(2, 8, 2),
					   rho_set = 10^-2)
# Benchmarking PATH vs. MLE (varying cell number)
simPATHvsMLEvarN <- simulate_PATH_MLE_grid(N_set = c(100, 500, 1000), 
					   n_set = 4,
                                           rho_set = 10^-2)

