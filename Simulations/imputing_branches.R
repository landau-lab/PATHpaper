# Applying PATH to reconstructed phylogenies
# with imputed and measured branch lengths

source("PATHpaper_functions.R")

# Simulating inference and phylogenetic reconstruction
Imputed_vs_measured_recon <- simulate_PATH_recon_sets(reps = 100, N_set = 1000, 
                                                      n_set = 3, rho_set = 10^-2, 
                                                      ncuts_set = c(5, 50, 100,
                                                                    500, 1000), 
                                                      q0 = 0.5, b0 = 1, d0 = 0, 
                                                      cutrate = 0.001,
                                                      cores = 5)

# Imputing branch lengths for range of parameters
pendant_heat_df <- expand.grid("gamma" = seq(0.1, 1, 0.02), 
                               "xi" = 10^-seq(0.01, 8, 0.01)) %>% 
                     as_tibble() %>% mutate("dist" = pend_len(gamma, xi)*2)

