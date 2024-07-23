# PATH examples

source("PATHpaper_functions.R")

# Simulate heritable, plastic, and three state phylogenies.

herit_tree <- ran_tree_with_states(logm(sq(0.9, 0.1, 0.1, 0.9)), 
                                   N = 200, rho = 1, lambda = 1, mu = 0)
plastic_tree <- ran_tree_with_states(sq(-10, 10, 10, -10), 
                                     N = 200, rho = 1, lambda = 1, mu = 0)
three_state_tree <- ran_tree_with_states(logm(sq(0.5, 0.5, 0,
                                                0.45, 0.5, 0.05,
                                                0, 0.1, 0.9)), N = 200, rho = 1,
                                        lambda = 1, mu = 0)



