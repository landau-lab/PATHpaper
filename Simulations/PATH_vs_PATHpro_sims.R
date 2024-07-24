# PATH vs. PATHpro inferences on subsampled forward-time simulated phylogenies

source("PATHpaper_functions.R")

PATHvProsim <- Reduce("rbind", 
                 lapply(seq(1.2, 2, 0.2), 
                   function(g) Reduce("rbind", 
                     pbmclapply(1:100, 
                     function(t) PATH_vs_PATHPro_sims(prolif = c(g, 1),
						      Q = sym(-0.1, 0.1, -0.1), 
						      N = 10^4, n = 10^3, 
						      death = c(0, 0))))))


