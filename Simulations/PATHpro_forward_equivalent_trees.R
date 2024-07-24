# Benchmark PATHpro vs. SSE-MLE on forward equivalent simulated phylogenies

source("../PATHpaper_functions.R")

# Model 1
dir1 <- "SIM1"
P1 <- sq(0.9, 0.1,
         0.1, 0.9)
g1 <- c(2, 1)
d1 <- c(0, 0)
rho1 <- 1e-2

# Model 2
dir2 <- "SIM2"
P2 = sq(0.9, 0.1, 0,
        0.05, 0.85, 0.1,
        0, 0.05, 0.95)
g2 <- c(1.5, 1.25, 1)
d2 <- c(0.1, 0.1, 0.1)
rho2 <- 1e-5

# Model 3
dir3 <- "SIM3"
P3 <- sq(0.93, 0.05, 0.01, 0.01,
         0.01, 0.89, 0.05, 0.05,
         0.01, 0.01, 0.93, 0.05,
         0.01, 0.01, 0.05, 0.93)
g3 = c(2.25, 2, 2, 2)
d3 <- c(0.25, 0.25, 0.25, 0.25)
rho3 <- 1e-4

# PATHpro and SSE-MLE inferences
sim1_pro <- SSE_vs_PATHpro_FE_dir(dir = dir1, Pmat = P1,  births = g1,
                                  deaths = d1, rho = rho1, ntrees = 100)
sim2_pro <- SSE_vs_PATHpro_FE_dir(dir = dir2, Pmat = P2, births = g2,
                                  deaths = d2, rho = rho2, ntrees = 100)
sim3_pro <- SSE_vs_PATHpro_FE_dir(dir = dir3, Pmat = P3, births = g3,
                                  deaths = d3, rho = rho3, ntrees = 100)
