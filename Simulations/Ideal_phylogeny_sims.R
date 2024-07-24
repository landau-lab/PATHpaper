# PATH applied to ideal phylogenies

source("PATHpaper_functions.R")

ideal_sims <- as_tibble(Reduce(rbind, pbmclapply(1:1000, simulate_PATH_ideal)))

