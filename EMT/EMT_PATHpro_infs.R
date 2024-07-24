# PATHpro inference on composite EMT phylogeny

source("PATHpaper_functions.R")

# Load clonal phylogenies
supertree_files <- lapply(grep("extended",
        list.files("EMTtrees",
        pattern = "Mouse1"), value = T), 
        function(x) readRDS(paste0("EMTtrees/", x)))

# Merge phylogenies with at least 10 cells
supertree_filter <- which(sapply(supertree_files, 
                                function(x) length(x@phylo$tip.label)) < 10)
supertrees <- lapply(supertree_files[-supertree_filter], cleanup_EMT_tree)

# Cell state transition inference
emA <- Reduce("c", sapply(supertrees, '[[', "pseudotime"))
emC <- Reduce("c", sapply(supertrees, '[[', "Phase"))
wST <- Reduce(bdiag, lapply(supertrees, function(x) one_node_depth(x@phylo)$W))

STstates <- data.frame("Var1" = ceiling(emA)) %>%
  mutate(state = case_when(
    between(Var1, -Inf, 13)  ~ "T1",
    between(Var1, 14, 24) ~ "T2",
    between(Var1, 25, 30) ~ "T3",
    Var1 >= 31             ~ "T4"
  )) %$% state

nST <- length(STstates)
Zsup <- catMat(STstates)

# Infer transition probabilities with PATHpro between EMT cell states:
# T1, T2, T3, and M.

# Arguments for PATHpro include phylogeny, cell states, mean patristic distance,
# and optionally cell state proliferation rates. 

# Cell state proliferation rates are estimated from the fraction of cycling cells
# in each state using scRNAseq data.
cyST <- sapply(c("T1", "T2", "T3", "T4"), 
              function(r) length(which(STstates == r & emC %in%  c("S", "G2M"))) /
                length(which(STstates == r)))

# Estimate patristic distances with tumor proliferation rate, cell sampling rate, 
# and cell division duration. 
emt_state_birthsSC <- log(cyST + 1)
emt_total_birthsSC <- sum(emt_state_birthsSC * colMeans(Zsup))

t1A <- est_inf_pars(Nsamp = nST, birth = emt_total_birthsSC, 
                    death = 0, age = 35, cycle_length = 0.5)
t2A <- est_inf_pars(Nsamp = nST, birth = emt_total_birthsSC, 
                    death = 0.01, age = 35, cycle_length = 0.5)
t3A <- est_inf_pars(Nsamp = nST, birth = emt_total_birthsSC, 
                    death = 0.02, age = 35, cycle_length = 0.5)

# PATHpro with resampling
super1 <- PATHpro.resample(Z = Zsup, W = wST, mspd = t1A, alph = 1/12,
                           pop_prolif_est = emt_total_birthsSC - 0)
super2 <- PATHpro.resample(Z = Zsup, W = wST, mspd = t2A, alph = 1/12,
                           pop_prolif_est = emt_total_birthsSC - 0.01)
super3 <- PATHpro.resample(Z = Zsup, W = wST, mspd = t3A, alph = 1/12,
                           pop_prolif_est = emt_total_birthsSC - 0.02)

