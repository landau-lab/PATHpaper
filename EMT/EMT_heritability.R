# EMT heritability and plasticity in a mouse model of PDAC

source("../PATHpaper_functions.R")

# Mouse 1 Clone 1 from Simeonov et al. 2021.
mouse1cl1 <- cleanup_EMT_tree(readRDS("EMTtrees/Mouse1.Clone1.extended.rds"))
EMTtree <- mouse1cl1@phylo
n_cells_pdac <- length(EMTtree$tip.label)

# weight matrix from phylogeny
w <- rowNorm(one_node_depth(EMTtree)$W)

# Bin pseudoEMT into discrete states
em <- ifelse(mouse1cl1@data$pseudotime[1:n_cells_pdac] <= 7, 
             7, ceiling(mouse1cl1@data$pseudotime[1:n_cells_pdac]))
embh <- ifelse(mouse1cl1@data$pseudotime[1:n_cells_pdac] <= 7, 
             7, ceiling(2*mouse1cl1@data$pseudotime[1:n_cells_pdac])/2)
emb2 <- ifelse(mouse1cl1@data$pseudotime[1:n_cells_pdac] <= 7, 
             7, ceiling(mouse1cl1@data$pseudotime[1:n_cells_pdac]/2)*2)
emb3 <- ifelse(mouse1cl1@data$pseudotime[1:n_cells_pdac] <= 7, 
             7, ceiling(mouse1cl1@data$pseudotime[1:n_cells_pdac]/3)*3)

# Bin into four EMT states: T1, T2, T3, and M
em2 <- data.frame("Var1"=em) %>%
  mutate(state2 = case_when(
    between(Var1, -Inf, 13)  ~ "T1",
    between(Var1, 14, 24) ~ "T2",
    between(Var1, 25, 30) ~ "T3",
    Var1 >= 31             ~ "M"
  )) %$% state2

# Organize data into data frame. 
EMTdf <- data.frame("state" = mouse1cl1@data$orig.ident[1:n_cells_pdac], 
                    "Cycle" = mouse1cl1@data$G2M.Score[1:n_cells_pdac],
                    "EMT" = mouse1cl1@data$pseudotime[1:n_cells_pdac],
                    "EMTstate" = factor(em2, ordered = TRUE, 
					levels = c("T1", "T2", "T3", "M")),
                    "label" = EMTtree$tip.label)

# Measuring heritability
# Measure phylogenetic auto-correlations of harvest site location (EMTlocI), 
# pseudotimeEMT (continuous) and G2M (EMTpseudo_cycleI), and
# binned EMT states (EMTpseudobinI).
# For each analysis, 4000 of the 7968 cells were resampled 1000 times.
samp_size <- 4000
resamp_reps <- 1000
EMTlocI <- bootstrap_xcor(resamp_reps, catMat(EMTdf$state), w, samp_size)
EMTpseudo_cycleI <- bootstrap_xcor(resamp_reps, as.matrix(EMTdf[c("EMT", "Cycle")]),
				   w, samp_size)
EMTpseudobinI <- bootstrap_xcor(resamp_reps, catMat(em), w, samp_size)

# Phylogenetic cross-correlation of EMT pseudotime bins
bin_xcor_df <- rbind(melt(xcor(catMat(factor((embh))), w)$Z) %>% 
                         mutate("Bin size" = 0.5), 
                     melt(xcor(catMat(factor((em))), w)$Z) %>% 
                         mutate("Bin size" = 1),
                     melt(xcor(catMat(factor((emb2))), w)$Z) %>% 
                         mutate("Bin size" = 2),
                     melt(xcor(catMat(factor((emb3))), w)$Z) %>% 
                         mutate("Bin size" = 3)) %>% tibble() 

# Location phylogenetic cross-correlations 
EMT_loc_xcor <- melt(xcor(catMat(EMTdf$state), w)$Z) %>% tibble

