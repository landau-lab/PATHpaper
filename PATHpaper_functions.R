# PATH paper functions

library(PATH)
library(expm)
library(Matrix)
library(pbmcapply)
library(magrittr)
library(dplyr)
library(tidyverse)
library(tidytree)
library(treeio)
library(castor)
library(reshape2)
library(fgsea)
library(phytools)
library(ape)
library(phangorn)
library(msigdbr)
library(diversitree)
library(Seurat)
library(SeuratWrappers)
library(CGE)
library(stats)


# Square matrix
sq <- function(...) {
  x <- c(...)
  t(matrix(x, sqrt(length(x)), sqrt(length(x))))
}

select <- dplyr::select

# Symmetric square matrix
sym <- function(...) {
  x <- c(...)
  l <- length(x)
  n <- (-1+sqrt(1+8*l))/2
  M <- t(matrix(0, n, n ))
  M[lower.tri(M, diag = TRUE)] <- x
  M <- t(M)
  M2 <- t(M); diag(M2) <- 0
  M <- M + M2
  return(M)
}

makeQ <- function(parameters, dimensions) {
  
  parameters <- abs(parameters)
  mat <- matrix(NA, dimensions, dimensions)
  seq1 <- seq(1, (dimensions-1)*dimensions, dimensions-1)
  seq2 <- seq(dimensions-1, (dimensions-1)*dimensions, dimensions-1)
  
  for(i in 1:dimensions) {
    mat[i,-i] <- parameters[seq1[i]:seq2[i]]
    mat[i,i] <- -sum(parameters[seq1[i]:seq2[i]])
  }
  
  return(mat)
}

# Optimize PATH inf
op_Pinf <- function(n, f, p) {
  loss_func <- function(pars, dims, freq, pat_dist) {
    Qguess <- makeQ(pars, dims)
    D <- diag(rowSums(freq))
    Pguess <- expm(Qguess*(pat_dist/2))
    X <- t(Pguess)%*%D%*%Pguess
    out <- euc(freq, X)
    return(list("out"=out, "pars"=pars))
  }
  op_search <- nlminb(runif(n^2, 0, 1), 
                      function(guess_pars) loss_func(pars=guess_pars, 
                                                     dims=n,
                                                     freq=f, pat_dist=p)$out, 
                      control = list(abs.tol=1e-2, rel.tol=1e-2, x.tol=1e-2))
  inf <- expm(makeQ(op_search$par, n))
  return(inf)
}

# PATH inf
P_inf.zw <- function(Z, W, mspd) {
  
  W <- rowNorm(W)
  W <- W/sum(W)
  Freq <- as.matrix(t(Z)%*%W%*%Z)
  u <- rowSums(Freq)
  
  Pt <- diag(1/u)%*%Freq
  rownames(Pt) <- colnames(Pt)
  
  Pinf <- tryCatch(CorrectP(expm(logm(Pt, method="Eigen")/mspd)),
                   error=function(e){op_Pinf(ncol(Freq), Freq, mspd)})
  
  return(list("Pt"=Pt, "P"=Pinf, "Q"=logm(Pinf)))
}

# MLE inf
mk_MLE <- function(phy, states="cell_state", nstates, max_iterations=1000) {
  cas <- tryCatch(castor::fit_mk(phy, 
                                 Nstates= nstates,
                                 rate_model="ARD", 
                                 tip_states=phy[[states]],
                                 optim_max_iterations=max_iterations, guess_transition_matrix),
                  error=function(e) {return(matrix(NA, nstates, nstates))} )
  return(list("cas"=cas, "Q"=cas$transition_matrix, "conv"=cas$converged))
}

melt_mat <- function(x) {
  reshape2::melt(x, value.name = deparse(substitute(x)))
}

# Simulate idealized phylogeny
simulate_PATH_ideal <- function(i, q_range = 0.5, n_cells = 2^10, n_states = 3) {
  # Generate random cell state transitions
  Q <- rQ(n_states, runif(1, 0, q_range))
  P <- expm(Q)
  # Generate ideal phylogeny
  tree <- pbtree(b = 1, d = 0, n = 2^round(log2(n_cells)), type="discrete")
  tree$states <- simulate_mk_model(tree, Q=Q)$tip_states
  # Run PATH
  Z <- catMat(tree$states, num_states = n_states)
  W <- one_node_depth(tree)
  path_run <- P_inf.zw(Z, W$W, W$mean.pat)
  Pt <- path_run$Pt
  Pinf <- path_run$P
  path.euclid <- sqrt(sum((Pinf-P)^2))
  XC <- xcor(Z, normalize(rowNorm(W$W)))
  Moran <- XC$phy_cor
  MZ <- XC$Z
  minZ <- min(diag(XC$Z), na.rm = TRUE)
  # Format into data frames
  P.df <- melt_mat(P)
  Pt.df <- melt_mat(Pt)
  Pinf.df <- melt_mat(Pinf)
  Moran.df <- melt_mat(Moran)
  MZ.df <- melt_mat(MZ)
  df <- Reduce(function(x,y) full_join(x,y, by=c("Var1", "Var2")), 
               list(P.df, Pinf.df, Pt.df, Moran.df, MZ.df)) %>% 
    mutate("replicate"=i, "nstates"=n_states, "Ntips"=n_cells,
           "rq.range"=q_range, "minZ"=minZ, "euclid.PATH"=path.euclid)
  return(df)
}

# Benchmark PATH and MLE on phylogenies from simulated somatic evolution
simulate_PATH_vs_MLE <- function(i, q_range, n_cells, n_states, sample_frac, 
                                 birth_rate = 1, death_rate = 0, runMLE = TRUE) {
  # Generate random cell state transitions
  Q <- rQ(n_states, runif(1, 0, q_range))
  P <- expm(Q)
  # Simulate somatic evolution
  tree <- PATH::ran_tree_with_states(Q0 = Q, N = n_cells, rho = sample_frac,
                                     lambda = birth_rate, mu = death_rate)
  # Run PATH
  s1 <- Sys.time()
  Z <- catMat(tree$states, num_states = n_states)
  W <- one_node_depth(tree)
  path_run <- P_inf.zw(Z, W$W, W$mean.pat)
  Pt <- path_run$Pt
  Pinf <- path_run$P
  s2 <- Sys.time()
  time1 <- as.numeric(difftime(s2, s1, units="secs"))
  path.euclid <- sqrt(sum((Pinf-P)^2))
  XC <- xcor(Z, normalize(rowNorm(W$W)))
  Moran <- XC$phy_cor
  MZ <- XC$Z
  minZ <- min(diag(XC$Z), na.rm = TRUE)
  # Format into data frames
  P.df <- melt_mat(P)
  Pt.df <- melt_mat(Pt)
  Pinf.df <- melt_mat(Pinf)
  Moran.df <- melt_mat(Moran)
  MZ.df <- melt_mat(MZ)
  df <- Reduce(function(x,y) full_join(x,y, by=c("Var1", "Var2")), 
               list(P.df, Pinf.df, Pt.df, Moran.df, MZ.df)) %>% 
    mutate("replicate"=i, "nstates"=n_states, "Ntips"=n_cells,
           "rho"=sample_frac, "rq.range"=q_range, "minZ"=minZ, 
           "euclid.PATH"=path.euclid, "time.PATH"=time1)
  # Run MLE 
  if(runMLE==TRUE) {
    s3 <- Sys.time()
    MLErun <- mk_MLE(phy = tree, states="states",
                     nstates = n_states, max_iterations = 10000)
    s4 <- Sys.time()
    time2 <- as.numeric(difftime(s4, s3, units = "secs"))
    PinfMLE <- expm(MLErun$Q)
    # Format
    PinfMLE.df <- melt_mat(PinfMLE)
    MLEconv <- MLErun$conv
    mle.euclid <- sqrt(sum((PinfMLE - P)^2))
    df <- full_join(df, PinfMLE.df, by=c("Var1", "Var2")) %>% 
      mutate("euclid.MLE" = mle.euclid,
             "time.MLE" = time2, "converge.MLE" = MLEconv)
  } 
  return(df)
}

simulate_PATH_MLE_grid <- function(reps = 1000, N_set, n_set, rho_set, q0 = 0.5, 
                                   b0 = 1, d0 = 0, runMLE = TRUE, cores = 8) {
  pars <- expand.grid("rep"=1:reps, "N"=N_set, "n"=n_set, "rho"=rho_set)
  out <- pbmclapply(1:nrow(pars),
                    function(j) simulate_PATH_vs_MLE(i = pars[j, "rep"],
                                                    q_range = q0,
                                                    n_cells = pars[j, "N"],
                                                    n_states = pars[j, "n"], 
                                                    sample_frac = pars[j, "rho"], 
                                                    birth_rate = b0, 
                                                    death_rate = d0, 
                                                    runMLE = runMLE), 
                    mc.cores = cores,
                    mc.preschedule = F)
  
  out <- Reduce(rbind, out) %>% 
    group_by(replicate, rho, nstates, Ntips) %>% 
    mutate(trial=cur_group_id())
  return(out)
}


eqQ <- function(r, n) {
  m <- matrix(r/(n-1), n, n)
  diag(m) <- -r
  return(m)
}

seq.sim <- function(tree, nsites, rate, n=2) {
  seq <- simulate_mk_model(tree, Q = eqQ(rate, n), Nsimulations = nsites, 
                           include_nodes = F)
  seq.mat <- matrix(seq$tip_states, length(tree$tip.label), nsites, byrow = TRUE)
  rownames(seq.mat) <- factor(tree$tip.label)
  return(seq.mat)
}

seq.hamm <- function(cuts) {
  proxy::dist(cuts, method = function(x,y) sum(x!=y), diag = T, upper = T)
}

# Pendant branch length est for birth-death process
pend_len <- function(gamma, xi) {
  # gamma := (birth - death), or net growth rate
  # xi := birth * sampling
  
  (xi - gamma + gamma * log(gamma/xi)) / (gamma - xi)^2
}

# Reconstruct phylogeny with UPGMA
recon_tree <- function(phy, cuts, r0) {
  names(phy$states) <- phy$tip.label
  hamm <- seq.hamm(cuts)
  times_from_hamm <- suppressWarnings(-0.5*log(1-2*(hamm/ncol(cuts)))/r0)
  times_from_hamm[which(is.nan(times_from_hamm) | is.infinite(times_from_hamm) | is.na(times_from_hamm), 
                        arr.ind = T)] <- max(times_from_hamm < Inf, na.rm = T)
  recon.phy <- phangorn::upgma(times_from_hamm)
  recon.phy$states <- phy$states[match(recon.phy$tip.label, phy$tip.label)]
  return(recon.phy)
}

# Simulate phylogenetic reconstruction and PATH inference
simulate_PATH_recon <- function(i, q_range, n_cells, n_states, sample_frac, 
                                birth_rate = 1, death_rate = 0, ncuts, r0=0.01) {
  # Generate random cell state transitions
  Q <- rQ(n_states, runif(1, 0, q_range))
  P <- expm(Q)
  # Simulate somatic evolution
  tree <- PATH::ran_tree_with_states(Q0 = Q, N = n_cells, rho = sample_frac,
                                     lambda = birth_rate, mu = death_rate)
  # Simulate tree reconstruction
  cuts <- seq.sim(tree = tree, nsites = ncuts, rate = r0, n = 2)
  recon.tree <- recon_tree(tree, cuts, r0)
  # Run PATH true
  Ztrue <- catMat(tree$states, num_states = n_states)
  Wtrue <- one_node_depth(tree)
  path_run_true <- P_inf.zw(Ztrue, Wtrue$W, Wtrue$mean.pat)
  path.true.euclid <- euc(P, path_run_true$P) 
  XCtrue <- xcor(Ztrue, rowNorm(Wtrue$W))
  # Run PATH recon
  Zrecon <- catMat(recon.tree$states, num_states = n_states)
  Wrecon <- one_node_depth(recon.tree)
  path_run_recon <- P_inf.zw(Zrecon, Wrecon$W, Wrecon$mean.pat)
  path.recon.euclid <- euc(P, path_run_recon$P) 
  XCrecon <- xcor(Zrecon, rowNorm(Wrecon$W))
  # Run PATH impute
  rhoestL <- ifelse(log10(1/sample_frac) - 1 < 0, 0, 
                    log10(1/sample_frac) - 1)
  rhoestU <- log10(1/sample_frac) + 1
  sample_frac_est <- 10^(-runif(1, rhoestL, rhoestU))
  imputed_pat <- pend_len(birth_rate - death_rate,
                          birth_rate * sample_frac_est)*2 - 1
  path_run_imp <- P_inf.zw(Zrecon, Wrecon$W, imputed_pat)
  path.imp.euclid <- euc(P, path_run_imp$P) 
  #
  MIdist <- euc(XCtrue$phy_cor, XCrecon$phy_cor)
  RFdist <- tree_distance(tree, recon.tree, normalized = TRUE, 
                          metric = "RFrooted")
  MPLdist <- tree_distance(tree, recon.tree, normalized = TRUE, 
                           metric = "MeanPathLengthDifference")
  # Format into data frames
  df <- data.frame("replicate" = i, "nstates" = n_states, "Ntips" = n_cells, 
                   "rho" = sample_frac, "ncuts" = ncuts, "rq.range" = q_range,
                   "euclid.PATHtrue" = path.true.euclid,
                   "euclid.PATHrecon" = path.recon.euclid,
                   "euclid.PATHimp" = path.imp.euclid,
                   "dist.MI" = MIdist, "dist.RF" = RFdist, 
                   "dist.MPL" = MPLdist)
  return(df)
}


simulate_PATH_recon_sets <- function(reps, N_set, n_set, 
                                     rho_set, ncuts_set,
                                     q0, b0, d0,
                                     cutrate, cores) {
  pars <- expand.grid("rep"=1:reps, "N"=N_set, "n"=n_set, "rho"=rho_set, "ncuts"=ncuts_set)
  out <- pbmclapply(1:nrow(pars),
               function(j) simulate_PATH_recon(i = pars[j, "rep"], q_range = q0,
                                               n_cells = pars[j, "N"],
                                               n_states = pars[j, "n"], 
                                               sample_frac = pars[j, "rho"], 
                                               birth_rate = b0, death_rate = d0, 
                                               ncuts = pars[j, "ncuts"], 
                                               r0 = cutrate), 
               mc.cores = cores,
               mc.preschedule = F)
  
  out <- Reduce(rbind, out) %>% 
    group_by(replicate, rho, nstates, Ntips, ncuts) %>% 
    mutate(trial=cur_group_id())
  return(out)
}

# Process EMT PDAC clonal phylogenies
cleanup_EMT_tree <- function(tree_data) {
  N <- length(tree_data@phylo$tip.label)
  na <- which(is.na(tree_data@data$pseudotime[1:N]))
  tree_data@data <- tree_data@data[1:N,]
  tree_data@data <- tree_data@data[-na,]
  tree_data@phylo <- drop.tip(tree_data@phylo, na)
  return(tree_data)
}

subsample_zw <- function(state_mat, weight_mat, subsample_size) {
  s <- sample(1:nrow(state_mat), subsample_size, replace = F)
  zn <- state_mat[s,]
  wn <- weight_mat[s,s]
  return(list("z"=zn, "w"=wn))
}

subsample_xcor <- function(state_mat, weight_mat, subsample_size) {
  sub <- subsample_zw(state_mat, weight_mat, subsample_size)
  xcor(sub$z, rowNorm(sub$w))$Z
}

# Resample phylogenies to get CIs for xcor
bootstrap_xcor <- function(reps = 100, state_mat, weight_mat, subsample_size) {
  melt(pbmclapply(1:reps, function(x) subsample_xcor(state_mat, weight_mat, subsample_size)))
}

# Est branch lengths for EMT Pro inf
est_inf_pars <- function(Nsamp, birth, death = 0, age, cycle_length) {
  Npop <- 1*exp((birth-death)*age/cycle_length)
  rhoest <- Nsamp/Npop
  stopifnot(rhoest < 1)
  mean(replicate(10, one_node_depth(PATH:::ran_tree(N=Nsamp, rho = rhoest, lambda = birth, mu = death))$mean.pat))
}

# Gene heritability
gene_autocor <- function(tree, rna_mat) {
  gene_names <- colnames(rna_mat)
  genes <- full_join(tree, as.data.frame(rna_mat) %>% 
                       rownames_to_column(var = "label")) %>% 
    filter(isTip==T) %>% get_tree_data() %>% 
    select(gene_names) %>% as.matrix()
  
  w <- inv_tree_dist(tree, node = TRUE, norm = FALSE)
  
  PATH:::parM(genes, w) %>% rownames_to_column(var = "gene")
}

# Gene xcor
gene_xcor <- function(tree, rna_mat) {
    gene_names <- colnames(rna_mat)
    genes <- full_join(tree, as.data.frame(rna_mat) %>% 
                           rownames_to_column(var = "label")) %>% 
        filter(isTip==T) %>% get_tree_data() %>% 
        select(gene_names) %>% as.matrix()
    
    w <- inv_tree_dist(tree, node = TRUE, norm = FALSE)
    
    melt(xcor(genes, w)$Z) %>% as_tibble()
}

# Gene module xcor
mod_xcor <- function(tree, mod, weight_function = function(t) exp_tree_dist(t)) {
  z <- full_join(tree, as.data.frame(mod) %>% 
              rownames_to_column(var = "label"), by = "label") %>% 
      filter(isTip==TRUE) %>% get_tree_data() %>% 
      select("NPC", "OPC", "AC", "MES") %>% as.matrix

  w <- weight_function(tree)
  xcor(z, w)$Z
}

# Apply PATH to GBM
GBM_PATH <- function(tree, mod) {
  tree <- castor::extend_tree_to_height(tree)$tree
  mp <- one_node_depth(tree)$mean.pat
  tree$edge.length <- tree$edge.length/mp * 2
  z <- full_join(tree, as.data.frame(mod) %>% 
                   rownames_to_column(var = "label"), by = "label") %>% 
    mutate(state = c("Stem", "Stem", "AC", "MES")[
      max.col(select(., c("NPC", "OPC", "AC", "MES")))]) %>% 
    filter(isTip==T) %>% get_tree_data() %>% select(state) %>% as.matrix %>%
    catMat(., state_order = c("Stem", "AC", "MES")) 
  W0 <- one_node_depth(tree)
  p <- P_inf.zw(z, W0$W, W0$mean.pat)
  return(p)
}

# Apply PATHpro to GBM
GBM_PATHpro <- function(tree, mod, depth = 1) {
  tree <- castor::extend_tree_to_height(tree)$tree	
  mp <- one_node_depth(tree)$mean.pat
  tree$edge.length <- tree$edge.length/mp * 2
  z <- full_join(tree, as.data.frame(mod) %>% 
                   rownames_to_column(var = "label"), by = "label") %>% 
    mutate(state = c("Stem", "Stem", "AC", "MES")[
      max.col(select(., c("NPC", "OPC", "AC", "MES")))]) %>% 
    filter(isTip==T) %>% get_tree_data() %>% select(state) %>% as.matrix
  p <- PATHpro(tree = tree, cell_states = z, depth = depth, cell_state_order = c("Stem", "AC", "MES"))
  return(p)
}

# PATH vs. PATHpro on forward time sims
PATH_vs_PATHPro_sim <- function(prolif = c(1,1), Q, N = 10000, n = 1000, death = c(0,0)) {
  tree <- tree.musse(pars = c(prolif, death, convertQ(Q)), max.taxa = N, 
                     x0 = sample(length(prolif), 1, 
                                 prob = PF_eig(t(diag(prolif) + Q))$vec))
  if(is.null(tree)) {
    return(NULL)
  }
  
  s <- sample(1:N, N-n)
  tree <- ape::drop.tip(tree, s)
  tree$tip.state <- tree$tip.state[tree$tip.label]
  z <- catMat(tree$tip.state)
  w0 <- one_node_depth(tree)
  w <- normalize(rowNorm(w0$W))
  mspd <- w0$mean.pat
  
  Pinf <- P_inf.zw(z, w, mspd)$P
  
  v1 <- PF_eig(t(diag(prolif-death) + Q))$val
  Proinf <- PATHpro(tree, tree$tip.state, 3, total_prolif = v1, alpha = NULL,
		    guess_list = list(runif(ncol(z)^2, 0, 0.1)))$P
  Ptrue <- expm(Q)
  PATHEuclid <- euc(Ptrue, Pinf)
  ProEuclid <- euc(Ptrue, Proinf)
  
  fracsdf <- data.frame(t(colMeans(z)))
  colnames(fracsdf) <- gsub("X", "state", colnames(fracsdf))
  prolifdf <- data.frame(t(prolif))
  colnames(prolifdf) <- gsub("X", "prolif", colnames(prolifdf))
  tb <- tibble("PATHEuclid"=PATHEuclid, "ProEuclid"=ProEuclid,
               prolifdf, fracsdf)
  return(tb)
}

# Retrieve forward-equivalent (FE) phylogeny sim
get_FEtree <- function(dir, tree_num = 1) {
  filepath <- paste0(dir, "/trees/")
  l <- list.files(filepath)
  stopifnot(length(l) > 0)
  tree_data <- treeio::read.nhx(paste0(filepath, l[tree_num]))
  tree_data@phylo <- collapse_monofurcations(tree_data@phylo, 
                                             force_keep_root = F)$tree
  tree <- tree_data@phylo
  w0 <- one_node_depth(tree)
  w <- w0$W
  pat <- w0$mean.pat
  states <- tree_data %>% filter(isTip == T) %>% get_tree_data() %$% state
  z <- catMat(states)
  u <- colMeans(z)
  out <- list("tree"=tree, "z"=z, "w"=w, "pat"=pat, "states"=states,
              "u"=u, "data"=tree_data)
  return(out)
}

# Benchmark PATHpro and SSE-MLE on FE trees
SSE_vs_PATHpro_FE <- function(x, P, birth, death, total_prolif = NULL, depth = 3, 
                              rho = NULL) {
  n <- ncol(P)
  N <- length(x$tree$tip.label)
  true_prolif <- birth - death
  s1 <- Sys.time()
  Pro <- PATHpro(x$tree, x$states, depth, total_prolif = total_prolif, alpha = depth/(n^2 - n))
  s2 <- Sys.time()
  Protime <- as.numeric(difftime(s2, s1, units = "secs"))
  Pro_eucP <- euc(P, Pro$P)
  Pro_eucG <- euc(true_prolif, Pro$gamma)
  
  s3 <- Sys.time()
  SSE <- fit_musse(x$tree, n , tip_pstates = x$states, sampling_fractions = rho, 
                   transition_rate_model = "ARD")
  s4 <- Sys.time()
  Psse <- expm(SSE$parameters$transition_matrix)
  gsse <- SSE$parameters$birth - SSE$parameters$death
  SSEtime <- as.numeric(difftime(s4, s3, units = "secs"))
  SSE_eucP <- euc(P, Psse)
  SSE_eucG <- euc(true_prolif, gsse)
  
  out <- c("Pro_eucP" = Pro_eucP, "SSE_eucP" = SSE_eucP,
           "Pro_eucG" = Pro_eucG, "SSE_eucG" = SSE_eucG,
           "Protime" = Protime, "SSEtime" = SSEtime)
  return(list("res" = out, 
              "PATH_P" = Pro$P, "SSE_P" = Psse,
              "PATH_g" = Pro$gamma, "SSE_g" = gsse))
}

SSE_vs_PATHpro_FE_dir <- function(dir, Pmat, births, deaths, rho = NULL, 
                                  depth = 3, ntrees = 100, cores = NULL) {
  if(is.null(cores)) {
    cores <- detectCores()
  }
  v <- CGE::PF_eig(t(diag(births - deaths) + logm(Pmat)))$val
  infs <- pbmclapply(1:ntrees,
                     function(i) SSE_vs_PATHpro_FE(x = get_FEtree(dir, i), 
                                                   P = Pmat, 
                                                   birth = births, 
                                                   death = deaths, 
                                                   total_prolif = v, 
                                                   depth = depth, rho = rho), 
                     mc.cores = cores, mc.preschedule = F)
  df <- as_tibble(t(sapply(infs, '[[', "res")))
  PATH_P <- lapply(infs, '[[', "PATH_P")
  PATH_prolif <- lapply(infs, '[[', "PATH_g")
  SSE_P <- lapply(infs, '[[', "SSE_P")
  SSE_prolif <- lapply(infs, '[[', "SSE_g")
  out <- list("df"=df, "PATH_P" = PATH_P, "SSE_P" = SSE_P, 
              "PATH_prolif" = PATH_prolif, "SSE_prolif" = SSE_prolif, 
              "True_P" = Pmat, "True_prolif" = births - deaths, "rho" = rho,
              "depth" = depth, "dir" = dir, "infs" = infs)
  return(out)
}

# List mean
matMean <- function(M, func = mean) {
    apply(simplify2array(M), 1:2, func, na.rm = TRUE)
}

# PATHpro with subsampling replicates
PATHpro.resample <- function(Z, W, mspd, 
                             sample_size = NULL, reps = 100, 
                             pop_prolif_est = NULL, 
                             state_prolif_ests = NULL, alph, cores = 8) {
  if(is.null(sample_size)) {
    sample_size <- round(0.5*nrow(Z))
  }
  
  subsamplePro <- function(Z, W, mspd, s, pop_pro, state_pros) {
    s <- sample(1:nrow(Z), s)
    Zs <- Z[s,]
    Ws <- W[s,s]
    n <- ncol(Z)
    PATHpro.XW(X = Zs, W = Ws, t = mspd, 
               total_prolif_rate_est = pop_pro,
	       alpha0 = alph,
               cell_prolifs = NULL)
  }
  
  L <- pbmclapply(1:reps, function(x) subsamplePro(Z, W, mspd, 
                                                   sample_size, 
                                                   pop_prolif_est,
                                                   state_prolif_ests), 
                  mc.cores = cores)
  
  Pmean <- matMean(lapply(L, '[[', "P"))
  Pmedian <- matMean(lapply(L, '[[', "P"), median)
  Pvar <- matMean(lapply(L, '[[', "P"), var)
  Qmean <- matMean(lapply(L, '[[', "Q"))
  Qvar <- matMean(lapply(L, '[[', "Q"), var)
  
  pro_mean <- apply(t(sapply(L, '[[', "gamma")), 2, mean)
  pro_var <- apply(t(sapply(L, '[[', "gamma")), 2, var)
  
  output <- list("P" = Pmean, "Pmed" = Pmedian, "Q" = Qmean, "prolif" = pro_mean,
                 "Pvar" = Pvar, "Qvar" = Qvar, "prolifvar" = pro_var, "sims" = L)
  return(output)
}
