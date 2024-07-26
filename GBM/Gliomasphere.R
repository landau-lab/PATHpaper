# Gliomasphere analysis

source("../PATHpaper_functions.R")

# Barcode matrix location, e.g.,
# no3_path <- "/PCT062_no03/sample_feature_bc_matrix/"
# no4_path <- "/PCT062_no04/sample_feature_bc_matrix/"
# clonal phylogenies in Newick format in "/clonetrees/"

# Code to process gliomasphere RNA and phylogeny data

# Use Seurat to process RNA seq data
library(Seurat)
process_invitro_seurat <- function(RNA_path, rep, mt.perc = 25, min.feat = 200) {
  
  wk_10X <- Read10X(RNA_path)
  wk_data <- CreateSeuratObject(wk_10X$`Gene Expression`)
  wk_data[["MultiplexingCapture"]] <- CreateAssayObject(wk_10X$`Multiplexing Capture`)
  wk_data[["percent.mt"]] <- PercentageFeatureSet(wk_data, pattern = "^MT-")
  wk_data <- subset(wk_data, subset = nFeature_RNA > min.feat & percent.mt < mt.perc) 
  wk_data <- NormalizeData(wk_data)
  wk_data <- FindVariableFeatures(wk_data)
  wk_data <- ScaleData(wk_data, features = rownames(wk_data))
  wk_data <- AddModuleScore(wk_data, features = c(neftel_genes, cc.genes.updated.2019), 
                            name = c(names(neftel_genes), "G2S", "G2M") , search = T)
  
  wk_data@misc <- data.frame("rep"=rep)
  return(wk_data)
}

# List clone files
list_clone_numbers <- function(rep, dir = "clonetrees/") {
  tr_subset <- paste0("wk4.subset.no", rep)
  sort(as.numeric(gsub(".nwk", "", 
		       gsub(paste0(tr_subset, "_"), "", 
			    grep(tr_subset, list.files(dir), value = T)))))
}

# Process trees
proc_invtro_nwk <- function(wk_data, clone, dir) {
  rep <- wk_data@misc$rep
  tr_fl <- paste0(dir, "wk4.subset.no", rep, "_", clone, ".nwk")
  tree.subset <- read.newick(tr_fl)
  mod.subset <- wk_data[[c("NPC1","OPC2","AC3","MES4", "G2S5", "G2M6")]][which(gsub("-1", "", colnames(wk_data)) %in% tree.subset$tip.label),]
  rownames(mod.subset) <- gsub("-1", "", rownames(mod.subset))
  tib.subset <- as.treedata(full_join(as_tibble(tree.subset), tibble(mod.subset, label=rownames(mod.subset))))
  tib.subset@data$rep <- rep
  tib.subset@data$clone <- clone

  tib.subset@data$state <- tib.subset@data %>% 
    select(NPC1, OPC2, AC3, MES4) %>% 
    mutate(state=max.col(.)) %>% 
    mutate(state.name = c("NPC", "OPC", "AC", "MES")[state]) %>% 
    select(state.name) %>% unlist
return(tib.subset)
}

get_annotated_clone_trees <- function(wk_data, tree_drc="clonetrees/") {
  rep <- wk_data@misc$rep
  clonetrees <- list()
  count <- 1
  for(i in list_clone_numbers(rep, tree_drc)) {
    clonetrees[[count]] <- proc_invitro_nwk(wk_data, i, tree_drc)
    count <- count + 1
  }
  return(clonetrees) 
}

# Get cells state and phylogenetic weight matrices for clone trees
clonetrees.WX <- function(clonetrees) {
	select_GBM_mods <- function(tree_data) {
		tree_data %>% filter(isTip==T) %>% 
			get_tree_data() %>% 
			select(NPC1, OPC2, AC3, MES4) %>% 
			as.matrix()
	}
    Wlist <- lapply(clonetrees, function(x) as.matrix(one_node_depth(x@phylo)$W*1))
    Xlist <- lapply(clonetrees, function(x) select_GBM_mods(x))
    W <- bdiag(Wlist)
    X <- Reduce("rbind", Xlist)
    drp <- which(apply(X, 1, function(x) any(is.na(x)))==T)
    if(length(drp) > 0) {
      W <- W[-drp,-drp]
      X <- X[-drp,]
    }
    return(list("W" = W, "X" = X))
}

# Calculate phylogenetic correlations
xcor_from_clonetrees_df <- function(wk_data, tree_direct = "clonetrees/") {
  rep <- wk_data@misc$rep
  ct <- get_annotated_clone_trees(wk_data, tree_drc = tree_direct)
  if(is.null(ct)==TRUE) {
    return()
  } else {
  w_and_x <- clonetrees.WX(ct)
  xr <- xcor(w_and_z$X, w_and_z$W)
  xrZ <- xr$Z
  xrI <- xr$phy_cor
  
  df1 <- melt(xrZ, value.name = "Z")
  df2 <- melt(xrI, value.name = "I")
  
  df <- full_join(df1, df2, by=c("Var1", "Var2"))
  df$rep <- rep
  return(df)
  }
}

sphere_df <- rbind(xcor_from_clonetrees_df(process_invitro_seurat(no3_path, 3)), 
		   xcor_from_clonetrees_df(process_invitro_seurat(no4_path, 4)))
