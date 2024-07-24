# GBM analysis

source("../PATHpaper_functions.R")

# MGH105 location heritability
MGH105_tree <- phytools::read.newick("MGH105_tree.nwk")
MGH105_location <- regmatches(MGH105_tree$tip.label, 
			      regexpr('MGH105([A-D])', 
				      MGH105_tree$tip.label)) %>% catMat()
MGH105_weights <- inv.tree.dist(MGH105_tree, node = TRUE, norm = FALSE)

# Permute cell locations to estimate phylogenetic auto-correlation sig.
MGH105_perm_df <- permute_trees(reps = 10^5, z = MGH105_location,
                                w = MGH105_weights)
MGH105_obs_df <- tibble(melt(xcor(MGH105_location, MGH105_weights)$Moran, 
                             value.name = "obs")) %>% filter(Var1 == Var2)

# GSEA
GL <- Reduce("rbind", 
             pbmclapply(MGH115_trees, 
                        function(t) gene_autocor(tree = t, rna_mat = genes115))) %>% 
       group_by(gene) %>% 
       summarise_at(vars(Z, I), mean, na.rm = TRUE) %>% 
       as_tibble() 

herit_gene_vec115 <- GL %>% arrange(desc(Z)) %$% set_names(Z, gene)

genesets <- msigdbr::msigdbr(species = "Homo sapiens", subcategory = "CGP")
pathways <- split(x = genesets$gene_symbol, f = genesets$gs_name)
pathways <- c(lapply(pathways, unique), as.list(neftel.genes))

gsea115out <- fgsea(pathways, herit_gene_vec115, maxSize = Inf, minSize = 20, 
                    scoreType = "std", nPermSimple = 10^5)
# ORA
gn <- GL %>%
  arrange(desc(Z)) %>%
  mutate(rank = rank(-Z)) %>%
  filter(rank <= 100 & Z >= 2) %$% gene

xcg_df <- lapply(MGH115_trees, 
       function(a) gene_xcor(a, genes115[,which(colnames(genes115) %in% gn)])) %>%
  melt %>% as_tibble %>% rename(gene1 = Var1, gene2 = Var2, Z = value, rep = L1) %>%
  group_by(gene1, gene2) %>%
  summarize(Z=mean(Z, na.rm=T))

h1 <- xcg_df %>% pivot_wider(id_cols = gene1, names_from = gene2, values_from = Z) %>% 
  ungroup() %>% select(-gene1) %>% as.matrix() %>%
  (function(dx) (max(dx)-dx) ) %>% 
  as.dist(.) %>% 
  hclust(method = "ward.D2")
h2 <- h1 %>% stats::cutree(k = 2)

xca_list <- lapply(1:2, function(x) {
  D <- fora(pathways, names(which(h2 == x)), 
            universe = GL$gene, 
            minSize = 20, maxSize = Inf);
  D$xcor.cluster <- x
  return(D)
} )

# Gene module (NPC-like, OPC-like, AC-like, MES-like) phylogenetic correlations
MGH115_xcor_df <- melt(lapply(MGH115_trees,
                              function(t) mod_xcor(t, modules115))) %>% tibble()

# Cell state transition inference
MGH115_PATH_inf <- pbmclapply(MGH115_trees, function(x) GBM_PATH(x, modules115))
MGH115_PATHpro_inf <- pbmclapply(MGH115_trees, function(k) GBM_PATHpro(k, modules115))

