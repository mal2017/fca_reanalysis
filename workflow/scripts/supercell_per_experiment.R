if (!require("SuperCell")) {
  remotes::install_github("GfellerLab/SuperCell",ref="v1.0")
}

library(SuperCell)
library(tidyverse)
library(scater)
library(scran)
library(remotes)

# See https://github.com/GfellerLab/SuperCell/blob/master/vignettes/a_SuperCell.Rmd
# for workflow I recreate here

sce_fl <- "results/import/fca_annotated_r_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5.usa.fca_annot.sce.rds"
sce_fl <- snakemake@input$sce
sce <- read_rds(sce_fl)

sce <- logNormCounts(sce)

# from osca:
# Our recommendation is to simply pick an arbitrary n 
# and proceed with the rest of the analysis, with the intention
# of testing other choices later, rather than spending much time 
# worrying about obtaining the “optimal” value.
dec <- modelGeneVar(sce)
hvg <- getTopHVGs(dec, n=1000)
hvg <- hvg[str_detect(hvg,"FBgn")]

gamma <- 20 # graining level
gamma <- snakemake@params$gamma

k.knn <- 5
k.knn <- snakemake@params$k_knn

SC <- SCimplify(logcounts(sce),  # gene expression matrix 
                #cell.split.condition=sces$time,
                k.knn = k.knn, # number of nearest neighbors to build kNN network
                gamma = gamma, # graining level
                genes.use = hvg,
)

# by default, this assumes input counts are normalized and log transformed
SC.GE <- supercell_GE(logcounts(sce), SC$membership)

SC$lineage <- supercell_assign(clusters = sce$annotation, # single-cell assigment to cell lines (clusters)
                               supercell_membership = SC$membership, # single-cell assignment to metacells
                               method = "jaccard")

#dimensionality reduction 
SC.PCA  <- supercell_prcomp(Matrix::t(SC.GE), # metacell gene expression matrix
                            genes.use = SC$genes.use, # genes used for the coarse-graining, but any set can be provided
                            supercell_size = SC$supercell_size, # sample-weighted pca
                            k = 20) 
## compute distance
D <- dist(SC.PCA$x)

## cluster metacells
SC.clusters    <- supercell_cluster(D = D, k = length(unique(sce$annotation)), supercell_size = SC$supercell_size) 
SC$clustering  <- SC.clusters$clustering

## mapping metacell cluster to cell line 
map.cluster.to.lineage    <- supercell_assign(supercell_membership = SC$clustering, clusters  = SC$lineage)

## clustering as cell line
SC$clustering_reordered   <- map.cluster.to.lineage[SC$clustering]

res <- list(SC = SC, SC.GE=SC.GE, SC.PCA=SC.PCA, k.knn=k.knn, gamma=gamma)

scsce <- supercell_2_sce(SC.GE, SC, fields=c("lineage"))

metadata(scsce)$superccell_k.knn <- k.knn
metadata(scsce)$supercell_gamma <- gamma

write_rds(res,snakemake@output$res_l)
write_rds(scsce,snakemake@output$sce)

# ------------------------------------------------------------------------------
# make scatterplot of 1 feature vs another
# plot_indiv_relationship <- function(x, y, caption="") {
#   supercell_GeneGenePlot(SC.GE,gene_x = x, gene_y = y,supercell_size = SC$supercell_size,color.use = "black") -> res
#   
#   g <- res$p$data |> 
#     as_tibble() |>
#     ggplot(aes(x, y, size=size)) +
#     ggdensity::geom_hdr_points(method="mvnorm") +
#     xlab(x) +
#     ylab(y) +
#     labs(title=sprintf("%s vs %s",x,y),subtitle = sprintf("raw weighted corr:%s; p=%s",round(res$w.cor[[1]],digits = 3),format.pval(res$w.pval[[1]],digits = 3)),caption = caption) +
#     geom_smooth(method="lm",se=F,color="red",linetype="dashed",linewidth=1) +
#     theme(aspect.ratio = 1)
#   
#   return(g)
# }
# 
# # pan FBgn0085432
# # toll-7 FBgn0034476
# # nkd FBgn0002945
# # egg FBgn0086908
# # piwi FBgn0004872
# # mael FBgn0016034
# # cg4115 FBgn0038017
# # nplp2 FBgn0287423
# library(patchwork)
# (plot_indiv_relationship("FBgn0085432","FBgn0034476") + labs(title="pan vs toll-7 (control, positively regulated by pan)")) +
# (plot_indiv_relationship("FBgn0085432","FBgn0264753") + labs(title="pan vs Rgk1 (control, positively regulated by pan)")) +
# (plot_indiv_relationship("FBgn0085432","1360") + labs(title="pan vs 1360 (TE also highly correlated in embryo)")) +
# (plot_indiv_relationship("FBgn0287423","FBgn0085432") + labs(title="pan vs Nplp2 (control, negatively regulated by pan)") ) + plot_layout(guides="collect") & theme(title = element_text(size=7))
# 
# 
# # ------------------------------------------------------------------------------
# 
# tfs <- read_tsv("~/work/tetf3/resources/Drosophila_melanogaster_TF.txt")
# tfs <- unique(tfs$Ensembl)
# 
# tes <- jsonlite::read_json("~/work/tetf3/upstream/te_element_lookup.json") %>%
#   names()
# 
# get_weighted_cor <- function(x, n) {
#   cm <- t(x) |>
#     as.matrix() |>
#     psych::cor.wt(w = n)
#   
#   cm$r
# }
# 
# # fun to get spqn-corrected corr coefs
# get_corrected_cor <- function(x, ave_GE) {
#   ngrp <- 20; sizegrp <- 200; refgrp <- 6
#   
#   cormat.spqn <- spqn::normalize_correlation(x, 
#                                              ave_exp = ave_GE,
#                                              ngrp = ngrp,
#                                              size_grp = sizegrp,
#                                              ref_grp = refgrp)
#   
#   rownames(cormat.spqn) <- rownames(x)
#   colnames(cormat.spqn) <- colnames(x)
#   cormat.spqn
# }
# 
# 
# get_p_from_r <- function(r, n) {
#   t_stat <- r *sqrt(n-2) / sqrt(1-r^2)
#   # ripped from source code of stats::cor.test
#   p <- 2 * min(pt(t_stat, n-2), pt(t_stat, n-2, lower.tail=FALSE))
#   p
# }
# 
# # annotates and reshapes correlatio coef mat
# cormat2tbl <- function(x,n) {
#   x |> 
#     as_tibble(rownames = "feature") |>
#     pivot_longer(-feature,names_to = "y", values_to = "coef") |>
#     filter(feature %in% tfs & y %in% tes & feature!=y) |>
#     mutate(p=map_dbl(coef,get_p_from_r,n=n)) |>
#     mutate(padj = p.adjust(p,method="BH"))
# }
# 
# 
# 
# # raw info for calculating results
# df <-  c(list(all_cells = SC$lineage)) |>
#   map(names) |>
#   map(as.integer) |>
#   enframe(name = "lineage", value = "SC_idx") |>
#   mutate(n=map_dbl(SC_idx, length)) |>
#   mutate(SC_sizes = map(SC_idx, ~{SC$supercell_size[.x]})) |>
#   mutate(GE_mat = map(SC_idx, ~{SC.GE[rownames(SC.GE) %in% c(tfs,tes),.x]}))
# 
# df <- filter(df, lineage %in% c("all_cells"))
# 
# # get weighted correlation
# df <- mutate(df, cor.raw = map2(GE_mat, SC_sizes, get_weighted_cor))
# 
# # get average expression used for spqn
# df <- mutate(df, ave_GE = map(GE_mat, ~log1p(rowMeans(expm1(.x)))))
# 
# # get spqn-corrected weighted correlation
# df <- mutate(df, cor.spqn = map2(cor.raw, ave_GE, get_corrected_cor))
# 
# # get pvals and coefs and annotations in tbl format
# df <- mutate(df, res.raw = map2(cor.raw, n, cormat2tbl),
#              res.spqn  = map2(cor.spqn, n, cormat2tbl) )
# 
# 
# df |>
#   dplyr::select(res.spqn) |>
#   unnest(res.spqn) |>
#   filter(str_detect(feature,"FBgn")) |>
#   filter(!str_detect(y,"FBgn")) |>
#   mutate(tf = if_else(feature=="FBgn0085432","pan","other")) |>
#   ggplot(aes(tf,coef)) +
#   geom_boxplot() +
#   ggpubr::stat_compare_means() +
#   ylab("weighted correlation") +
#   labs(title="TE/TF correlations")
