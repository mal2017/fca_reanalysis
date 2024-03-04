library(SuperCell)
library(tidyverse)
library(scater)
library(scran)
library(rtracklayer)
library(batchelor)

# See https://github.com/GfellerLab/SuperCell/blob/master/vignettes/a_SuperCell.Rmd
# for workflow I recreate here

#sce_fls <- Sys.glob("~/amarel-matt/fca_reanalysis/results/import/fca_annotated_r_objs/FCA*fca_annot.sce.rds")[1:3] %>%
#  set_names(.,str_extract(.,"FCA\\d+"))

sce_fls <- snakemake@input$sce %>%
  set_names(.,str_extract(.,"FCA\\d+"))

# make unique cell names
sces <- map(sce_fls,~{x<-read_rds(.x)
colnames(x) <- paste(x$FCAID,colnames(x),sep="_")
x <- logNormCounts(x)
x})

# combine into one sce
sces <- do.call(cbind,sces)

# normalize
sces <- multiBatchNorm(sces,batch = sces$FCAID)

# from osca:
# Our recommendation is to simply pick an arbitrary n 
# and proceed with the rest of the analysis, with the intention
# of testing other choices later, rather than spending much time 
# worrying about obtaining the “optimal” value.
dec <- modelGeneVar(sces)
hvg <- getTopHVGs(dec, n=1000)
hvg <- hvg[str_detect(hvg,"FBgn")]

#gamma <- 20 # graining level
gamma <- snakemake@params$gamma

#k.knn <- 5
k.knn <- snakemake@params$k_knn

# for local testing
#sces <- sces[,sample(1:ncol(sces),size = round(0.2*ncol(sces)),replace = F)]

SC <- SCimplify(logcounts(sces),  # gene expression matrix 
                #cell.split.condition=sces$time,
                k.knn = k.knn, # number of nearest neighbors to build kNN network
                gamma = gamma, # graining level
                genes.use = hvg,
)

# by default, this assumes input counts are normalized and log transformed
SC.GE <- supercell_GE(logcounts(sces), SC$membership)

SC$lineage <- supercell_assign(clusters = sces$annotation, # single-cell assigment to cell lines (clusters)
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
SC.clusters    <- supercell_cluster(D = D, k = length(unique(sces$annotation)), supercell_size = SC$supercell_size) 
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