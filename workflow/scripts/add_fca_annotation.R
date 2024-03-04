library(scater)
library(scran)
library(scuttle)
library(tidyverse)

fca_id <- snakemake@params$fca_id

#j_fl <- "resources/icmuj3_parsing_output.loom_col_attrs_Embedding_SCENIC_AUC_UMAP.json"
j_fl <- snakemake@input$j
j <- jsonlite::read_json(j_fl,simplifyVector = F)

a_fl <- "resources/icmuj3_parsing_output.loom_col_attrs_R_annotation.tsv"
a_fl <- snakemake@input$a
a <- read_tsv(a_fl,col_names = c("cell","annotation"))

df <- tibble(cell=unlist(j$cells),
       UMAP1=unlist(j$values[[1]]),
       UMAP2=unlist(j$values[[2]])) |>
  full_join(a,by="cell")

df <- df |>
  mutate(barcode = str_extract(cell,".+(?=-)")) |>
  mutate(hash=str_extract(cell,"(?<=-).+(?=__)")) |>
  mutate(FCAID=str_extract(cell,"FCA\\d+")) |>
  mutate(tissue = map2_chr(cell,FCAID,~str_extract(.x,sprintf("(?<=%s_).+",.y))))

# confirm this looks like UMAP shown by the ASAP portal for proj ASAP49/icmuj3
#ggplot(arrange(df,desc(annotation)),aes(UMAP1,UMAP2,color=annotation)) +
#  geom_point() +
#  guides(color="none")

#sce_fl <- "results/downstream/raw_r_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5.usa.sce.rds"
sce_fl <- snakemake@input$sce
sce <- read_rds(sce_fl)

# the cell is composed of the following BARCODE-HASH?__FCAID_SEX_TISSUE
# the important ones needed to uniquely join are the barcode and the FCAID
colData(sce) <- sce |> colnames() |>
  tibble(barcode=_) |>
  mutate(FCAID=fca_id) |>
  mutate(barcode2=barcode) |>
  left_join(df,by=c("barcode","FCAID")) |>
  column_to_rownames("barcode2") |>
  DataFrame()

# remove cells not included in FCA's analysis
sce <- sce[,!is.na(sce$cell) & !is.na(sce$UMAP1)]

reducedDim(sce,"UMAP") <- colData(sce)[,c("UMAP1","UMAP2")]

write_rds(sce,snakemake@output$sce)