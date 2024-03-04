library(scater)
library(scran)
library(tidyverse)
library(ggVennDiagram)
library(patchwork)

# script to check that my assumptions about the 'runs' within each ENA "library"
# correspond to same physical library

l001 <- read_rds("results/import/fca_annotated_r_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5.usa.fca_annot.sce.rds")

l003 <- read_rds("results/import/fca_annotated_r_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5_sample2.usa.fca_annot.sce.rds")


g <- list(L001=l001$barcode,
     L003=l003$barcode) |>
  ggVennDiagram() + labs(title="after intersection and annotation w/ FCA metadata")


l001_noanno <- read_rds("results/import/raw_r_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5.usa.sce.rds")

l003_noanno <- read_rds("results/import/raw_r_objs/FCA14_Female_head_adult_5dWT_Luo_sample5_S5_sample2.usa.sce.rds")


g_noanno <- list(L001=colnames(l001_noanno),
     L003=colnames(l003_noanno)) |>
  ggVennDiagram() + labs(title="straight out of alevin-fry")


g_p <- g_noanno + g + patchwork::plot_annotation(title="barcode overlaps for library FCA14 (run L001, L003)")

ggsave("config/library_assumption_check/library_assumption_check.svg",g_p)
