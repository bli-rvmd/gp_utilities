## Provided a list of gene names, here's script to query depmap CRISPRGeneEffect.csv downloaded from https://depmap.org/portal/data_page/?tab=currentRelease
## and to annotate model/cell lines from Model.csv, also downloaded from https://depmap.org/portal/data_page/?tab=currentRelease

library(dplyr)

df_cr <- read.csv("~/Downloads/CRISPRGeneEffect.csv", header = T)
colnames(df_cr) <- sub("\\..*", "", colnames(df_cr))
colnames(df_cr)[1] <- "ModelID"

library(readxl)

gene_names <- as.vector(read_excel("../data/shSHOC2 off-target gene list.xlsx", col_names = FALSE)[[1]])

df_cr_selected <- df_cr[, colnames(df_cr) %in% c("ModelID", gene_names)]

df_model <- read.csv("~/Downloads/Model.csv", header = T)

df_res <- df_cr_selected %>%
  left_join(
    df_model %>%
      select(ModelID, StrippedCellLineName), 
    by = "ModelID"
  ) %>%
  select(ModelID, StrippedCellLineName, everything())

write.csv(df_res, file = "../data/shSHOC2_crispr_off-target_gene_effect.csv", row.names = F)