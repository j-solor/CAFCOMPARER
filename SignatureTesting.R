library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readxl)
library(GSVA)
library(pheatmap)
################################################################################
sample_ID = "private" # private | public
signature_selection = NULL  # Luo_NatCom2022 | Foster_CancerCell2022 | Huang_CancerCell2022 | Verginadis_NatureCellBio2022 | Grauel_NatComms2020 | Carpenter_CancerDiscovery2023
include_tech_annot = FALSE # TRUE | FALSE
whatever_genes_filename = "Receptors_HGNC.csv" # Receptors_info.csv "whatever csv ythat uses tab as separator and has columns of groups (avoid spaces in names)
################################################################################
#' Filter a dataframe to keep genes with at least a defined % of non 0 expression samples
#'@description
#' `Exclude_0s` requires a data.frame with genes as rows,
#'  and its ids as rownames. returns a dataframe

Exclude_0s <- function(df, threshold){
  
  n0s <- apply(df, 1, function(x) {
    sum(
      x > 0
    )/length(x)
  })
  
  
  df.f <- df[which(n0s > threshold), ]
  return(df.f)
}

ID_converter <- function(df, # dataset with ProbeIDs as rownames
                         annotation_table, # Bioconductor annotation table ie. AnnotationDbi::select(hgu219.db, probes, c("SYMBOL", "ENSEMBL", "GENENAME"))
                         old_IDs, # Current IDs ("SYMBOL", "ENSEMBL", "GENENAME")
                         new_IDs # Desired IDs ("SYMBOL", "ENSEMBL", "GENENAME")
)
  #TBI: choose function of aggregation
{
  final_df <- merge(df,annotation_table, by.x=0, by.y=old_IDs)
  
  # Stats of conversion
  non_agg <- nrow(final_df)
  non_agg_uniq <- length(unique(final_df[new_IDs]))
  non_agg_nas <- sum(is.na(final_df[new_IDs]))
  non_agg_nonas <- non_agg-non_agg_nas
  
  # Conversion
  final_df <- aggregate(final_df,
                        final_df[new_IDs],
                        FUN = mean) %>% 
    column_to_rownames(new_IDs) %>% 
    subset(select = colnames(df))
  
  agg <- nrow(final_df)
  
  print(paste(100*(non_agg-agg)/non_agg,
              "% of the originally merged df have been agregated/removed"))
  
  print(paste(100*(non_agg_nas)/non_agg,
              "% of the originally merged df have been removed due no", new_IDs))
  
  print(paste(100*(non_agg_nonas - agg)/non_agg,
              "% of the non NAs df have been aggregated"))
  
  return(final_df)
}
################################################################################
setwd("/home/jacobo/Documents/02_TRANSDUCER/04_Deconvolution_TE/03_Sigantures/cell_types/fibroblast/PJ2003085/260523_CAFComparer/")

# Data loading
## CAF
acafs.info <- read_tsv("data/CAF_metadata.tsv")
acafs.raw <- read_tsv("data/CAF_rawcounts.tsv") %>% 
  dplyr::select(EnsemblID, acafs.info$rawID) %>%
  dplyr::relocate(EnsemblID, acafs.info$rawID)


## Signatures
if (is.null(signature_selection)){
  signatures = tibble(signature = character(), value = character())
  sign_list <- excel_sheets("data/CAF_signatures.xlsx")
for (sign in sign_list){
  signatures_temp <- read_xlsx("data/CAF_signatures.xlsx", sheet = sign) %>%
    rename_with( ~ paste0(sign,"_", .x)) %>%
    pivot_longer(cols = everything(), names_to = "signature")
  
  signatures <- bind_rows(signatures, signatures_temp) 
}

}else {
  signatures <-  read_xlsx("data/CAF_signatures.xlsx", sheet = signature_selection) %>% 
    pivot_longer(cols = everything(), names_to = "signature")
  
}

list_of_signatures <- split(signatures, f = signatures$signature) %>%
  map(~ .$value)

# Preprocessing
## choose sample IDs
all(colnames(acafs.raw)[-1] == acafs.info$rawID) %>% stopifnot()
if (sample_ID == "public") {
  names(acafs.raw)[-1] <-  acafs.info$Official_name
  acafs.info <- column_to_rownames(acafs.info, "Official_name")

} else if (sample_ID == "private"){
  names(acafs.raw)[-1] <-  acafs.info$Name
  acafs.info <- column_to_rownames(acafs.info, "Name")
}
acafs.raw <-column_to_rownames(acafs.raw, "EnsemblID")

## Normalize
acafs.tmm <- Exclude_0s(acafs.raw, 0.5) %>% DGEList() %>% 
  calcNormFactors(method = "TMM") %>% cpm() 

## convert ID
annot <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(acafs.tmm), columns="SYMBOL", keytype="ENSEMBL")
cafs.choose.sym <- ID_converter(df = acafs.tmm,annotation_table = annot,
                                old_IDs = "ENSEMBL", new_IDs = "SYMBOL")

# Analysis
## Barplot of whatever gene
gene = "SHMT1"
cafs.gene.expression <- as_tibble(cafs.choose.sym, rownames = "gene") %>% 
  pivot_longer(cols = -gene ,names_to = "CL",values_to = "expression" ) %>%
  pivot_wider(id_cols = "CL", names_from = "gene", values_from = "expression")

dplyr::mutate(cafs.gene.expression, CL = fct(CL, levels = arrange(cafs.gene.expression, desc(get(gene)))$CL)) %>%
  ggplot(aes(x = get(gene), y = CL)) + 
  #aes_string(x = gene, y = "CL") +
  geom_bar(stat = "identity")


## GSEA of signatures
gsvaRes <- gsva(cafs.choose.sym %>% data.matrix(), list_of_signatures)

### pheatmap of all the CAF subtypes
if (include_tech_annot == TRUE){
  full_colanot <- c("proliferation", "passes", "Inmortalized", "tissue_origin",
                    "sex", "age", "tumor_differentiation", "initial_sample",
                    "treatment_neoadjuvant", "OS")
  tech_annot = acafs.info[full_colanot]
  
}else {tech_annot = NULL}


pheatmap(gsvaRes,
         cluster_cols = T, cluster_rows = T, annotation_col = tech_annot)


### Specific pheatmaps per signature set (output in files)
for (sign in sign_list){
  sign_sub <- dplyr::filter(signatures, str_detect(signature, sign)) %>%
    dplyr::filter(value %in% rownames(cafs.choose.sym)) %>% 
    group_by(value) %>% 
    mutate(signature = ifelse(base::duplicated(value,fromLast=T), "multiple", signature)) %>%
    distinct(value, .keep_all = T,) %>%
    dplyr::arrange(signature) %>%
    column_to_rownames("value")
  
  gsvaRes_sub <- as.data.frame(t(gsvaRes[unique(sign_sub$signature)[!unique(sign_sub$signature) == 'multiple'],]))
  
  cafs.choose.sym[rownames(sign_sub),] %>% 
    pheatmap(cellwidth=15, cellheight=15, filename = paste0("output/",sign,".png"),
             cluster_cols = T, cluster_rows = F, scale = "row", show_rownames = T, annotation_col = gsvaRes_sub,
             annotation_row = sign_sub)
}

## Whatever list of genes Enrichment and pheatmap
whatever_genes <- read_tsv(paste0("data/",whatever_genes_filename)) %>%
  pivot_longer(cols = everything(), names_to = "signature") %>%
  dplyr::filter(!is.na(value))

list_of_whatever_genes <- split(whatever_genes, f = whatever_genes$signature) %>%
  map(~ .$value)

### GSEA
gsvaRes_whatever_genes <- gsva(data.matrix(cafs.choose.sym), list_of_whatever_genes)

### Plot
gene_anotation <- dplyr::filter(whatever_genes, value %in% rownames(cafs.choose.sym)) %>% 
  group_by(value) %>% 
  mutate(signature = ifelse(base::duplicated(value,fromLast=T), "multiple", signature)) %>%
  distinct(value, .keep_all = T,) %>%
  dplyr::arrange(signature) %>%
  column_to_rownames("value")

cafs.choose.sym[rownames(gene_anotation),] %>% 
  pheatmap(cellwidth=15, cellheight=15, filename = paste0("output/",whatever_genes_filename,".png"),
           cluster_cols = T, cluster_rows = F, scale = "row", show_rownames = T, annotation_col = as.data.frame(t(gsvaRes_whatever_genes)),
           annotation_row = gene_anotation)
