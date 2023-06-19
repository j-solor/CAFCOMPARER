library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readxl)
library(GSVA)
library(pheatmap)
library(biomaRt)

################################################################################
sample_ID = "private" # private | public
signature_selection = NULL  # Luo_NatCom2022 | Foster_CancerCell2022 | Huang_CancerCell2022 | Verginadis_NatureCellBio2022 | Grauel_NatComms2020 | Carpenter_CancerDiscovery2023
include_tech_annot = FALSE # TRUE | FALSE
whatever_genes_filename = "Receptors_HGNC.csv" # Receptors_info.csv "whatever csv ythat uses tab as separator and has columns of groups (avoid spaces in names)
################################################################################
setwd("/home/margaux.dore/Documents/CAFCOMPARER")

# Functions 
source("src/Exclude_0s.R")
source("src/ID_converter.R")
source("src/mouseID_to_humanID.R")

# Data loading
## CAF
acafs.info <- read_tsv("data/CAF_metadata.tsv") # metadata about CAFs 
acafs.raw <- read_tsv("data/CAF_rawcounts.tsv") %>% # count table
  dplyr::select(EnsemblID, acafs.info$rawID) %>%
  dplyr::relocate(EnsemblID, acafs.info$rawID)

## Signatures
if (is.null(signature_selection)){ # initialize at NULL
  signatures = tibble(signature = character(), value = character())
  sign_list <- excel_sheets("data/CAF_signatures.xlsx") # list of sheet's names
for (sign in sign_list){
  signatures_temp <- read_xlsx("data/CAF_signatures.xlsx", sheet = sign) %>%
    rename_with( ~ paste0(sign,"_", .x)) %>%
    pivot_longer(cols = everything(), names_to = "signature")
  
  signatures <- bind_rows(signatures, signatures_temp) # table with gene name in each sheet
}

} else {
  signatures <-  read_xlsx("data/CAF_signatures.xlsx", sheet = signature_selection) %>% 
    pivot_longer(cols = everything(), names_to = "signature")
  
}

## Call the mouseID_to_humanID to only have human genes in signatures
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl") # the mart
signatures_human <- mouseID_to_humanID(signatures, ensembl)

list_of_signatures <- split(signatures_human, f = signatures_human$signature) %>% # list of lists of gene names in each sheet
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

acafs.raw <-column_to_rownames(acafs.raw, "EnsemblID") # count table (EnsemblID)

## Normalize
acafs.tmm <- Exclude_0s(acafs.raw, 0.5) %>% DGEList() %>% # normalized count table
  calcNormFactors(method = "TMM") %>% cpm() 

## convert ID
annot <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(acafs.tmm), columns="SYMBOL", keytype="ENSEMBL") # Annotation table
cafs.choose.sym <- ID_converter(df = acafs.tmm,annotation_table = annot,
                                old_IDs = "ENSEMBL", new_IDs = "SYMBOL")

# Analysis
## Barplot of whatever gene
gene = "COL1A1"
cafs.gene.expression <- as_tibble(cafs.choose.sym, rownames = "gene") %>% 
  pivot_longer(cols = -gene ,names_to = "CL",values_to = "expression" ) %>%
  pivot_wider(id_cols = "CL", names_from = "gene", values_from = "expression")

dplyr::mutate(cafs.gene.expression, CL = fct(CL, levels = arrange(cafs.gene.expression, desc(get(gene)))$CL)) %>%
  ggplot(aes(x = get(gene), y = CL)) + 
  geom_bar(stat = "identity") +
  labs(x = "COL1A1", y = "CAFs") +
  theme_bw()

## GSEA of signatures
gsvaRes <- gsva(cafs.choose.sym %>% data.matrix(), list_of_signatures, min.sz = 5)

### pheatmap of all the CAF subtypes
if (include_tech_annot == TRUE){
  full_colanot <- c("proliferation", "passes", "Inmortalized", "tissue_origin",
                    "sex", "age", "tumor_differentiation", "initial_sample",
                    "treatment_neoadjuvant", "OS")
  tech_annot = acafs.info[full_colanot]
  
} else {tech_annot = NULL}

pheatmap(gsvaRes,
         cluster_cols = T, cluster_rows = T, annotation_col = tech_annot)

### Specific pheatmaps per signature set (output in files)

#### Convert the mouse gene names into their human ortholog

for (sign in sign_list){
  sign_sub <- dplyr::filter(signatures_human, str_detect(signature, sign)) %>% # filter the signatures == sign
    dplyr::filter(value %in% rownames(cafs.choose.sym)) %>% # only keep the values present in cafs.choose.sym
    group_by(value) %>% 
    mutate(signature = ifelse(base::duplicated(value,fromLast=T), "multiple", signature)) %>%
    distinct(value, .keep_all = T,) %>%
    dplyr::arrange(signature) %>%
    column_to_rownames("value")
  
  gsvaRes_sub <- as.data.frame(t(gsvaRes[unique(sign_sub$signature)[!unique(sign_sub$signature) == 'multiple'],]))
  
  if((gsvaRes_sub %>% summarize(count = n())) == 1) {
    rownames(gsvaRes_sub) <- sign
  }
  
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

### Violin plot

tb_cafs_GE <- cafs.gene.expression %>%
  pivot_longer(cols = -CL, names_to = "genes")

whatever_genes_violin <- left_join(whatever_genes, tb_cafs_GE, by = c("value" = "genes")) %>%
  dplyr::rename(genes = value, value = value.y, Signature = signature) 
#%>% dplyr::filter(!is.na(CL))

levels <- group_by(whatever_genes_violin, CL) %>% summarise(mean=mean(value)) %>% arrange(dplyr::desc(mean))
tb_organized <- whatever_genes_violin %>% dplyr::mutate(CL = fct(CL, levels = levels$CL))

ggplot(tb_organized, aes(x=CL, y=value, fill=Signature)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width=0.1) +
  labs(x = "CAFs", y = "Genes")

#### Boxplot 

ggplot(tb_organized, aes(x=CL, y=value, fill=Signature)) +
  geom_boxplot() +
  labs(x = "CAFs", y = "Genes")
