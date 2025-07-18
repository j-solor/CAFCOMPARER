library(tidyverse)
library(edgeR)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(readxl)
library(GSVA)
library(pheatmap)
library(biomaRt)
library(ComplexUpset) #devtools::install_github("krassowski/complex-upset")
library(Hmisc)
library(ggpubr)
library(sva)

################################################################################
sample_ID = "private" # private | public
signature_selection = NULL  # Luo_NatCom2022 | Foster_CancerCell2022 | Huang_CancerCell2022 | Verginadis_NatureCellBio2022 | Grauel_NatComms2020 | Carpenter_CancerDiscovery2023
include_tech_annot = FALSE # TRUE | FALSE
whatever_genes_filename = "Receptors_HGNC.csv" # Receptors_info.csv "whatever csv ythat uses tab as separator and has columns of groups (avoid spaces in names)
################################################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # automatic set to script location

# Functions 
source("src/Exclude_0s.R")
source("src/ID_converter.R")
source("src/mouseID_to_humanID.R")
source("src/whatever_genes_analyser.R")
source("src/dfs_corrplot.R")

# Data loading
## Expression Data
acafs.info <- read_tsv("data/CAF_metadata.tsv") # metadata about CAFs 
acafs.raw <- read_tsv("data/CAF_rawcounts.tsv") %>% # count table
  dplyr::select(EnsemblID, acafs.info$rawID) %>%
  dplyr::relocate(EnsemblID, acafs.info$rawID)

## Signatures
if (is.null(signature_selection)){ # initialize at NULL
  signatures = tibble(signature = character(), value = character())
  sign_list <- excel_sheets("data/CAF_signatures.xlsx")[-1] # list of sheet's names
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

## Translate mice signatures to human
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

acafs.raw <- column_to_rownames(acafs.raw, "EnsemblID") # count table (EnsemblID)

## Normalize
acafs.tmm <- Exclude_0s(acafs.raw, 0.5) %>% DGEList() %>% # normalized count table
  calcNormFactors(method = "TMM") %>% cpm() 

## Batch effect correction
acafs.tmm.nobatch <- acafs.tmm %>% as.matrix() %>% 
  ComBat(batch = acafs.info$extraction_group, 
         par.prior=TRUE, 
         prior.plots=FALSE)

## convert ID
annot <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(acafs.tmm.nobatch), columns="SYMBOL", keytype="ENSEMBL") # Annotation table
cafs.choose.sym <- ID_converter(df = acafs.tmm.nobatch,annotation_table = annot,
                                old_IDs = "ENSEMBL", new_IDs = "SYMBOL")

################################################################################
# 1.  Analyze the enrichment of each CAF line in the different CAF subtypes
################################################################################
## GSEA of signatures
gsvaRes <- gsva(param = gsvaParam(exprData = cafs.choose.sym %>% data.matrix(), geneSets = list_of_signatures, minSize = 5))

### pheatmap of all the CAF subtype enrichments
if (include_tech_annot == TRUE){
  full_colanot <- c("proliferation", "passes", "Inmortalized", "tissue_origin",
                    "sex", "age", "tumor_differentiation", "initial_sample",
                    "treatment_neoadjuvant", "OS")
  tech_annot = acafs.info[full_colanot]
  
} else {tech_annot = NULL}

pheatmap(gsvaRes, filename = "output/GSEA_heatmap.pdf",
         cluster_cols = T, cluster_rows = T, annotation_col = tech_annot)

### Detailed pheatmaps per signature set (output in files)
for (sign in sign_list) {
  sign_sub <- dplyr::filter(signatures_human, str_detect(signature, sign)) %>% # filter the signatures == sign
    dplyr::filter(value %in% rownames(cafs.choose.sym)) %>% # only keep the values present in cafs.choose.sym
    group_by(value) %>% 
    mutate(signature = ifelse(base::duplicated(value,fromLast=T), "multiple", signature)) %>%
    distinct(value, .keep_all = T,) %>%
    dplyr::arrange(signature) %>%
    column_to_rownames("value")
  
  unique_sign <- unique(sign_sub$signature)[!unique(sign_sub$signature) == 'multiple']
  gsvaRes_filter <- as_tibble(gsvaRes, rownames = "signatures") %>%
    filter(signatures %in% unique_sign)
  
  if (nrow(gsvaRes_filter) != 0) {
    annotation_col <- gsvaRes_filter %>%
      column_to_rownames(var = "signatures") %>%
      t() %>%
      as.data.frame()
    
  } else {
    annotation_col = NULL
  }
  
  cafs.choose.sym[rownames(sign_sub),] %>% 
    pheatmap(cellwidth=15, cellheight=15, filename = paste0("output/",sign,".pdf"),
             cluster_cols = T, cluster_rows = F, scale = "row", show_rownames = T, annotation_col = annotation_col,
             annotation_row = sign_sub)
}


################################################################################
# Perform a comparative analysis of the different CAF signatures.
################################################################################
## Whatever list of genes Enrichment and pheatmap
whatever_genes_analyser(file_path = paste0("data/", whatever_genes_filename), 
                        expression_table = cafs.choose.sym, 
                        scale_option = "row",
                        output_file = paste0("output/", str_remove(whatever_genes_filename, "\\..*$"), "_heatmap.pdf"))
                        
ggsave(filename = paste0("output/", str_remove(whatever_genes_filename, "\\..*$"),"_boxplot.png")) # save the boxplot

# Comparative Analysis of CAF subtype signatures
## Quantitative analysis using Corrplots
gsvaRes_corr <- as_tibble(gsvaRes, rownames = "signature") %>%
  pivot_longer(names_to = "CL", cols = -signature) %>%
  pivot_wider(names_from = "signature") %>%
  column_to_rownames("CL") %>%
  formatted_cors(cor.stat = "pearson")

clust <- dplyr::select(gsvaRes_corr, measure1, measure2, r) %>% 
  pivot_wider(names_from=measure2, values_from = r) %>%
  column_to_rownames("measure1") %>%
  dist() %>%
  hclust()

ggplot(gsvaRes_corr, aes(measure1, measure2, fill=r, label=round(r_if_sig,2))) +
  geom_tile() +
  labs(x = NULL, y = NULL, fill = "Pearson's\nAbsolute\nCorrelation", title="Correlations Signatures",
       subtitle="Only significant correlation coefficients shown (95% I.C.)") +
  scale_fill_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
  geom_text() +
  theme_classic() +
  scale_x_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  scale_y_discrete(expand=c(0,0),limits = clust$labels[clust$order]) +
  ggpubr::rotate_x_text(angle = 90)


## Qualitative analysis using UpsetPlots 
### Data frame compatible with UpsetPlots
signatures_upset <- pivot_wider(signatures_human, 
                                names_from = "signature", 
                                values_from = "signature",
                                values_fill = 0,
                                values_fn = function(x) 1) %>%
  column_to_rownames(var = "value") %>%
  as.data.frame()

### List of the sets 
signature_sets <- colnames(signatures_upset)

### Colors 
list_colors <- hcl.colors(length(sign_list), palette = "viridis", alpha = NULL, rev = FALSE, fixup = TRUE)
color_articles <- unlist(map2(sign_list, list_colors, ~ rep(.y, sum(str_detect(signature_sets, .x)))))
list_colors_articles <- lapply(1:length(signature_sets), function(i) list(set = signature_sets[i], color = color_articles[i]))

### Upset Plot
ComplexUpset::upset(signatures_upset, signature_sets,
                    group_by = "sets",
                    mode = "exclusive_intersection", # or inclusive_intersection
                    min_size = 5, # minimum number of items in the displayed intersections 
                    min_degree = 2, # minimum number of groups involved in an intersection
                    n_intersections = 25, # number of intersections displayed on the plot (maximum)
                    keep_empty_groups=TRUE,
                    base_annotations=list(
                      'Intersection ratio'= intersection_ratio( # display the ratio
                        text=list(
                          vjust=0.5,
                          hjust=-0.1,
                          angle=90))),
                    set_sizes=(
                      upset_set_size()
                      + geom_text(aes(label=after_stat(count)), hjust=1, stat='count')), # display the number of genes in each set
                    queries = c(lapply(list_colors_articles, function(x) upset_query(group = x$set, color = x$color)), # colors
                                lapply(list_colors_articles, function(x) upset_query(set = x$set, fill = x$color))))

ggsave("output/Upset_plot.pdf", width = 22, height = 13) # save the Upset plot as pdf

################################################################################
# 3. Explore the enrichment and expression of other genes
################################################################################
## Barplot of whatever gene
gene = "COL1A1" # replace this by whatever gene symbol
cafs.gene.expression <- as_tibble(cafs.choose.sym, rownames = "gene") %>% 
  pivot_longer(cols = -gene ,names_to = "CL",values_to = "expression" ) %>%
  pivot_wider(id_cols = "CL", names_from = "gene", values_from = "expression")

dplyr::mutate(cafs.gene.expression, CL = fct(CL, levels = arrange(cafs.gene.expression, desc(get(gene)))$CL)) %>%
  ggplot(aes(x = get(gene), y = CL)) + 
  geom_bar(stat = "identity") +
  labs(x = gene, y = "CAFs") +
  theme_bw()

ggsave(filename = paste0("output/", gene,".png"))


## Whatever list of genes Enrichment and pheatmap
whatever_genes_analyser(file_path = paste0("data/", whatever_genes_filename), 
                        expression_table = cafs.choose.sym, 
                        scale_option = "row",
                        output_file = paste0("output/", str_remove(whatever_genes_filename, "\\..*$"), "_heatmap.pdf"))
                        
ggsave(filename = paste0("output/", str_remove(whatever_genes_filename, "\\..*$"),"_boxplot.png")) # save the boxplot
