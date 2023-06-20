
whatever_genes_function <- function(file_path, expression_table, scale_option, output_file) {
  
  # Main
  whatever_genes <- read_tsv(file_path) %>%
    pivot_longer(cols = everything(), names_to = "signature") %>%
    dplyr::filter(!is.na(value))
  
  list_of_whatever_genes <- split(whatever_genes, f = whatever_genes$signature) %>%
    map(~ .$value)
  
  gsvaRes <- gsva(data.matrix(expression_table), list_of_whatever_genes)
  
  # Plots
  ## Heatmap
  gene_anotation <- dplyr::filter(whatever_genes, value %in% rownames(expression_table)) %>% 
    group_by(value) %>% 
    mutate(signature = ifelse(base::duplicated(value,fromLast=T),
                              "multiple", signature)) %>%
    distinct(value, .keep_all = T,) %>%
    dplyr::arrange(signature) %>%
    column_to_rownames("value")
  
  expression_table[rownames(gene_anotation),] %>% 
    pheatmap(cellwidth=15, cellheight=15, filename = output_file,
             cluster_cols = T, cluster_rows = F, scale = scale_option,
             show_rownames = T, annotation_col = as.data.frame(t(gsvaRes)),
             annotation_row = gene_anotation)
  
  ## Boxplot
  expression_table_boxplot <- as_tibble(expression_table, rownames = "gene") %>% # tb_cafs_GE
    pivot_longer(cols = -gene , names_to = "CL",values_to = "value" )
  
  whatever_genes_boxplot <- left_join(whatever_genes, expression_table_boxplot, 
                                      by = c("value" = "gene")) %>% 
    dplyr::rename(gene = value, value = value.y, Signature = signature) %>%
    dplyr::filter(!is.na(CL))
  
  levels <- group_by(whatever_genes_boxplot, CL) %>% summarise(mean=mean(value)) %>%
    arrange(dplyr::desc(mean))
  whatever_genes_organized <- whatever_genes_boxplot %>%
    dplyr::mutate(CL = fct(CL, levels = levels$CL)) # tb_organized
  
  dev.off()
  
  boxplot <- ggplot(whatever_genes_organized, aes(x=CL, y=value, fill=Signature)) + # boxplot_whatever_genes
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = "CAFs", y = "Gene expression")
  
  print(boxplot)
  
  genes <- whatever_genes %>% dplyr::pull(value)
  heatmap_genes <- rownames(expression_table[rownames(gene_anotation),])
  missing_genes <- setdiff(genes, heatmap_genes)
  
  if(length(missing_genes) > 0) {
    cat("The following genes were not in the expression table :",
        paste(missing_genes, collapse = ", "))
  }
  
}
