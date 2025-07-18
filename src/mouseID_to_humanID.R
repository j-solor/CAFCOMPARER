#' Translate mouse gene IDs to their human ortholog gene IDs  
#'
#' @description
#' `mouseID_to_humanID` requires a data frame with mouse gene names in one column
#' and replaces them with the gene names corresponding to their human ortholog 
#' using biomart.
#' 
#' @usage 
#' @param df data frame with mouse gene names in a column called value and the 
#' signature corresponding to each gene in a column called signature.
#' This data frame can be taken from the excel file CAF_signatures.xlsx.
#' @param mart the name of the mart that will be used to find the human orthologs 
#' with biomart.  
#' 
#' @return
#' `mouseID_to_humanID` returns a data.frame similar to `df` where the mouse gene
#' names are replaced with the gene names of their human ortholog. 
#' Warning messages indicates the mouse genes that were not found in the database,
#' the mouse genes without human orthologs or the mouse genes with multiple human
#' orthologs.

mouseID_to_humanID <- function(df, mart) {
  
  # Get the mouse gene names from the df in argument 
  mouse_genes <- c()
  for (val in df$value) {
    if(!is.na(val) && str_detect(val, "[[:lower:]]")) {
      mouse_genes <- c(mouse_genes, val) 
    }
  }
  
  # Run biomart
  human_id = getBM(attributes = c("ensembl_gene_id", "external_gene_name","hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
                   filters = "external_gene_name",
                   values = mouse_genes,
                   mart = mart)
  
  # Mouse genes not present in the database
  not_in_db <- setdiff(mouse_genes, human_id$external_gene_name)
  
  # Mouse genes without human orthologs
  without_ortho <- dplyr::filter(human_id, hsapiens_homolog_ensembl_gene == "") %>% 
    dplyr::pull(external_gene_name)
  
  # Filter the mouse genes without human orthologs 
  human_id_final <- dplyr::filter(human_id, hsapiens_homolog_ensembl_gene != "") 
  
  # Mouse genes with multiple orthologs
  multiple_ortho <- dplyr::filter(human_id_final, duplicated(external_gene_name)) %>% 
    dplyr::pull(external_gene_name)
  
  # Modify the df in argument : convert the mouse ID into human ID
  res <- left_join(df, human_id_final, by = c("value" = "external_gene_name"), relationship = "many-to-many") %>%
    mutate(value = ifelse(is.na(hsapiens_homolog_associated_gene_name), value, hsapiens_homolog_associated_gene_name)) %>%
    dplyr::select(signature, value) %>%
    dplyr::filter(!is.na(value))
  
  # Filter fom the result table the mouse genes that were not in the database
  final <- anti_join(res, data.frame(value = c(not_in_db, without_ortho)), by = "value")
  
  # Warning messages 
  if(length(not_in_db) > 0) {
    cat("The following mouse genes were not found in the ensembl database :", paste(not_in_db, collapse = ", "), "\n\n")
  }
  
  if(length(without_ortho) > 0) {
    cat("The following mouse genes don't have human orthologs :", paste(without_ortho, collapse = ", "), "\n\n")
  }
  
  if(length(multiple_ortho) > 0) {
    cat("The following mouse genes have multiple human orthologs :", paste(multiple_ortho, collapse = ", "))
  }
  
  return(final)
}
