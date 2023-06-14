#' Convert the IDs of a gene expression table based inan annotation table 
#'
#' @description
#' `ID_converter` takes a dataframe with some define IDs as rownames and translate these to others based in an annotation table and strings determining the new and old ID.
#' 
#' @usage 
#' @param df dataset with ProbeIDs as rownames
#' @param annotation_table  Bioconductor annotation table, or any other source 
#' containing the original and the desired IDs correspondance ie. 
#' AnnotationDbi::select(hgu219.db, probes, c("SYMBOL", "ENSEMBL", "GENENAME"))
#' @param old_IDs string matching a name of `annotation_table` and rownames of `df`
#' @param new_IDs string matching a name of `annotation_table` and desired new 
#' rownames of `df`.
#'
#'@return
#'`ID_converter` will return a data.frame with `new_IDs` equivalents as the 
#'rownames. £ messages indicating statistics of this conversion will be produced
#' indicating the % of rownames merged or withourt equivalents.

ID_converter <- function(df, 
                         annotation_table, 
                         old_IDs,
                         new_IDs 
)
  #TBI: choose function of aggregation
{
  final_df <- merge(df,annotation_table, by.x=0, by.y=old_IDs)
  
  # Stats of conversion
  non_agg <- nrow(final_df) # nb of genes 
  non_agg_uniq <- length(unique(final_df[new_IDs])) # nb of unique 
  non_agg_nas <- sum(is.na(final_df[new_IDs])) # nb of gene with no new_ID 
  non_agg_nonas <- non_agg-non_agg_nas # 
  
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
