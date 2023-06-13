#' Convert gene IDs into desired gene IDs  
#'
#' @description
#' `ID_converter` requires a data.frame with genes as rows,
#'  and its ids as rownames. returns a dataframe

ID_converter <- function(df, # dataset with ProbeIDs as rownames
                         annotation_table, # Bioconductor annotation table ie. AnnotationDbi::select(hgu219.db, probes, c("SYMBOL", "ENSEMBL", "GENENAME"))
                         old_IDs, # Current IDs ("SYMBOL", "ENSEMBL", "GENENAME")
                         new_IDs # Desired IDs ("SYMBOL", "ENSEMBL", "GENENAME")
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
