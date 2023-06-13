#' Filter a dataframe to keep genes with at least a defined % of non 0 expression samples
#'
#' @description
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
