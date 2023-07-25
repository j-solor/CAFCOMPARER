# Customizable corrplot (modified from https://www.khstats.com/blog/corr-plots/corr-plots/)
cors <- function(df, cor.stat) {
  M <- Hmisc::rcorr(as.matrix(df), type = cor.stat)
  Mdf <- map(M, ~data.frame(.x))
  return(Mdf)
}

formatted_cors <- function(df, cor.stat){
  cors(df, cor.stat) %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA)) 
}


dfs_corrplot <- function(df_x, df_y, abscorr = FALSE) {
  
  stopifnot(rownames(df_x)==rownames(df_y))
  joint_df <- df_y %>% bind_cols(df_x)
  
  df_corr <- rcorr(data.matrix(joint_df), type = "spearman")
  df_corr$r <- df_corr$r[colnames(df_x), colnames(df_y)]
  df_corr$P <- df_corr$P[colnames(df_x), colnames(df_y)]
  
  if (abscorr == TRUE){
    cmap = c("white", "white","red")
    scale_lims = c(0,1)
  } else {
    cmap = c("#7D0E29", "white", "#004376")
    scale_lims = c(-1,1)
  }
  
  ggcorrplot(df_corr$r, p.mat = df_corr$P) +
    scale_fill_gradient2(limit = scale_lims, low = cmap[1], mid = cmap[2], high =  cmap[3])
  
}