
thekappa <- function (df1, df2, colname, isweighted=FALSE) {
  # This function calculate the Kappa for one column in df1 and df2. It returns
  # a Kappa object in library(vcd)
  #
  # Inputs: 
  #    df1, df2; dataframes
  #    colname: the name of the column to be inspected
  #    isweighted: whether save the weighted results. Defaults to FALSE.
  #
  # Output:
  #    Kappa object in library(vcd)
  
  w_strs = c("Unweighted", "Weighted")
  df_anainfo <- tibble(Colname = colname,
                       method = w_strs[isweighted+1])
  
  thelevels <- sort(unique(tolower(c(pull(df1, colname), pull(df2, colname)))))
  
  # manual add one pair of consistent responses (to be removed)
  tmp_df <- tibble(col1 = thelevels, col2 = thelevels)
  
  conf_matrix <- bind_cols(transmute(df1, col1=tolower(df1[[colname]])),
                           transmute(df2, col2=tolower(df2[[colname]]))) %>% 
    bind_rows(tmp_df) %>% 
    group_by(col1, col2) %>% 
    summarize(count= n(), .groups="drop") %>% 
    pivot_wider(names_from=col1, values_from=count, values_fill=0) %>% 
    arrange(col2) %>% 
    column_to_rownames("col2")
  
  # remove the manually added 1s
  tmp_mat <- as.matrix(conf_matrix) - diag(length(thelevels))
  # calculate Kappa
  tmp_K <- Kappa(tmp_mat)
  
  if (nrow(tmp_mat)*ncol(tmp_mat)==1) {
    # when both raters only choose the same 1
    outdf <- tibble(
      value = 1,
      ASE = 0,
      lwr = 1,
      upr = 1
    )
    
  } else {

      outdf <- bind_cols(
      as_tibble(t(tmp_K[[w_strs[isweighted+1]]])), # Kappa and SE
      as_tibble(confint(tmp_K))[isweighted+1,] # CI
    ) 
  }
  
  # save the output
  outdf <- bind_cols(df_anainfo, outdf) %>% 
    mutate(conf_mat = list(tmp_mat))
  
  return(outdf)
}
