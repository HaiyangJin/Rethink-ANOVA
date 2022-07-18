
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

  conf_matrix <- cbind(transmute(df1, col1=df1[[colname]]),
        transmute(df2, col2=df2[[colname]])) %>% 
    group_by(col1, col2) %>% 
    summarize(count= n(), .groups="drop") %>% 
    pivot_wider(names_from=col2, values_from=count, values_fill=0) 
  
  tmp_matrix <- select(conf_matrix, -col1) 
  
  tmp_K <- tmp_matrix %>% 
    as.matrix() %>% 
    Kappa()
  
  if (nrow(tmp_matrix)*ncol(tmp_matrix)==1) {
    
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
  
  outdf <- bind_cols(df_anainfo, outdf) %>% 
    mutate(conf_mat = list(conf_matrix))
  
  return(outdf)
}
