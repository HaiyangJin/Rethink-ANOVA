# function for simulating data
simu_null <- function(iter=100, N=30, rho=.5, sd=1, n_core=2,
                      file_cache=NULL){
  # iter: iteration number for each sample size
  # N: sample size in within-subjects design. sample size in each condition in between-subject design.
  # rho: the correlations among conditions
  # sd: standard deviations for each condition
  # n_core: number of cores to run the simulation
  # file: file name of the cache file. If NULL, no file will be cached.
  
  # source the function for single iteration
  simu_null_single <- function(N_s, rho_s, sd_s){
    library(tidyverse)
    library(emmeans)
    
    mean <- rep(3, 4)
    if (rho_s == 0){
      # if it is between-subjects design
      df_null <- tibble(
        Subj = 1:(N_s*4), 
        Congruency = rep(c("congruent", "incongruent"), each = N_s*2),
        Alignment = rep(c("aligned", "misaligned"), times = N_s*2),
        d = rnorm(N_s*4, mean[1], sd_s) 
      )
      
      aov_tmp <- aov(d ~ Congruency * Alignment, data = df_null)
      
    } else {
      # if it is within-subject design
      C <- matrix(rho_s, nrow = 4, ncol = 4)
      diag(C) <- 1
      
      df_null <- mvtnorm::rmvnorm(N_s, mean, sigma = C*sd_s^2) %>% 
        as_tibble(.name_repair = NULL) %>% 
        transmute(Subj = 1:n(),
                  congruent_aligned = V1,  incongruent_aligned = V2,
                  congruent_misaligned = V3, incongruent_misaligned = V4) %>% 
        pivot_longer(contains("_"), names_to = c("Congruency", "Alignment"), 
                     values_to = "d", names_sep = "_")
      
      aov_tmp <- afex::aov_4(d ~ Congruency * Alignment + (Congruency * Alignment | Subj),
                             data = df_null) 
    }
    
    # simple effect 
    tmp_simple <- contrast(emmeans(aov_tmp, ~Congruency|Alignment), "pairwise") %>% 
      as_tibble() %>% 
      filter(Alignment == "aligned")
    
    tmp_inter <- contrast(emmeans(aov_tmp, ~Congruency+Alignment),
                          interaction="pairwise") %>% 
      as_tibble()
    
    tmp_main1 <- contrast(emmeans(aov_tmp, ~Congruency), "pairwise") %>% 
      as_tibble()
    
    tmp_main2 <- contrast(emmeans(aov_tmp, ~Alignment), "pairwise") %>% 
      as_tibble()
    
    p_simple <- tmp_simple$p.value
    p_inter <- tmp_inter$p.value
    p_main1 <- tmp_main1$p.value
    p_main2 <- tmp_main2$p.value
    
    rawp <- tibble(N = N_s,
                   p_simple = p_simple,
                   p_inter = p_inter,
                   p_main1 = p_main1,
                   p_main2 = p_main2,
                   rho = rho_s,
                   sd = sd_s,
                   N_row = length(unique(df_null$Subj)))
    
    return(rawp)
  }
  
  if (!is.null(file_cache) && file.exists(as.character(file_cache))){
    df_simu <- read_rds(file_cache)
  } else {
    # set parallel processing
    Ns_iter <- rep(N, times=iter)
    ls_tibble <- pbapply::pblapply(Ns_iter, simu_null_single,  
                                   rho_s=rho, sd_s=sd, cl=n_core)
    
    df_simu <- reduce(ls_tibble, rbind) %>% 
      dplyr::mutate(iter=1:n())
    
    if (!is.null(file_cache)){
      write_rds(df_simu, file=file_cache)
    }
  }
  
  return(df_simu)
}

simu_alpha <- function(df_simu, alpha=.05, disp=TRUE){
  # function to calculate the Type I error rate. 
  # df_simu: dataframe obtained from simu_null
  # alpha: alpha level used to calculate the Type I error rate.
  # disp: if true, the output will be a wide format. If false, it will be long.
  
  df_sig <- df_simu %>% 
    mutate(sig_simple = p_simple <= alpha,
           sig_inter = p_inter <= alpha,
           sig_both = sig_simple * sig_inter) %>% 
    select(N, starts_with("sig_")) %>% 
    pivot_longer(starts_with("sig_"), names_to = "effects", values_to = "isSig") %>% 
    group_by(N, effects) %>% 
    summarize(N_iter = n(),
              N_sig = sum(isSig),
              Type_I = mean(isSig),
              .groups="keep")
  
  if (disp){
    df_sig <- df_sig %>% 
      select(-c(N_sig)) %>% 
      pivot_wider(c(N, N_iter), names_from = "effects", values_from = "Type_I")
  }
  
  return(df_sig)
}