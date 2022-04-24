
##### Omnibus ANOVA #####
sim_omnibus <- function(N_subj = 30, iter = 100, n_core=2, file_cache = NULL,
                         N_IV = 2, seed = 2022){
  # N_subj: number of participants per condition (between-subject design)
  # iter: number of simulation/iteration
  # n_core: number of cores to be used for simulation
  # file_cache: file name of the cache file. If NULL, no file will be cached.
  # N_IV: number of independent variables (each has two levels)
  # seed: seed used for simulation
  
  set.seed(seed)
  
  sim_omni_single <- function(N_IV){
    
    IV_list <- letters[1:N_IV]
    
    level_list <- lapply(IV_list, function(x){
      return(paste0(x,1:2))
    })
    
    df_IV <- expand.grid(level_list, KEEP.OUT.ATTRS = F)
    names(df_IV) <- toupper(IV_list)
    
    df_IV_subj <- do.call("rbind", replicate(N_subj, df_IV, simplify = FALSE))
    
    # simulate data
    df_omni_sim <- df_IV_subj %>% 
      mutate(Subj = 1:n(),
             DV = rnorm(n())) %>% 
      as_tibble()
    
    # formula
    formulastr <- paste0("DV ~ ", 
                         paste(sprintf("%s * ", toupper(IV_list[1:N_IV-1])), collapse = ""), 
                         toupper(IV_list[N_IV]),
                         " + (1|Subj)")
    
    aov_omni <- aov_4(as.formula(formulastr), data = df_omni_sim)
    
    # calculate omnibus F-test
    sumsq_all <- aov_omni$Anova$`Sum Sq`[-1]
    dfn_all <- aov_omni$Anova$Df[-1]
    
    sumsq_IV <- sum(sumsq_all[-length(sumsq_all)])
    sumsq_re <- sumsq_all[length(sumsq_all)]
    
    dfn_IV <- sum(dfn_all[-length(dfn_all)])
    dfn_re <- dfn_all[length(dfn_all)]
    
    # F-stat and p-value
    F_omni <- sumsq_IV/dfn_IV/(sumsq_re/dfn_re)
    p_omni <- pf(F_omni, dfn_IV, dfn_re, lower.tail = F)
    
    p.value <- c(aov_omni$Anova$`Pr(>F)`[c(-1, -length(aov_omni$Anova$`Pr(>F)`))], p_omni)
    effnames <- c(row.names(aov_omni$Anova)[c(-1, -length(aov_omni$Anova$`Pr(>F)`))], "omniF") 
    
    p_df <- tibble(effnames = effnames,
                   p.value = p.value) %>% 
      mutate(N_IV = N_IV)
    
    return(p_df)
  }
  
  # run simulation in parallel
  if (!is.null(file_cache) && file.exists(as.character(file_cache))){
    df_simu <- read_rds(file_cache)
    # message("Load simulation results from local cache files successfully.")
    
  } else {
    # set parallel processing
    Ns_IV <- rep(N_IV, times=iter)
    ls_tibble <- pbapply::pblapply(Ns_IV, sim_omni_single, cl=n_core)
    df_simu <- bind_rows(ls_tibble, .id = "iter")
    
    if (!is.null(file_cache)){
      write_rds(df_simu, file=file_cache)
    }
  }
  
  return(df_simu)
}

sig_omnibus <- function(df_simu_p, alphas=0.05){
  
  sig_omnibus_single <- function(df_simu_p, thealpha){
    
    df_simu_sig_long <- df_simu_p %>% 
      mutate(sig = p.value < thealpha,
             alpha = thealpha)
    
    df_simu_sig_omni <- df_simu_sig_long %>% 
      filter(effnames != "omniF") %>% 
      group_by(iter, N_IV, alpha) %>% 
      summarize(N_sig = sum(sig), 
                effnames = "sig_any",
                alpha = thealpha,
                sig = N_sig > 0, .groups="drop") %>% 
      bind_rows(df_simu_sig_long)
    
    df_simu_sig <- df_simu_sig_long %>% 
      filter(effnames == "omniF") %>% 
      select(-p.value) %>% 
      pivot_wider(names_from = effnames, values_from = sig) %>% 
      right_join(df_simu_sig_omni, by=c("iter", "N_IV", "alpha")) %>% 
      # filter(effnames != "omniF") %>% 
      mutate(sig_with_omni = omniF * sig)
    
    return(df_simu_sig)
  }
  
  # combine results with multiple alphas
  ltmp <- lapply(alphas, sig_omnibus_single, df_simu_p=df_simu_p)
  
  return(bind_rows(ltmp))
}

##### Main effect and post-hoc analysis #####
# Simulation for Type I error in main effect and post-hoc analysis
sim_main_posthoc <- function (N_subj = 30, iter = 100, n_core=2, file_cache = NULL,
                              N_con = 3, seed = 2022,
                              adjusts = c("none", "tukey", "scheffe", "sidak", "bonferroni", "dunnettx")
) {
  # N_subj: number of participants per condition (between-subject design)
  # iter: number of simulation/iteration
  # n_core: number of cores to be used for simulation
  # file_cache: file name of the cache file. If NULL, no file will be cached.
  # N_con: number of conditions for the one-way ANOVA
  # seed: seed used for simulation
  # adjusts: methods to be used in emmeans() to apply multiple comparison corrections
  
  # this function require library(emmeans)
  
  N_adjust <- length(adjusts)
  N_row <- iter * N_adjust
  set.seed(seed)
  
  # function to perform a single simulation
  sim_posthoc_single <- function(N_subj, N_con, adjusts) {
    
    # simulate data
    df_main_sim <- tibble(
      IV = rep(as.character(1:N_con), each=N_subj),
      DV = rnorm(N_subj * N_con),
      Subj = 1:(N_subj * N_con)
    )
    
    # perform ANOVA
    aov_main <- aov_4(DV ~ IV + (1|Subj), df_main_sim)
    # significance of the main effect
    p_main <- aov_main$Anova$`Pr(>F)`[2]
    
    # perform post-hoc (regardless of main effect significance)
    emm <- emmeans(aov_main, ~ IV)
    df_list <- vector(mode = "list", length = N_adjust)
    
    # this post hoc analysis results
    thisposthoc <- contrast(emm, "pairwise", adjust = "none")
    
    for (iadjust in 1:N_adjust) {
      
      # apply different multiple comparison corrections
      contra_sim <- summary(thisposthoc, adjust = adjusts[iadjust])
      
      df_contra_sim <- as_tibble(contra_sim) %>% 
        select(contrast, p.value) %>% 
        pivot_wider(names_from = contrast, names_prefix = "p_", values_from = p.value) %>% 
        mutate(adjust = adjusts[iadjust],
               p_main = p_main)
      
      df_list[[iadjust]] = df_contra_sim
    }
    
    one_main_posthoc <- bind_rows(df_list)
    
    return(one_main_posthoc)
  }
  
  # run simulation in parallel
  if (!is.null(file_cache) && file.exists(as.character(file_cache))){
    df_simu <- read_rds(file_cache)
    # message("Load simulation results from local cache files successfully.")
    
  } else {
    # set parallel processing
    Ns_iter <- rep(N_subj, times=iter)
    ls_tibble <- pbapply::pblapply(Ns_iter, sim_posthoc_single,  
                                   N_con=N_con, adjusts=adjusts, cl=n_core)
    df_simu <- bind_rows(ls_tibble, .id = "iter")
    
    if (!is.null(file_cache)){
      write_rds(df_simu, file=file_cache)
    }
  }
  
  return(df_simu)
}


sig_main_posthoc <- function(df_simu_p, alphas=0.05) {
  # df_simu_p: the output from p_main_posthoc()
  # alphas: alpha to be applied for claiming significant results
  
  # apply one alpha
  sig_main_posthoc_single <- function(df_simu_p, thealpha){
    
    df_simu_sig_long <- df_simu_p %>% 
      pivot_longer(contains(" - "), names_to = "contrast", values_to = "p.value") %>% 
      mutate(sig = p.value < thealpha,
             sig_main = p_main < thealpha,
             alpha = thealpha,
             contrast = str_remove(contrast, "p_")) 
    
    df_sim_sig_any <- df_simu_sig_long %>% 
      group_by(iter, adjust, alpha) %>% 
      summarize(N_sigpost = sum(sig), .groups = "drop") 
    
    df_simu_sig <- df_simu_sig_long  %>% 
      pivot_wider(id_cols = c(iter, adjust, alpha, sig_main), names_from = contrast, 
                  names_prefix = "sig_", values_from = sig) %>% 
      left_join(df_sim_sig_any, by=c("iter", "adjust", "alpha"))
    
    return(df_simu_sig) 
  }
  
  # combine results with multiple alphas
  ltmp <- lapply(alphas, sig_main_posthoc_single, df_simu_p=df_simu_p)
  
  return(bind_rows(ltmp))
}


##### Interaction and simple effect analysis #####
# Simulation for Type I error in main effect and post-hoc analysis
sim_inter_simple <- function (N_subj = 30, iter = 100, n_core=2, file_cache = NULL,
                              seed = 2022,
                              adjusts = c("none", "scheffe", "sidak", "bonferroni", "dunnettx")
) {
  # N_subj: number of participants per condition (between-subject design)
  # iter: number of simulation/iteration
  # n_core: number of cores to be used for simulation
  # file_cache: file name of the cache file. If NULL, no file will be cached.
  # seed: seed used for simulation
  # adjusts: methods to be used in emmeans() to apply multiple comparison corrections
  #           "tukey" was not used it is not applicable in this situation.
  
  # this function require library(emmeans)
  
  # apply a 2*2 (between-subject ) design
  N_adjust <- length(adjusts)
  N_row <- iter * N_adjust
  set.seed(seed)
  
  simu_simple_single <- function(N_subj, adjusts){
    
    df_null <- tibble(
      Subj = 1:(N_subj*4), 
      IV_A = rep(c("a1", "a2"), each = N_subj*2),
      IV_B = rep(c("b1", "b2"), times = N_subj*2),
      DV = rnorm(N_subj*4) 
    )
    
    aov_tmp <- aov_4(DV ~ IV_A * IV_B + (1|Subj), data = df_null)
    
    # interaction
    tmp_inter <- contrast(emmeans(aov_tmp, ~IV_A+IV_B),
                          interaction="pairwise") %>% 
      as_tibble()
    p_inter <- tmp_inter$p.value
    
    # main effects
    tmp_mainA <- contrast(emmeans(aov_tmp, ~IV_A), "pairwise") %>% 
      as_tibble()
    tmp_mainB <- contrast(emmeans(aov_tmp, ~IV_B), "pairwise") %>% 
      as_tibble()
    p_mainA <- tmp_mainA$p.value
    p_mainB <- tmp_mainB$p.value
    
    # simple effect (with/without multiple comparison corrections)
    tmp_simple_AB <- contrast(emmeans(aov_tmp, ~IV_A|IV_B), "pairwise", adjust="none")
    tmp_simple_BA <- contrast(emmeans(aov_tmp, ~IV_B|IV_A), "pairwise", adjust="none") 
    
    tmp_simple_all <- contrast(emmeans(aov_tmp, ~IV_A + IV_B), "pairwise", adjust="none")
    
    df_listAB <- vector(mode = "list", length = N_adjust)
    df_listBA <- vector(mode = "list", length = N_adjust)
    df_listall <- vector(mode = "list", length = N_adjust)
    
    for (iadjust in 1:N_adjust) {
      
      tmp_simple_AB_adjust <- summary(tmp_simple_AB[1:2], adjust = adjusts[iadjust])
      tmp_simple_BA_adjust <- summary(tmp_simple_BA[1:2], adjust = adjusts[iadjust])
      
      tmp_simple_all_adjust <- summary(tmp_simple_all[c(1:2,5:6)], adjust = adjusts[iadjust])
      
      p_simpleAB <- as_tibble(tmp_simple_AB_adjust) %>% 
        mutate(contrast = paste0(IV_B, "_(", contrast, ")")) %>% 
        select(contrast, p.value) %>% 
        pivot_wider(names_from = contrast, names_prefix = "p_", values_from = p.value) %>% 
        mutate(adjust = adjusts[iadjust])
      df_listAB[[iadjust]] = p_simpleAB
      
      p_simpleBA <- as_tibble(tmp_simple_BA_adjust) %>% 
        mutate(contrast = paste0(IV_A, "_(", contrast, ")")) %>% 
        select(contrast, p.value) %>% 
        pivot_wider(names_from = contrast, names_prefix = "p_", values_from = p.value) %>% 
        mutate(adjust = adjusts[iadjust])
      df_listBA[[iadjust]] = p_simpleBA
      
      p_simpleall <- as_tibble(tmp_simple_all_adjust) %>% 
        select(contrast, p.value) %>% 
        pivot_wider(names_from = contrast, names_prefix = "p_", values_from = p.value) %>% 
        mutate(adjust = adjusts[iadjust])
      df_listall[[iadjust]] = p_simpleall
    }
    
    # save all the p-values
    one_inter_simple <- left_join(bind_rows(df_listAB),
                                  bind_rows(df_listBA)) %>% 
      left_join(bind_rows(df_listall)) %>% 
      mutate(N = N_subj,
             p_inter = p_inter,
             p_mainA = p_mainA,
             p_mainB = p_mainB)
    
    return(one_inter_simple)
  }
  

  # run simulation in parallel
  if (!is.null(file_cache) && file.exists(as.character(file_cache))){
    df_simu <- read_rds(file_cache)
    message("Load simulation results from local cache files successfully.")
    
  } else {
    # set parallel processing
    Ns_iter <- rep(N_subj, times=iter)
    ls_tibble <- pbapply::pblapply(Ns_iter, simu_simple_single,  
                                   adjusts=adjusts, cl=n_core)
    df_simu <- bind_rows(ls_tibble, .id = "iter")
    
    if (!is.null(file_cache)){
      write_rds(df_simu, file=file_cache)
    }
  }
  
  return(df_simu)
}


sig_inter_simple <- function(df_simu_p, alphas=0.05) {
  # df_simu_p: the output from p_inter_simple()
  # alphas: alpha to be applied for claiming significant results
  
  # apply one alpha
  sig_main_posthoc_single <- function(df_simu_p, thealpha){
    
    df_simu_sig_long <- df_simu_p %>% 
      pivot_longer(contains(" - "), names_to = "contrast", values_to = "p.value") %>% 
      mutate(sig = p.value < thealpha,
             sig_inter = p_inter < thealpha,
             alpha = thealpha,
             contrast = str_remove(contrast, "p_")) 
    
    # whether any of the four simple effects are significant 
    # (when the multiple comparison corrections are applied for the four comparisons)
    df_sim_sig_4any <- df_simu_sig_long %>% 
      filter(substr(contrast, 3,3) != "_") %>% 
      group_by(iter, adjust, alpha) %>% 
      summarize(N_sigsimple4 = sum(sig), .groups = "drop") 
    
    # whether either simple effect is significant (corrected with two)
    df_sim_sig_2any <- df_simu_sig_long %>% 
      filter(substr(contrast, 3,3) == "_") %>% 
      mutate(IV_by = toupper(substr(contrast, 1, 1))) %>% 
      # separate(contrast, into = c("IV_by", "contrast"), sep = "_") %>% 
      group_by(iter, IV_by, adjust, alpha) %>% 
      summarize(N_sigsimple2 = sum(sig), .groups = "drop") %>% 
      pivot_wider(c(iter, adjust, alpha), names_from = IV_by, names_prefix = "N_sigsimple2by",
                  values_from = N_sigsimple2)
    
    df_simu_sig <- df_simu_sig_long  %>% 
      pivot_wider(id_cols = c(iter, adjust, alpha, sig_inter), names_from = contrast, 
                  names_prefix = "sig_", values_from = sig) %>% 
      left_join(df_sim_sig_4any, by=c("iter", "adjust", "alpha")) %>% 
      left_join(df_sim_sig_2any, by=c("iter", "adjust", "alpha"))
    
    return(df_simu_sig) 
  }
  
  # combine results with multiple alphas
  ltmp <- lapply(alphas, sig_main_posthoc_single, df_simu_p=df_simu_p)
  
  return(bind_rows(ltmp))
}
