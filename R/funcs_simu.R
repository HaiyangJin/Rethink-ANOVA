
# Omnibus ANOVA #####
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
    
    # simulate null data
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
    
    p.value <- aov_omni$Anova$`Pr(>F)`[c(-1, -length(aov_omni$Anova$`Pr(>F)`))]
    effnames <- row.names(aov_omni$Anova)[c(-1, -length(aov_omni$Anova$`Pr(>F)`))]
    
    p_df <- tibble(effnames = effnames,
                   p.value = p.value) %>% 
      mutate(N_IV = N_IV,
             omniF = p_omni)
    
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

sig_omnibus <- function(df_simu_p, alphas=0.05, isBonferroni=TRUE){
  
  sig_omnibus_single <- function(df_simu_p, thealpha, bonf){
    
    corr_str <- c("uncorrected", "Bonferroni")
    
    tmp_end <- max(df_simu_p$N_IV)
    if (bonf) {
      tmp_x <- 1
      num_effects <- c()
      for (i in 1:tmp_end) { 
        # 杨辉三角 (aka Pascal's triangle)
        tmp_x <- c(0, tmp_x) + c(tmp_x, 0)
        num_effects[i] <- sum(tmp_x[2:length(tmp_x)])
      }
      
    } else {
      num_effects <- rep(1, times=tmp_end)
    }
    
    df_simu_sig_long <- df_simu_p %>% 
      mutate(p.value = pmin(p.value * num_effects[N_IV], 1), # cap to 1
             sig = p.value < thealpha,
             sig_omniF = omniF < thealpha,
             alpha = thealpha,
             adjust = corr_str[bonf+1])
    
    df_simu_sig <- df_simu_sig_long %>% 
      group_by(iter, N_IV, alpha) %>% 
      summarize(sig_any = sum(sig) > 0, .groups="drop") %>% 
      left_join(df_simu_sig_long, by = c("iter", "N_IV", "alpha")) %>% 
      mutate(`sig_omni&any` = sig_omniF & sig_any) %>% 
      select(iter, N_IV, effnames, alpha, adjust, 
             p.value, sig, sig_any, omniF, sig_omniF, `sig_omni&any`)
    
    return(df_simu_sig)
  }
  
  # combine results with multiple alphas (without Bonferroni)
  ltmp <- lapply(alphas, sig_omnibus_single, 
                 df_simu_p=df_simu_p, bonf = FALSE)
  
  if (isBonferroni) {
    ltmpb <- lapply(alphas, sig_omnibus_single, 
                    df_simu_p=df_simu_p, bonf = TRUE)
    
    df_sig_omni <- bind_rows(ltmp, ltmpb)
  } else {
    df_sig_omni <- bind_rows(ltmp)
  }
  
  return(df_sig_omni)
}

# Main effect and post-hoc analysis #####
# Simulation for Type I error in main effect and post-hoc analysis
sim_main_posthoc <- function (N_subj = 30, iter = 100, n_core=2, file_cache = NULL,
                              N_levels = 3, seed = 2022,
                              adjusts = c("none", "tukey", "scheffe", "sidak", "bonferroni", "dunnettx")
) {
  # N_subj: number of participants per condition (between-subject design)
  # iter: number of simulation/iteration
  # n_core: number of cores to be used for simulation
  # file_cache: file name of the cache file. If NULL, no file will be cached.
  # N_levels: number of conditions for the one-way ANOVA
  # seed: seed used for simulation
  # adjusts: methods to be used in emmeans() to apply multiple comparison corrections
  
  # this function require library(emmeans)
  
  N_adjust <- length(adjusts)
  N_row <- iter * N_adjust
  set.seed(seed)
  
  # function to perform a single simulation
  sim_posthoc_single <- function(N_subj, N_level, adjusts) {
    
    # simulate data
    df_main_sim <- tibble(
      IV = rep(as.character(1:N_level), each=N_subj),
      DV = rnorm(N_subj * N_level),
      Subj = 1:(N_subj * N_level)
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
        mutate(adjust = adjusts[iadjust],
               p_main = p_main,
               N_level = N_level)
      
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
    Ns_iter <- rep(N_levels, times=iter)
    ls_tibble <- pbapply::pblapply(Ns_iter, sim_posthoc_single,  
                                   N_subj=N_subj, adjusts=adjusts, cl=n_core)
    df_simu <- bind_rows(ls_tibble, .id = "iter") %>% 
      mutate(iter = as_factor(iter),
             contrast = as_factor(contrast),
             adjust = if_else(adjust=="none", "uncorrected", adjust),
             adjust = as_factor(adjust))
    
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
      mutate(sig = p.value < thealpha,
             sig_main = p_main < thealpha,
             alpha = thealpha,
             contrast = str_remove(contrast, "p_")) 
    
    df_simu_sig <- df_simu_sig_long %>% 
      group_by(N_level, iter, adjust, alpha) %>% 
      summarize(sig_anypost = sum(sig) > 0,
                .groups = "drop") %>% 
      left_join(df_simu_sig_long, 
                by = c("N_level", "iter", "adjust", "alpha")) %>% 
      mutate(`sig_main&posthoc` = sig_main & sig_anypost) %>% 
      select(iter, N_level, contrast, alpha, adjust, 
             p.value, sig, sig_anypost, p_main, everything())
    
    return(df_simu_sig) 
  }
  
  # combine results with multiple alphas
  ltmp <- lapply(alphas, sig_main_posthoc_single, df_simu_p=df_simu_p)
  
  return(bind_rows(ltmp))
}


# Interaction and simple effect analysis #####
# Simulation for Type I error in interaction and simple effect analysis
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
      
      # a group of two
      tmp_simple_AB_adjust <- summary(tmp_simple_AB[1:2], adjust = adjusts[iadjust])
      # a group of two
      tmp_simple_BA_adjust <- summary(tmp_simple_BA[1:2], adjust = adjusts[iadjust])
      
      # a group of 4 tests
      tmp_simple_all_adjust <- summary(tmp_simple_all[c(1:2,5:6)], adjust = adjusts[iadjust])
      
      p_simpleAB <- as_tibble(tmp_simple_AB_adjust) %>% 
        mutate(contrast = paste0(IV_B, "_(", contrast, ")")) %>% 
        select(contrast, p.value) %>% 
        mutate(adjust = adjusts[iadjust],
               group = "byB_2")
      df_listAB[[iadjust]] = p_simpleAB
      
      p_simpleBA <- as_tibble(tmp_simple_BA_adjust) %>% 
        mutate(contrast = paste0(IV_A, "_(", contrast, ")")) %>% 
        select(contrast, p.value) %>% 
        mutate(adjust = adjusts[iadjust],
               group = "byA_2")
      df_listBA[[iadjust]] = p_simpleBA
      
      p_simpleall <- as_tibble(tmp_simple_all_adjust) %>% 
        select(contrast, p.value) %>% 
        mutate(adjust = adjusts[iadjust],
               group = "byAB_4")
      df_listall[[iadjust]] = p_simpleall
    }
    
    # save all the p-values
    one_inter_simple <- bind_rows(df_listAB, df_listBA, df_listall) %>% 
      mutate(N = N_subj, p_inter = p_inter,
             p_mainA, p_mainB) # add other information
    
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
    df_simu <- bind_rows(ls_tibble, .id = "iter") %>% 
      mutate(iter = as_factor(iter),
             contrast = as_factor(contrast),
             group = as_factor(group),
             adjust = if_else(adjust=="none", "uncorrected", adjust),
             adjust = as_factor(adjust))
    
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
  sig_inter_simple_single <- function(df_simu_p, thealpha){
    
    df_simu_sig_long <- df_simu_p %>% 
      mutate(sig_inter = p_inter < thealpha,
             sig_simple = p.value < thealpha,
             alpha = thealpha,
             contrast = str_remove(contrast, "p_")) 
    
    # whether any of the four/two simple effects are significant 
    # (when the multiple comparison corrections are applied for the four/two comparisons)
    df_sim_sig_any <- df_simu_sig_long %>% 
      group_by(group, iter, adjust, alpha) %>% 
      summarize(sig_any_simple = sum(sig_simple)>0, .groups = "drop")
    
    df_simu_sig <- df_simu_sig_long  %>% 
      select(iter, adjust, group, alpha, N, contrast, sig_simple, sig_inter) %>% 
      inner_join(df_sim_sig_any, by=c("iter", "adjust", "group", "alpha")) %>% 
      mutate(`sig_inter&any` = sig_inter & sig_any_simple)
    
    return(df_simu_sig) 
  }
  
  # combine results with multiple alphas
  ltmp <- lapply(alphas, sig_inter_simple_single, df_simu_p=df_simu_p)
  
  return(bind_rows(ltmp))
}



# Interaction and simple interaction analysis #####
# Simulation for Type I error in interaction and simple effect analysis
sim_simpinter <- function(N_subj = 30, iter = 100, n_core=2, file_cache = NULL,
                             seed = 2022, nlevel_A = 2, nlevel_B = 3,
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
  
  # apply a 2*3 (between-subject ) design
  set.seed(seed)
  
  levels_A <- paste0("a", 1:nlevel_A) 
  levels_B <- paste0("b", 1:nlevel_B)
  
  simu_simpinter_single <- function(N_subj, adjusts){
    
    N_adjust <- length(adjusts)
    
    df_null <- tibble(
      Subj = 1:(N_subj*length(levels_A)*length(levels_B)), 
      IV_A = rep(levels_A, each = N_subj*length(levels_B)),
      IV_B = rep(levels_B, times = N_subj*length(levels_A)),
      DV = rnorm(N_subj*length(levels_A)*length(levels_B)) 
    )
    
    aov_tmp <- aov_4(DV ~ IV_A * IV_B + (1|Subj), data = df_null)
    
    # p value of the 2*3 interaction
    p_inter <- aov_tmp$anova_table$`Pr(>F)`[3]
    
    # simple interaction analysis (with/without multiple comparison corrections)
    simpinter <- contrast(emmeans(aov_tmp, ~IV_A+IV_B),
                             interaction="pairwise", adjust="none") 
    
    df_listsim <- vector(mode = "list", length = N_adjust)
    
    for (iadjust in 1:N_adjust) {
      
      # simple interaction of all
      tmp_simpinter_adjust <- summary(simpinter, adjust = adjusts[iadjust])
      
      p_simpleinter <- as_tibble(tmp_simpinter_adjust) %>% 
        select(IV_A_pairwise, IV_B_pairwise, p.value) %>% 
        mutate(adjust = adjusts[iadjust])
      
      df_listsim[[iadjust]] = p_simpleinter
    }
    
    # save all the p-values
    one_simpinter <- df_listsim %>% 
      bind_rows() %>% 
      mutate(N = N_subj, p_inter = p_inter) # add other information
    
    return(one_simpinter)
  }
  
  
  # run simulation in parallel
  if (!is.null(file_cache) && file.exists(as.character(file_cache))){
    df_simu <- read_rds(file_cache)
    message("Load simulation results from local cache files successfully.")
    
  } else {
    # set parallel processing
    Ns_iter <- rep(N_subj, times=iter)
    ls_tibble <- pbapply::pblapply(Ns_iter, simu_simpinter_single,  
                                   adjusts=adjusts, cl=n_core)
    df_simu <- bind_rows(ls_tibble, .id = "iter") %>% 
      mutate(iter = as_factor(iter),
             IV_A_pairwise = as_factor(IV_A_pairwise),
             IV_B_pairwise = as_factor(IV_B_pairwise),
             adjust = if_else(adjust=="none", "uncorrected", adjust),
             adjust = as_factor(adjust))
    
    if (!is.null(file_cache)){
      write_rds(df_simu, file=file_cache)
    }
  }
  
  return(df_simu)
}


sig_simpinter <- function(df_simu_p, alphas=0.05) {
  # df_simu_p: the output from sim_simpinter()
  # alphas: alpha to be applied for claiming significant results
  
  # apply one alpha
  sig_simpinter_single <- function(df_simu_p, thealpha){
    
    df_simu_sig_long <- df_simu_p %>% 
      mutate(sig_inter = p_inter < thealpha,
             sig_simpinter = p.value < thealpha,
             alpha = thealpha,
             A = str_remove_all(IV_A_pairwise, " "),
             B = str_remove_all(IV_B_pairwise, " "),
             A = factor(A),
             B = factor(B)) %>% 
      select(-c(IV_A_pairwise, IV_B_pairwise))
    
    # whether any of the four/two simple effects are significant 
    # (when the multiple comparison corrections are applied for the four/two comparisons)
    df_siminter_sig_any <- df_simu_sig_long %>% 
      group_by(iter, adjust, alpha) %>% 
      summarize(sig_any_simpinter = sum(sig_simpinter)>0, .groups = "drop")
    
    df_simu_sig <- df_simu_sig_long  %>% 
      select(iter, adjust, alpha, N, A, B, sig_simpinter, sig_inter) %>% 
      inner_join(df_siminter_sig_any, by=c("iter", "adjust", "alpha")) %>% 
      mutate(`sig_inter&any` = sig_inter & sig_any_simpinter)
    
    return(df_simu_sig) 
  }
  
  # combine results with multiple alphas
  ltmp <- lapply(alphas, sig_simpinter_single, df_simu_p=df_simu_p)
  
  return(bind_rows(ltmp))
}




