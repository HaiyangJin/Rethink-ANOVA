


### Simulation for main effect and post-hoc analysis
sim_main_posthoc <- function (N_subj = 30, N_sim = 1000, alpha = 0.05, adjusts = c("none", "tukey", "scheffe", "sidak", "bonferroni", "dunnettx"), seed = 2022) {
  
  N_con <- 3
  N_adjust <- length(adjusts)
  N_row <- N_sim * N_adjust
  # N_con: number of conditions
  # N_subj: number of participants per condition (between-subject design)
  # N_sim: number of simulation
  # alpha: alpha to control Type I error rate
  # seed: seed used for simulation
  
  # this function require library(emmeans)
  
  set.seed(seed)
  # NaN vector to save simulation results
  p_main <- rep(NaN, N_sim)
  p_12 <- rep(NaN, N_row)
  p_13 <- rep(NaN, N_row)
  p_23 <- rep(NaN, N_row)
  sig_main <- rep(NaN, N_sim)
  sig_12 <- rep(NaN, N_row)
  sig_13 <- rep(NaN, N_row)
  sig_23 <- rep(NaN, N_row)
  iteration_sim <- rep(NaN, N_row)
  adjust_post <- rep('', N_row)
  
  # Each simulation
  for (i in 1:N_sim) {
    # simulate data
    df_main_sim <- tibble(
      IV = rep(as.character(1:N_con), each=N_subj),
      DV = rnorm(N_subj * N_con)
    )
    
    # perform ANOVA
    aov_main <- aov(DV ~ IV, df_main_sim)
    main_sum <- summary(aov_main)
    # significance of the main effect
    p_main[i] <- main_sum[[1]][["Pr(>F)"]][1]
    sig_main[i] <- main_sum[[1]][["Pr(>F)"]][1] < alpha
    
    # perform post-hoc (regardless of main effect significance)
    emm <- emmeans(aov_main, ~ IV)
    
    for (iadjust in 1:N_adjust) {
      
      contra_sim <- contrast(emm, "pairwise", adjust = adjusts[iadjust])
      
      df_contra_sim <- as_tibble(contra_sim)
      
      post_index <- (i-1) * N_adjust + iadjust
      
      iteration_sim[post_index] <- i
      adjust_post[post_index] <- adjusts[iadjust]
      p_12[post_index] <- df_contra_sim$p.value[1]
      p_13[post_index] <- df_contra_sim$p.value[2]
      p_23[post_index] <- df_contra_sim$p.value[3]
      sig_12[post_index] <- df_contra_sim$p.value[1] < alpha
      sig_13[post_index] <- df_contra_sim$p.value[2] < alpha
      sig_23[post_index] <- df_contra_sim$p.value[3] < alpha
      
    }
    
  }
  
  out <- merge(tibble(p_main, sig_main,
                      iteration_sim = 1:N_sim),
               tibble(iteration_sim, p_12, p_13, p_23, 
                      sig_12, sig_13, sig_23, adjust_post))
  return(out)
}
