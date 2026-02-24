library(tidyverse)

setting <- 1
n_trials <- 5
n_estimators <- 3 # naive, standardized, ipw
n_sim = 1000
#epsilon <- .0000001
#setwd(paste("/Volumes/homes/setting",setting,"_nt",n_t,sep=""))

#wd <- paste0("home/srosin/sim_code/results/setting", setting, "_ntrials", n_trials)
outdir <- paste0("H:/sim_code/results/setting", setting, "_ntrials", n_trials)
wd <- paste0("H:/sim_code/results/setting", setting, "_ntrials", n_trials, "/raw_results/")
setwd(wd)

length(list.files())
#n_sim = 5000
#n_sim = length(list.files())

truevalues <- read.csv(paste0("H:/sim_code/truevalues/setting", 
                              setting, 
                              "_ntrials", 
                              n_trials,
                              "/truevalues.csv"))

truth_rho_p <- truevalues$rho_p
truth_beta <- truevalues$beta
truth_rho_s <- truevalues$rho_s

dscap_est_files <- list.files()[str_detect(list.files(), "cor_est")][1:n_sim]
vcov_files <- list.files()[str_detect(list.files(), "vcov")][1:n_sim]
lrt_files <- list.files()[str_detect(list.files(), "lrt")][1:n_sim]
trt_effect_files <- list.files()[str_detect(list.files(), "trt_effect")][1:n_sim]
bootstrap_files <- list.files()[str_detect(list.files(), "bootstrap")][1:n_sim]

# bias --------------------------------------------------------------------

dscap_est_list <- vector(mode= "list", length= length(dscap_est_files))
for(i in 1:length(dscap_est_files)){
  dscap_est_df <- read.csv(dscap_est_files[i])
  dscap_est_df$iteration <- i
  dscap_est_list[[i]] <- dscap_est_df
}

dscap_ests <- do.call(bind_rows, dscap_est_list)
names(dscap_ests) <- c("rho_p", "rho_s", "beta", "type", "iteration")

bias_df <- pivot_longer(dscap_ests,
                        cols = c(rho_p, rho_s, beta),
                        names_to = "parameter")|> 
  rowwise() |> 
  mutate(truevalue = truevalues[,parameter],
         bias = value - truevalues[,parameter])

mean_bias_df <- bias_df |> 
  group_by(type, parameter) |> 
  dplyr::summarize(mean_estimate = mean(value),
                   truevalue = mean(truevalue),
                   mean_bias = mean(bias),
                   median_bias = median(bias),
                   empirical_se = sd(value),
                   .groups = 'drop') |> 
  arrange(parameter, type)

mean_bias_df |> 
  select(-c(mean_estimate, truevalue)) 

write.csv(mean_bias_df, paste0(outdir, "/bias.csv"), row.names=F)

# individual treatment effects --------------------------------------------

trt_effect_list <- vector(mode= "list", length= length(trt_effect_files))
for(i in 1:length(trt_effect_files)){
  trt_effect_df <- read.csv(trt_effect_files[i])
  trt_effect_df$iteration <- i
  trt_effect_list[[i]] <- trt_effect_df
}

trt_effects <- do.call(bind_rows, trt_effect_list)[,c(1:8,20)]


## naive estimator -------------------------------------------------------------------

potential_outcome_mean_ests_naive_Anot0 <- trt_effects |> 
  filter(type == "naive") |> 
  select(-c(VE_est, mean_diff_S_est, mean_Y_0, mean_S_0)) |> 
  rename(A = trial, 
         Y_a = mean_Y_1, 
         S_a = mean_S_1) |> 
  mutate(trial = A)

potential_outcome_mean_ests_naive_A0 <- trt_effects |> 
  filter(type == "naive") |> 
  select(-c(VE_est, mean_diff_S_est, mean_Y_1, mean_S_1)) |> 
  mutate(A = 0) |> 
  rename(Y_a = mean_Y_0, 
         S_a = mean_S_0)

potential_outcome_mean_ests_naive <- bind_rows(
  potential_outcome_mean_ests_naive_A0,
  potential_outcome_mean_ests_naive_Anot0
)

trt_effect_bias_df_naive <- pivot_longer(potential_outcome_mean_ests_naive,
                                         cols = c(Y_a, S_a),
                                         names_to = "parameter")|> 
  rename(estimate = value) |> 
  mutate(truevalue = case_when(
    A == 0 & parameter == "Y_a" ~ truevalues$ey0,
    A == 1 & parameter == "Y_a" ~ truevalues$ey1,
    A == 2 & parameter == "Y_a" ~ truevalues$ey2,
    A == 3 & parameter == "Y_a" ~ truevalues$ey3,
    A == 4 & parameter == "Y_a" ~ truevalues$ey4,
    A == 5 & parameter == "Y_a" ~ truevalues$ey5,
    A == 0 & parameter == "S_a" ~ truevalues$es0,
    A == 1 & parameter == "S_a" ~ truevalues$es1,
    A == 2 & parameter == "S_a" ~ truevalues$es2,
    A == 3 & parameter == "S_a" ~ truevalues$es3,
    A == 4 & parameter == "S_a" ~ truevalues$es4,
    A == 5 & parameter == "S_a" ~ truevalues$es5
  ), 
  bias = (estimate - truevalue) / truevalue)

mean_trt_effect_bias_df_naive <- trt_effect_bias_df_naive |> 
  group_by(parameter, A, type) |> 
  dplyr::summarize(mean_estimate = mean(estimate),
                   median_estimate = median(estimate),
                   truevalue = mean(truevalue),
                   mean_relative_bias = round(mean(bias), 3),
                   empirical_se = round(sd(estimate), 3),
                   .groups = 'drop') |> 
  mutate(type = factor(type, levels = c("naive", "ipw", "standardized"))) |> 
  arrange(parameter, A, type)

View(mean_trt_effect_bias_df_naive |> filter(A != 1))
write.csv(mean_trt_effect_bias_df_naive, "po_means_bias_ese_naive.csv", row.names=F)


## nonnaive ----------------------------------------------------------------
potential_outcome_mean_ests_nonplacebo_nonnaive <- trt_effects |> 
  select(-c(VE_est, mean_diff_S_est, mean_Y_0, mean_S_0)) |> 
  rename(A = trial, 
         Y_a = mean_Y_1, 
         S_a = mean_S_1)

potential_outcome_mean_ests_placebo_nonnaive <- trt_effects |>
  filter(type != "naive") |> 
  #mutate(sim = rep(1:n_sim, each = n_trials * n_estimators)) |> 
  group_by(iteration, type) |> 
  dplyr::summarize(Y_a = mean(mean_Y_0), # the estimates of E[Y0] and E[S0] are the same across trial-rows in the dataframe,
                   S_a = mean(mean_S_0),        # so could take mean, or min, or max
                   .groups = 'drop') |> 
  mutate(A = 0) #|> 
#select(-sim)

potential_outcome_mean_ests_nonnaive <- bind_rows(
  potential_outcome_mean_ests_nonplacebo_nonnaive, 
  potential_outcome_mean_ests_placebo_nonnaive)

trt_effect_bias_df <- pivot_longer(potential_outcome_mean_ests_nonnaive,
                                   cols = c(Y_a, S_a),
                                   names_to = "parameter")|> 
  rename(estimate = value) |> 
  mutate(truevalue = case_when(
    A == 0 & parameter == "Y_a" ~ truevalues$ey0,
    A == 1 & parameter == "Y_a" ~ truevalues$ey1,
    A == 2 & parameter == "Y_a" ~ truevalues$ey2,
    A == 3 & parameter == "Y_a" ~ truevalues$ey3,
    A == 4 & parameter == "Y_a" ~ truevalues$ey4,
    A == 5 & parameter == "Y_a" ~ truevalues$ey5,
    A == 0 & parameter == "S_a" ~ truevalues$es0,
    A == 1 & parameter == "S_a" ~ truevalues$es1,
    A == 2 & parameter == "S_a" ~ truevalues$es2,
    A == 3 & parameter == "S_a" ~ truevalues$es3,
    A == 4 & parameter == "S_a" ~ truevalues$es4,
    A == 5 & parameter == "S_a" ~ truevalues$es5
  ), 
  bias = (estimate - truevalue) / truevalue)

mean_trt_effect_bias_df <- trt_effect_bias_df |> 
  group_by(parameter, A, type) |> 
  dplyr::summarize(mean_estimate = mean(estimate),
                   median_estimate = median(estimate),
                   truevalue = mean(truevalue),
                   mean_relative_bias = round(mean(bias), 3),
                   empirical_se = round(sd(estimate), 3),
                   .groups = 'drop') |> 
  mutate(type = factor(type, levels = c("naive", "ipw", "standardized"))) |> 
  arrange(parameter, A, type)

View(mean_trt_effect_bias_df |> filter(A != 1))
write.csv(mean_trt_effect_bias_df, "po_means_bias_ese.csv", row.names=F)

# taus and deltas ---------------------------------------------------------

taus_deltas <- trt_effects |> 
  select(c(trial, VE_est, mean_diff_S_est, type, iteration)) |> 
  rename(A = trial, 
         tau_a = VE_est, 
         delta_a = mean_diff_S_est)

tau_delta_bias_df <- pivot_longer(taus_deltas,
                                  cols = c(tau_a, delta_a),
                                  names_to = "parameter")|> 
  rename(estimate = value) |> 
  mutate(truevalue = case_when(
    A == 1 & parameter == "tau_a" ~ truevalues$tau.1.,
    A == 2 & parameter == "tau_a" ~ truevalues$tau.2.,
    A == 3 & parameter == "tau_a" ~ truevalues$tau.3.,
    A == 4 & parameter == "tau_a" ~ truevalues$tau.4.,
    A == 5 & parameter == "tau_a" ~ truevalues$tau.5.,
    A == 1 & parameter == "delta_a" ~ truevalues$delta.1.,
    A == 2 & parameter == "delta_a" ~ truevalues$delta.2.,
    A == 3 & parameter == "delta_a" ~ truevalues$delta.3.,
    A == 4 & parameter == "delta_a" ~ truevalues$delta.4.,
    A == 5 & parameter == "delta_a" ~ truevalues$delta.5.
  ), 
  bias = (estimate - truevalue) / truevalue)

mean_tau_delta_bias_df <- tau_delta_bias_df |> 
  group_by(parameter, A, type) |> 
  dplyr::summarize(mean_estimate = mean(estimate),
                   median_estimate = median(estimate),
                   truevalue = mean(truevalue),
                   mean_relative_bias = round(mean(bias), 3),
                   empirical_se = round(sd(estimate), 3),
                   .groups = 'drop') |> 
  mutate(type = factor(type, levels = c("naive", "ipw", "standardized"))) |> 
  arrange(parameter, A, type)

