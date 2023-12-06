library(tidyverse)
library(tidyplus)
library(comprehenr)
library(gridExtra)
library(cowplot)
library(stargazer)
library(fGarch)
library(parallel)
library(latex2exp)

# function: generating sample
# arguments: sample size, female feminist percentage, male feminist percentage, female percentage
#  deltas, own wealth weight mean and standard deviations, 
#  group wealth weight mean and standard deviations (sd), female wealth mean and sd,
#  male wealth mean and sd, female feminist weight mean and sd,
#  male feminist weight mean and sd.
gen_sample <- function(n, ffem_pct, mfem_pct, f_pct=0.5,
                       delta0, delta1, delta2, delta3,
                       oww_mean, oww_sd, gww_mean, gww_sd,
                       fw_mean, fw_sd, mw_mean, mw_sd,
                       ffem_w_mean, ffem_w_sd, mfem_w_mean, mfem_w_sd){
  tibble(agent = 1:n) %>% 
    mutate(delta0 = delta0, delta1 = delta1, 
           delta2 = delta2, delta3 = delta3,
           is_female = if_else(runif(n,min=0,max=1) < f_pct, T, F),
           is_feminist = if_else(
             is_female,
             if_else(runif(n,min=0,max=1) < ffem_pct, T, F),
             if_else(runif(n,min=0,max=1) < mfem_pct, T, F)
           ),
           own_wealth_weight = rnorm(n,mean = oww_mean, sd = oww_sd) %>%
             map(\(x) {
               min(x, 1) %>%
                 max(0)
             }) %>% unlist(),
           group_wealth_weight = rnorm(n,mean = gww_mean, sd = gww_sd) %>%
             map(\(x) {
               min(x, 1) %>%
                 max(0)
             }) %>% unlist(),
           # case_when(rnorm(n,mean = 0.5, sd = 1)<0 ~ 0,
           #                             rnorm(n,mean = 0.5, sd = 1)>1 ~ 1,
           #                             T~rnorm(n,mean = 0.5, sd = 1))
           agent = map(agent, \(x){str_c("a", as.character(x))}) %>% unlist()) %>% 
    mutate(wealth = if_else(is_female, rnorm(n, mean = fw_mean, sd = sqrt(fw_sd)), 
                            rnorm(n, mean = mw_mean, sd = sqrt(mw_sd))),
           feminist_weight = if_else(is_female, 
                                     rnorm(n,mean = ffem_w_mean, sd = ffem_w_sd), # increase order of magnitude to be comparable with weighted wealth
                                     rnorm(n,mean = mfem_w_mean, sd = mfem_w_sd)),
           own_wealth_weight_adj = if_else((own_wealth_weight + group_wealth_weight) == 0,
                                           0,
                                           own_wealth_weight/(own_wealth_weight + group_wealth_weight)),
           group_wealth_weight_adj =if_else((own_wealth_weight + group_wealth_weight) == 0,
                                            0,
                                            group_wealth_weight/(own_wealth_weight + group_wealth_weight)))
}

# function: generating encounter pairs
# arguments: meeting percentage (percentage of sample becoming the dictator), sample dataframe.
meet <- function(pct,df){ 
  num_actives = 0
  num = nrow(df)
  num_actives = round(num * pct)
  
  df_active <- df %>% 
    slice_sample(n = num_actives) %>%  
    mutate(ind = sample(seq(num_actives), num_actives),
           active_agent = agent)
  df_passive <- df %>% 
    anti_join(df_active, by = "agent") %>% 
    slice_sample(n = num_actives) %>% 
    mutate(ind = sample(seq(num_actives), num_actives),
           passive_agent = agent)
  
  df_encounters <- df_active %>% 
    inner_join(df_passive, by = "ind") %>% 
    select(active_agent, passive_agent)
  
  return(df_encounters)
}

# function: calculating identity distance (equation 1)
# arguments: female identity of dictator and counterpart (active and passive agents),
#  feminist identity of dictator and counterpart (active and passive agents),
#  deltas.
distance <- function(is_female_active_agent, is_female_passive_agent, 
                     is_feminist_active_agent, is_feminist_passive_agent,
                     delta0, delta1, delta2, delta3){
  same_h = if_else(is_female_active_agent == is_female_passive_agent, 1, 0)
  same_l = if_else(is_feminist_active_agent == is_feminist_passive_agent, 1, 0)
  
  dist = delta0 * same_h * same_l +
    delta1 * (1 - same_h)* (1 - same_l) +
    delta2 * same_h * (1 - same_l) +
    delta3 * same_l * (1 - same_h)
  
  return(dist)
}

# function: updating identity by wealth
# arguments: sample dataframe
update_sample <- function(df){
  
  max = max(df$wealth)
  min = min(df$wealth)
  
  df <- df %>% 
    mutate(norm_wealth = (wealth - min)/(max - min)) %>% 
    mutate(norm_wealth = norm_wealth*100) 
  
  wealth_df <- df %>% 
    group_by(is_feminist, is_female) %>% 
    summarise(mean_wealth = mean(norm_wealth)) %>% 
    ungroup()
  wealth_df <- expand_grid(tibble(is_female = c(T,F)),
                           tibble(is_feminist = c(T,F))) %>% 
    left_join(wealth_df) %>% 
    mutate(mean_wealth  = replace_na(mean_wealth, -1)) # if no member of the group exists, flag by -1 to be replaced by own wealth
  
  fem_F_w = wealth_df %>% filter(is_feminist,is_female) %>% select(mean_wealth) %>% unlist()
  fem_M_w = wealth_df %>% filter(is_feminist,!is_female) %>% select(mean_wealth) %>% unlist()
  non_fem_F_w = wealth_df %>% filter(!is_feminist,is_female) %>% select(mean_wealth) %>% unlist()
  non_fem_M_w = wealth_df %>% filter(!is_feminist,!is_female) %>% select(mean_wealth) %>% unlist()
  
  new_df <- df %>% 
    mutate(fem_wealth = if_else(is_female, fem_F_w, fem_M_w),
           non_fem_wealth = if_else(is_female, non_fem_F_w, non_fem_M_w),
           fem_wealth = if_else(fem_wealth <= -1, norm_wealth, fem_wealth),
           non_fem_wealth = if_else(non_fem_wealth <= -1, norm_wealth, non_fem_wealth),
           fem_component = feminist_weight * 1,
           u_fem = own_wealth_weight_adj * norm_wealth + 
             group_wealth_weight_adj * fem_wealth + 1 * tanh(fem_component), 
           u_non_fem = own_wealth_weight_adj * norm_wealth +     
             group_wealth_weight_adj * non_fem_wealth + fem_component * 0) %>% 
    mutate(is_feminist = u_fem > u_non_fem) %>% 
    select(-norm_wealth)

  return(new_df)
}

# function: updating feminist weight with network impacts
# arguments: sample dataframe, previous encouter history, log of previous dataframes
network_impact <- function(df, ect_log, dfs){
  df <- ect_log %>% 
    # group_by(active_agent, passive_agent, round) %>% 
    left_join(dfs %>% 
                select(passive_agent = agent, is_feminist, is_female, round),
              by = c("passive_agent", "round")) %>% 
    # ungroup() %>% 
    mutate(pa = T) %>% 
    group_by(active_agent) %>%
    summarise(fem_score = sum(is_feminist)/sum(pa)) %>% 
    arrange(fem_score) %>%
    rename("agent" = "active_agent") %>% 
    right_join(df, by = "agent") %>% 
    mutate(fem_score = replace_na(fem_score, 0)) %>% 
    mutate(feminist_weight = feminist_weight * (1 + fem_score)) 
  
  return(df)
}

delta0 = 1
delta1 = 5
delta2 = 7
delta3 = 3  

# function: "dictator game" giving
# arguments: sample dataframe, number of rounds, basic colummns needed in dataframe,
#  stake percentage, higher-outgroup coefficient, tanh coefficient
DG <- function(df, N, delta0, delta1, delta2, delta3, basecol,
               stake_pct, out_coef, tanh_coef){ 
  df <- df %>% 
    select(-contains("delta")) %>% 
    mutate(delta0 = delta0,
           delta1 = delta1,
           delta2 = delta2,
           delta3 = delta3) %>% 
    select(agent, contains("delta"), everything())
  
  encounter_log <- tibble()
  df_log <- df %>% mutate(round = 0, percent = -99)
  for(i in 1:N){
    pct <- runif(1, min = 0.3, max = 0.5) 
    encounter <- tibble()
    encounter_outcome <- tibble()
    encounter_types <- tibble()
    
    if(i > 1){
      df <- network_impact(df, encounter_log, df_log)
      df <- df %>%
        select(basecol)
    }
    encounter <- meet(pct, df)
    
    encounter_types <- encounter %>%
      pivot_longer(cols = everything(),
                   names_to = "type",
                   values_to = "agent") %>% 
      left_join(df %>% 
                  select(agent, is_female, is_feminist),
                by = "agent") %>%
      pivot_wider(names_from = "type",
                  values_from = c("agent","is_female","is_feminist")) %>% 
      unnest() %>% 
      left_join(df %>% 
                  select(agent_active_agent = agent, starts_with("delta")),
                by = "agent_active_agent")

    # behavioral outcomes
    if(i > 1){
      last_df <- df_log %>% 
        filter(round == i - 2) %>% 
        select(basecol)
      
      
      count1 <- df %>%
        group_by(is_female, is_feminist) %>% 
        summarise(count = n()) %>% 
        ungroup() 
      
      count2 <- last_df %>%
        group_by(is_female, is_feminist) %>% 
        summarise(count_last = n()) %>% 
        ungroup() 
      
      delta_mod <- count1 %>% 
        left_join(count2, by = c("is_female", "is_feminist")) %>% 
        right_join(expand_grid(is_female = c(T, F),
                               is_feminist = c(T, F)),
                   by = c("is_female", "is_feminist")) %>%
        replace_na(list("count" = 0, "count_last" = 0)) %>%
        transmute(is_female, is_feminist,
                  group_rate_multiplier = if_else(count_last != 0,
                                                  count / count_last,
                                                  1))
      
      encounter_types <- encounter_types %>% 
        left_join(delta_mod, by = c("is_female_active_agent" = "is_female", 
                                    "is_feminist_active_agent" = "is_feminist")) %>%
        mutate(delta3 = group_rate_multiplier/out_coef * delta3, 
               delta2 = group_rate_multiplier * delta2,
               delta1 = group_rate_multiplier/out_coef * delta1,
               delta0 = group_rate_multiplier * delta0) %>%  
        select(-group_rate_multiplier)
      
      df <- df %>% 
        left_join(delta_mod, by = c("is_female", "is_feminist")) %>%
        mutate(delta3 = group_rate_multiplier/out_coef * delta3, 
               delta2 = group_rate_multiplier * delta2,
               delta1 = group_rate_multiplier/out_coef * delta1,
               delta0 = group_rate_multiplier * delta0) %>%  
        select(-group_rate_multiplier)
    }
    
    encounter_outcome <- encounter_types %>%
      mutate(dist = pmap(
        encounter_types %>%
          select(starts_with("is"),
                 starts_with("delta")),
        .f = distance) %>% unlist()) %>% 
      left_join(df %>% select(agent, wealth), 
                by = c("agent_active_agent" = "agent")) %>% 
      transmute(active_agent = agent_active_agent,
                passive_agent = agent_passive_agent,
                dist,
                raw_giving_rate = rsnorm(nrow(encounter_types), mean = 30, sd = 25, xi = 1.5),
                raw_giving_rate = map(raw_giving_rate, \(x){max(x,0) %>% min(100)}) %>% unlist(), 
                stake_rate = runif(nrow(encounter_types), min = 0, max = stake_pct), 
                stake = wealth * stake_rate,
                raw_giving = raw_giving_rate * stake/100) %>% 
      mutate(giving = (1-tanh(dist/tanh_coef))*raw_giving, 
             keeping = 0 - giving) 
    
    # updating type
    df <- df %>% 
      left_join(encounter_outcome %>% 
                  transmute(agent = active_agent,
                            gain = keeping) %>% 
                  bind_rows(
                    encounter_outcome %>% 
                      transmute(agent = passive_agent,
                                gain = giving)
                  ),
                by = "agent") %>% 
      mutate(gain = replace_na(gain, 0)) %>% 
      mutate(wealth = wealth + gain) %>% 
      update_sample() 
    
    # logging df
    df_log <- df_log %>% 
      bind_rows(df %>% 
                  select(basecol) %>% 
                  mutate(round = i, percent = pct))
    
    # logging encounter
    encounter_log <- encounter_log %>% 
      bind_rows(encounter_outcome %>% 
                  mutate(round = i,
                         percent = pct))
    
    print("====end of round======")
  }
  return(list(latest_df = df,encounters = encounter_log, dfs = df_log)) 
} 
