library(dplyr)
setwd("~/Documents/GitHub/paralytic_polio_estimates")

source("simulation_functions.R")


nsims <- 1000

i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(i)

final <- random_simulation(n=nsims,index_start=(i-1)*nsims + 1,
                               incu_mean_prior_mean=16,incu_mean_prior_var=5,
                               incu_var_prior_mean=10,incu_var_prior_var=3,
                               infect_mean_prior_mean=7,infect_mean_prior_var=5,
                               infect_var_prior_mean=5,infect_var_prior_var=3,
                               tmax=180,
                               prop_immune_mean = 0.8,prop_immune_var=0.05,
                               prob_paralysis_var = 1/2000/10,
                               ini_infs=10,
                               R0_min=0,R0_max=10)


trajectories <- final[[1]]
pars <- final[[3]]
final_sizes <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),final_size=final[[4]])
paralysis_totals <- final[[7]]
data_consistent <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),consistent=final[[8]])
paralysis_totals <- paralysis_totals %>% left_join(data_consistent) %>% filter(consistent==TRUE)

## Flag which runs are consistent with the data
use_samps <- paralysis_totals %>% filter(consistent == TRUE) %>% pull(sim)

## Get outbreak length
run_times <- trajectories %>% filter(sim %in% use_samps) %>% 
    group_by(sim) %>% filter(t == max(t)) %>%
    rename(tmax = t) %>%
    rename(inc_end = inc) %>%
    mutate(ongoing = inc_end>0)

## Any infections in the X days prior to the end
x <- 10
run_times_long <- trajectories %>% filter(sim %in% use_samps) %>% 
    group_by(sim) %>% filter(t == max(t)) %>%
    rename(tmax = t) %>%
    select(tmax, sim) %>%
    left_join(trajectories) %>%
    group_by(sim) %>%
    filter(t >= (tmax - x)) %>%
    summarize(sum_final_infs = sum(inc)) %>%
    mutate(ongoing_10= sum_final_infs != 0)

## Get parameters
pars <- pars %>% 
    filter(sim %in% use_samps) %>% 
    mutate(Re=R0*prop_immune) %>%
    left_join(run_times) %>% 
    left_join(run_times_long) %>%
    left_join(final_sizes)

save(pars, file=paste0("sims/simulation_",i,".RData"))

