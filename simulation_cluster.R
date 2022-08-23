library(dplyr)
setwd("~/paralytic_polio_estimates")
#setwd("~/Documents/GitHub/paralytic_polio_estimates")

source("simulation_functions.R")


nsims <- 1000

i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(i)
#i <- 1
set.seed(i)
res <- random_simulation(n=nsims,index_start=(i-1)*nsims + 1,
                               incu_mean_prior_mean=16,incu_mean_prior_var=5,
                               incu_var_prior_mean=10,incu_var_prior_var=3,
                               infect_mean_prior_mean=7,infect_mean_prior_var=5,
                               infect_var_prior_mean=5,infect_var_prior_var=3,
                               tmax=180,
                               prop_immune_mean = 0.8,prop_immune_var=0.05,
                               prob_paralysis_var = 1/2000/10,
                               ini_infs=10,
                               R0_min=0,R0_max=10)


trajectories <- res$res
pars <- res$pars
final_sizes <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),
                          final_size=res$final_size)
paralysis_totals <- res$n_paralysis
data_consistent <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),
                              consistent=res$data_are_consistent)
paralysis_totals <- paralysis_totals %>% left_join(data_consistent) %>% filter(consistent==TRUE)

## Flag which runs are consistent with the data
use_samps <- paralysis_totals %>% filter(consistent == TRUE) %>% pull(sim)
print(use_samps)
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
if(length(use_samps) > 1){
    res2 <- restart_simulations_table(use_sims = use_samps,
                              pars=res$pars,
                              susceptibles=res$susceptibles,
                              paralysis=res$paralysis,
                              incidence=res$incidence,
                              t_starts=res$tmax_vector,
                              tmax=365,nruns=1)
    traj_future <- res2$res
    traj_future_par <- res2$res_par
    save(traj_future, file=paste0("sims_traj/traj_",i,".RData"))
    save(traj_future_par, file=paste0("sims_para/para_",i,".RData"))
}

#traj_future_par %>% ggplot() + geom_line(aes(x=t,y=para,group=interaction(sim,rep)))
