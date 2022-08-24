library(dplyr)
setwd("~/paralytic_polio_estimates")
#setwd("~/Documents/GitHub/paralytic_polio_estimates")

source("simulation_functions.R")


nsims <- 10000

i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#i <- 1
print(i)
set.seed(i)
res <- random_simulation_twoimmune(n=nsims,index_start=(i-1)*nsims + 1,
                               incu_mean_prior_mean=16,incu_mean_prior_var=5,
                               incu_var_prior_mean=10,incu_var_prior_var=3,
                               infect_mean_prior_mean=7,infect_mean_prior_var=5,
                               infect_var_prior_mean=5,infect_var_prior_var=3,
                               tmax=180, P=350000,
                               prob_paralysis_var = 1e-8,
                               ini_infs=10,
                               R0_par1=0,R0_par2=10,
                               prop_immune_pars = c(15,10,2),
                                rel_R0_mean = 0.75)


trajectories <- res$simulation_results
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

trajectories <- trajectories %>% filter(sim %in% use_samps)

## Get outbreak length
run_times <- trajectories %>% filter(sim %in% use_samps) %>% 
    group_by(sim) %>% filter(t == max(t)) %>%
    rename(tmax = t) %>%
    rename(inc_end = inc) %>%
    mutate(ongoing = inc_end>0)

t_starts <- trajectories %>% filter(para == 1) %>% group_by(sim) %>% filter(t == min(t)) %>% rename(tstart=t) %>% select(sim, tstart)

## Any infections in the X days prior to the end
x <- 10
run_times_long <- trajectories %>% filter(sim %in% use_samps) %>% 
    group_by(sim) %>% filter(t == max(t)) %>%
    rename(tmax = t) %>%
    select(tmax, sim) %>%
    left_join(trajectories) %>%
    group_by(sim) %>%
    filter(t >= (tmax - x)) %>%
    dplyr::summarize(sum_final_infs = sum(inc),
                     sum_final_para = sum(para)) %>%
    mutate(ongoing_10= sum_final_infs != 0,
           ongoing_10_para=sum_final_para != 0)

## Get parameters
pars <- pars %>% 
    filter(sim %in% use_samps) %>% 
    mutate(Re=R0*prop_immune_groups.1) %>%
    left_join(run_times) %>% 
    left_join(run_times_long) %>%
    left_join(final_sizes) %>%
    left_join(t_starts) %>%
    mutate(date_start = as.Date("2022-07-18")-tstart)

save(pars, file=paste0("sims/simulation_",i,".RData"))

## Restart simulations where they left off with/without intervention strategy
if(length(use_samps) > 1){
    vacc_strats <- matrix(c(0,0,
                            0.5,0.5),ncol=2,byrow=TRUE)
    res2 <- restart_simulations_table_twoimmune(use_sims = use_samps,
                              pars=res$pars,
                              final_conditions = res$final_conditions,
                              t_starts=res$tmax_vector,
                              vaccinate_proportion = vacc_strats,
                              P=350000,
                              tmax=500,nruns=1)
    save(res2, file=paste0("sims_traj/traj_",i,".RData"))
    
}

#traj_future_par %>% ggplot() + geom_line(aes(x=t,y=para,group=interaction(sim,rep)))
