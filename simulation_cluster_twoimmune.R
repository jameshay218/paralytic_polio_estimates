library(dplyr)
setwd("~/paralytic_polio_estimates")
#setwd("~/Documents/GitHub/paralytic_polio_estimates")

source("simulation_functions_twoimmune.R")

nsims <- 50000
#nsims <- 1000

i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#i <- 1
print(i)
set.seed(i)

res <- random_simulation_twoimmune(n=nsims,
                                   observed_data=c(1, rep(0,48)),
                                   
                                   index_start=(i-1)*nsims + 1,
                               incu_mean_prior_mean=14,
                               incu_mean_prior_var=3,
                               incu_var_prior_mean=15,
                               incu_var_prior_var=1,
                               gen_interval_susc_shape_par1=3.87,
                               gen_interval_susc_shape_par2=0.52,
                               gen_interval_susc_rate_par1=23.86,
                               gen_interval_susc_rate_par2=52.71,
                               
                               gen_interval_partial_shape_par1=2.60,
                               gen_interval_partial_shape_par2=0.61,
                               gen_interval_partial_rate_par1=6.82,
                               gen_interval_partial_rate_par2=19.00,
                               
                               tmax=180, P=325000,
                               prob_paralysis_mean=0.0005,
                               prob_paralysis_ps_par1 = -4,
                               prob_paralysis_ps_par2 = -1,
                               prob_paralysis_var = 2e-8,
                               #prob_paralysis_ps_var = 1e-11,
                               ini_infs=1,
                               R0_dist="truncnorm",
                               R0_par1=4.9,R0_par2=2,
                               rel_R0_par1=7.04,
                               rel_R0_par2=65.16,
                               prop_susceptible_par1 = 39.34,
                               prop_susceptible_par2 = 144.52,
                               prop_refractory_par1 = 4.66,
                               prop_refractory_par2 = 31.64
                               )


trajectories <- res$simulation_results

pars <- res$pars

final_sizes <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),
                          final_size=res$final_size)
paralysis_totals <- res$n_paralysis
data_consistent <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),
                              consistent=res$data_are_consistent)
paralysis_totals <- paralysis_totals %>% left_join(data_consistent) %>% filter(consistent==TRUE)


## Get outbreak length
run_times <- trajectories %>%
    group_by(sim) %>% filter(t == max(t)) %>%
    rename(tmax = t) %>%
    rename(inc_end = inc) %>%
    mutate(ongoing = inc_end>0)

t_starts <- trajectories %>% filter(para == 1) %>% group_by(sim) %>% filter(t == min(t)) %>% rename(tstart=t) %>% select(sim, tstart)


## Any infections in the X days prior to the end
x <- 10
run_times_long <- trajectories %>%
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
    mutate(Re=R0*prop_immune_groups.1) %>%
    left_join(run_times) %>% 
    left_join(run_times_long) %>%
    left_join(final_sizes) %>%
    left_join(t_starts) %>%
    mutate(date_start = as.Date("2022-06-22")-tstart)

## Make t a date
trajectories <- trajectories %>% left_join(pars %>% select(sim, date_start)) %>%
    mutate(t = as.Date(date_start + t)) %>%
    drop_na() 


## Flag which runs are consistent with the data
use_samps <- paralysis_totals %>% filter(consistent == TRUE) %>% pull(sim)
print(use_samps)

trajectories <- trajectories %>% filter(sim %in% use_samps)

## Rockland county trajectories only valid if they have some cases in May, June, July and August
trajectories <- trajectories %>% mutate(month = round_date(t, "month"))
trajectories %>% filter(month %in% as.Date(c("2022-05-01","2022-06-01","2022-07-01","2022-08-01"))) %>% group_by(sim, month) %>% dplyr::summarize(monthly_inc=sum(inc)) %>% mutate(keep=monthly_inc > 0) %>% select(sim, keep, month) %>% pivot_wider(id_cols = sim,values_from="keep",names_from="month") %>%filter(`2022-05-01`==TRUE,`2022-06-01`==TRUE,`2022-07-01`==TRUE,`2022-08-01`==TRUE) %>% pull(sim) -> use_samps

trajectories <- trajectories %>% filter(sim %in% use_samps)
pars <- pars %>%
    filter(sim %in% use_samps)

save(pars, file=paste0("sims/simulation_",i,".RData"))

## Rerun for NYC-like place
nyc <- NULL
## Need to improve this, but can seed at 1 infection per day for 2 weeks from start. Seeding would likely bey
for(j in 1:nrow(pars)){
    tmp <- run_simulation_twoimmune(R0=pars$R0[j], 
                                    rel_R0=0.089,
                                    P=8336817,
                                    ini_infs=rep(1,14),
                             observed_data=NULL,
                             continue_run=TRUE,
                             tmax=500,
                             infect_rate=0.313,
                             infect_shape=3.96,
                             infect_partial_shape=2.48,
                             infect_partial_rate=0.243,
                             incu_scale=pars$incu_scale[j],
                             incu_shape=pars$incu_shape[j],
                             prob_paralysis_s = pars$prob_paralysis_s[j], 
                             prob_paralysis_ps = pars$prob_paralysis_ps[j],
                             prop_immune_groups = c(0.158,(1-0.158-0.169),0.169))#as.numeric(pars[j,c("prop_immune_groups.1","prop_immune_groups.2","prop_immune_groups.3")]))
    
    dat <- tmp$dat
    dat$sim <- pars$sim[j]
    nyc[[j]] <- dat
}
nyc <- do.call("bind_rows",nyc)
save(nyc, file=paste0("sims_nyc/nyc_",i,".RData"))


## Restart simulations where they left off with/without intervention strategy
if(length(use_samps) > 1){
    vacc_strats <- matrix(c(0,0,
                            0.2,0.2,
                            0.4,0.4,
                            0.6,0.6,
                            0.8,0.8),ncol=2,byrow=TRUE)
    res2 <- restart_simulations_table_twoimmune(use_sims = use_samps,
                              pars=res$pars,
                              final_conditions = res$final_conditions,
                              t_starts=res$tmax_vector,
                              vaccinate_proportion = vacc_strats,
                              P=325000,
                              tmax=500,nruns=1)
    save(res2, file=paste0("sims_traj/traj_",i,".RData"))
    
}
