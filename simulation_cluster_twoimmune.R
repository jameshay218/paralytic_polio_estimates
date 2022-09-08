library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)

setwd("~/paralytic_polio_estimates")
#setwd("~/Documents/GitHub/paralytic_polio_estimates")

priors <- read_csv("pars/priors for model.csv")

source("simulation_functions_twoimmune.R")

nsims <- 20000

i <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
#i <- 1
print(i)
set.seed(i)

## We have two scenarios: the high immunity and low immunity scenarios. We call the same function for each scenario, but filter to use the entries in `priors` corresponding to the correct simulation.
scenarios <- c("rockland_high_coverage","rockland_low_coverage")
save_wds <- c("sims","sims_low_coverage")
for(index in 1:2){
    print(paste0("Scenario: ", scenarios[index]))
    ## Filter the table to the correct scenario
    tmp_pars <- priors %>% filter(scenario == scenarios[index])
    
    ## Main function call. There are a lot of pipes here, but they are simply
    ## extracting the correct parameter value for each argument from 
    ## the prior table
    res <- random_simulation_twoimmune(n=nsims, ## How many draws?
                                       ## Maximum length of the simulation,
                                       ## population size and seed size
                                       tmax=180, P=340000,
                                       ini_infs=1,
                                       
                                       ## Vector of observed data
                                       observed_data=c(1, rep(0,49)),
                                       ## Record simulation index. This just makes sure that the starting index is higher if i, set on L17, is not 1
                                       index_start=(i-1)*nsims + 1,
                                       
                                       
                                   incu_mean_prior_mean=14,
                                   incu_mean_prior_var=3,
                                   incu_var_prior_mean=15,
                                   incu_var_prior_var=1,
                                   gen_interval_susc_shape_par1=tmp_pars %>% 
                                       filter(`model.parameter` == "susceptible_generation_interval_shape") %>% pull(par1),
                                   gen_interval_susc_shape_par2=tmp_pars %>% 
                                       filter(`model.parameter` == "susceptible_generation_interval_shape") %>% pull(par2),
                                   gen_interval_susc_rate_par1=tmp_pars %>% 
                                       filter(`model.parameter` == "susceptible_generation_interval_rate") %>% pull(par1),
                                   gen_interval_susc_rate_par2=tmp_pars %>% 
                                       filter(`model.parameter` == "susceptible_generation_interval_rate") %>% pull(par2),
                                   gen_interval_partial_shape_par1=tmp_pars %>% 
                                       filter(`model.parameter` == "partial_generation_interval_shape") %>% pull(par1),
                                   gen_interval_partial_shape_par2=tmp_pars %>% 
                                       filter(`model.parameter` == "partial_generation_interval_shape") %>% pull(par2),
                                   gen_interval_partial_rate_par1=tmp_pars %>% 
                                       filter(`model.parameter` == "partial_generation_interval_rate") %>% pull(par1),
                                   gen_interval_partial_rate_par2=tmp_pars %>% 
                                       filter(`model.parameter` == "partial_generation_interval_rate") %>% pull(par2),
                                   
                                   ## Paralysis probability priors
                                   prob_paralysis_mean=0.0005,
                                   prob_paralysis_ps_par1 = -4,
                                   prob_paralysis_ps_par2 = -1,
                                   prob_paralysis_var = 2e-8,
                                   #prob_paralysis_ps_var = 1e-11,
                                   R0_dist="truncnorm",
                                   R0_par1=4.9,R0_par2=2,
                                   rel_R0_par1=tmp_pars %>% filter(`model.parameter` == "relative_infectiousness") %>% pull(par1),
                                   rel_R0_par2=tmp_pars %>% filter(`model.parameter` == "relative_infectiousness") %>% pull(par2),
                                   prop_susceptible_par1 = tmp_pars %>% filter(`model.parameter` == "prop_susceptible") %>% pull(par1),
                                   prop_susceptible_par2 = tmp_pars %>% filter(`model.parameter` == "prop_susceptible") %>% pull(par2),
                                   prop_refractory_par1 = tmp_pars %>% filter(`model.parameter` == "prop_refractory") %>% pull(par1),
                                   prop_refractory_par2 = tmp_pars %>% filter(`model.parameter` == "prop_refractory") %>% pull(par2)
                                   )
    
    ## Pull out all trajectories, number of paralysis cases etc
    trajectories <- res$simulation_results
    pars <- res$pars
    final_sizes <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),
                              final_size=res$final_size)
    paralysis_totals <- res$n_paralysis
    
    ## Create a data frame which flags each simulation ID as consistent
    ## with the data or not, then filter to only those which are consistent
    data_consistent <- data.frame(sim=((i-1)*nsims + 1):((i-1)*nsims + nsims),
                                  consistent=res$data_are_consistent)
    paralysis_totals <- paralysis_totals %>% left_join(data_consistent) %>% filter(consistent==TRUE)
    
    
    ## Find how long each simulation ran for and flag if it still had
    ## ongoing incidence by the end.
    run_times <- trajectories %>%
        group_by(sim) %>% filter(t == max(t)) %>%
        rename(tmax = t) %>%
        rename(inc_end = inc) %>%
        mutate(ongoing = inc_end>0)
    
    ## For each trajectory, find the first day of paralysis onset and call this `t_start`. This will be used to set the dates later on.
    t_starts <- trajectories %>% filter(para == 1) %>% group_by(sim) %>% filter(t == min(t)) %>% rename(tstart=t) %>% select(sim, tstart)
    
    
    ## Check if there were ANY new infections in the x days prior
    ## to the end of the simulation. if so, we flag this trajectory
    ## as "ongoing"
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
    
    
    ## Merge in all of the time-related parameters calculated above with
    ## the model.parameters
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
    use_samps <- paralysis_totals %>% filter(consistent == TRUE) %>% 
        pull(sim)
    print(use_samps)
    
    ## Only keep trajectories which are consistent with the data
    trajectories <- trajectories %>% filter(sim %in% use_samps)
    
    ## WASTEWATER CONSISTENCY CHECK
    ## Rockland county trajectories only valid if they have some cases in May, June, July and August
    ## Get date by month
    trajectories <- trajectories %>% mutate(month = round_date(t, "month"))
    
    ## For each month in May through August, find the total incidence in 
    ## each month, flag if it was >0, and then only keep trajectories with
    ## incidence in each month.
    use_samps <- trajectories %>% 
        filter(month %in% as.Date(c("2022-05-01","2022-06-01","2022-07-01","2022-08-01"))) %>% 
        group_by(sim, month) %>% dplyr::summarize(monthly_inc=sum(inc)) %>% 
        mutate(keep=monthly_inc > 0) %>% 
        select(sim, keep, month) %>% 
        pivot_wider(id_cols = sim,values_from="keep",names_from="month") %>%
        filter(`2022-05-01`==TRUE & `2022-06-01`==TRUE & `2022-07-01`==TRUE & `2022-08-01`==TRUE) %>% 
        pull(sim)
    
    trajectories <- trajectories %>% filter(sim %in% use_samps)
    pars <- pars %>%
        filter(sim %in% use_samps)
    
    ## Save the parameter draws
    if(!dir.exists(save_wds[index])) dir.create(save_wds[index])
    save(pars, file=paste0(save_wds[index],"/simulation_",i,".RData"))

    ## Restart simulations where they left off with/without intervention strategy
    if(length(use_samps) > 1){
        ## The trajectories will be extended for each row of `vacc_strats`,
        ## but for the initial version we only care about the first row (where no vaccination is performed).
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
                                  P=340000,tmax=500, ## Pop size and simulation duration.
                                  ## If desired, can run multiple restarts
                                  ## for each trajectory.
                                  nruns=1)
        
        ## Save the outputs
        if(!dir.exists(paste0(save_wds[index],"_traj/"))) dir.create(paste0(save_wds[index],"_traj/"))
        
        save(res2, file=paste0(save_wds[index],"_traj/traj_",i,".RData"))
    }
}

###################################
## NYC SIMULATIONS
###################################
## For each prior draw which was consistent with the Rockland County data, 
## start a new simulation instead using parameters for population immunity/size for NYC. The for loop simply goes through each row of the saved parameters and creates one new trajectory per row.
nyc <- NULL
load(paste0(save_wds[index],"/simulation_",i,".RData"))
tmp_pars <- priors %>% filter(scenario == "NYC")
for(j in 1:nrow(pars)){
    print(j)
    ## Sample generation interval parameters for susceptible
    infect_rate_par1 <- tmp_pars %>% filter(`model.parameter` == "susceptible_generation_interval_rate") %>% pull(par1)
    infect_rate_par2 <- tmp_pars %>% filter(`model.parameter` == "susceptible_generation_interval_rate") %>% pull(par2)
    
    infect_rate_nyc <- rbeta(1, infect_rate_par1, infect_rate_par2)
    
    
    infect_shape_par1 <- tmp_pars %>% filter(`model.parameter` == "susceptible_generation_interval_shape") %>% pull(par1)
    infect_shape_par2 <- tmp_pars %>% filter(`model.parameter` == "susceptible_generation_interval_shape") %>% pull(par2)
    
    infect_shape_nyc <- truncnorm::rtruncnorm(1, a=0, mean=infect_shape_par1, sd=infect_shape_par2)
    
    ## Sample generation interval parameters for partially immune
    infect_rate_ps_par1 <- tmp_pars %>% filter(`model.parameter` == "partial_generation_interval_rate") %>% pull(par1)
    infect_rate_ps_par2 <- tmp_pars %>% filter(`model.parameter` == "partial_generation_interval_rate") %>% pull(par2)
    
    infect_rate_ps_nyc <- rbeta(1, infect_rate_ps_par1, infect_rate_ps_par2)
    
    infect_shape_ps_par1 <- tmp_pars %>% filter(`model.parameter` == "partial_generation_interval_shape") %>% pull(par1)
    infect_shape_ps_par2 <- tmp_pars %>% filter(`model.parameter` == "partial_generation_interval_shape") %>% pull(par2)
    
    infect_shape_ps_nyc <- truncnorm::rtruncnorm(1, a=0, mean=infect_shape_ps_par1, sd=infect_shape_ps_par2)
    
    ## Sample relative infectiousness
    rel_infect_par1 <- tmp_pars %>% filter(`model.parameter` == "relative_infectiousness") %>% pull(par1)
    rel_infect_par2 <- tmp_pars %>% filter(`model.parameter` == "relative_infectiousness") %>% pull(par2)
    
    rel_infect_nyc <- rbeta(1, rel_infect_par1, rel_infect_par2)
    
    ## Sample proportion immune
    prop_susceptible_par1 <- tmp_pars %>% filter(`model.parameter` == "prop_susceptible") %>% pull(par1)
    prop_susceptible_par2 <- tmp_pars %>% filter(`model.parameter` == "prop_susceptible") %>% pull(par2)
    
    prop_susceptible_nyc <- rbeta(1, prop_susceptible_par1, prop_susceptible_par2)
    
    prop_refractory_par1 <- tmp_pars %>% filter(`model.parameter` == "prop_refractory") %>% pull(par1)
    prop_refractory_par2 <- tmp_pars %>% filter(`model.parameter` == "prop_refractory") %>% pull(par2)
    
    prop_refractory_nyc <- rbeta(1, prop_refractory_par1, prop_refractory_par2)
    
    prop_immune_nyc <- c(prop_refractory_nyc, 1-prop_refractory_nyc-prop_susceptible_nyc,prop_susceptible_nyc)
    prop_immune_nyc <- prop_immune_nyc/sum(prop_immune_nyc)
    
    tmp <- run_simulation_twoimmune(R0=pars$R0[j], 
                                    rel_R0=rel_infect_nyc,
                                    P=8500000, ## Population size of NYC
                                    ini_infs=rep(1,14), ## Vector of daily seeding from t0
                                    observed_data=NULL, ## Flags to ensure the simulation just runs for 500 time steps
                                    continue_run=TRUE,
                                    tmax=500,
                                    infect_rate=infect_rate_nyc,
                                    infect_shape=infect_shape_nyc,
                                    infect_partial_shape=infect_shape_ps_nyc,
                                    infect_partial_rate=infect_rate_ps_nyc,
                                    incu_scale=pars$incu_scale[j],
                                    incu_shape=pars$incu_shape[j],
                                    prob_paralysis_s = pars$prob_paralysis_s[j], 
                                    prob_paralysis_ps = pars$prob_paralysis_ps[j],
                                    prop_immune_groups = prop_immune_nyc)
    
    dat <- tmp$dat
    dat$sim <- pars$sim[j]
    nyc[[j]] <- dat
}
nyc <- do.call("bind_rows",nyc)
if(!dir.exists(paste0("sims_nyc/"))) dir.create(paste0("sims_nyc/"))
save(nyc, file=paste0("sims_nyc/nyc_",i,".RData"))

