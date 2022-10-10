library(dplyr)
library(tidyr)
library(tidyverse)
library(lubridate)

#setwd("~/paralytic_polio_estimates")



setwd("C:/Users/Mary/OneDrive - Harvard University/Github copy/paralytic_polio_estimates")

priors <- read.csv("pars/priors - alt low coverage scenario.csv")

source("simulation_functions_twoimmune.R")

pop_size <- 91000

nsims <- 5e4

{
    i <- 3
    print(i)
    set.seed(i)
    
    save_wd <- "sims alt low coverage"
    
    print("Scenario: alternative low coverage (rockland county)")
    ## Filter the table to the correct scenario
    tmp_pars <- priors
    
    
    
    ## Main function call. There are a lot of pipes here, but they are simply
    ## extracting the correct parameter value for each argument from 
    ## the prior table
    res <- random_simulation_twoimmune(n=nsims, ## How many draws?
                                       ## Maximum length of the simulation,
                                       ## population size and seed size
                                       tmax=180, P=pop_size,
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
                                       prop_refractory_par2 = tmp_pars %>% filter(`model.parameter` == "prop_refractory") %>% pull(par2))
    
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
    if(!dir.exists(save_wd)) dir.create(save_wd)
    save(pars, file=paste0(save_wd,"/simulation_",i,".RData"))
    
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
                                                    P=pop_size,tmax=500, ## Pop size and simulation duration.
                                                    ## If desired, can run multiple restarts
                                                    ## for each trajectory.
                                                    nruns=1)
        
        ## Save the outputs
        if(!dir.exists(paste0(save_wd,"_traj/"))) dir.create(paste0(save_wd,"_traj/"))
        
        save(res2, file=paste0(save_wd,"_traj/traj_",i,".RData"))
    }
}
