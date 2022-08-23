
run_simulation <- function(R0=2, 
                           observed_data=c(1, rep(0,30)),
                           infect_scale=0.71, infect_shape=9.8, incu_scale=0.625,incu_shape=25.6,
                           tmax=100, ini_infs=5, n_indiv=10000,max_infectious_period=50,
                           prob_paralysis=1/2000, P=10000000, S_ini=1,
                           inf_max=P){
    incubation_period <- function(n){extraDistr::rdgamma(n,scale=incu_scale,shape=incu_shape)}
    infectiousness <- function(t){extraDistr::ddgamma(t, scale=infect_scale, shape=infect_shape)}
    first_paralysis <- tmax
    susceptible <- rep(0, tmax)
    susceptible[1] <- floor(P*S_ini)
    
    new_infections <- rep(0, tmax)
    new_infections[1] <- ini_infs
    susceptible[2] <- susceptible[1] - ini_infs
    paralysis_incidence <- rep(0, tmax)
    
    cur_inf_number <- ini_infs + 1
    changed_first <- FALSE
    t <- 2
    
    total_paralysis_cases <- 0
    
    ## Check if simulated paralysis incidence is consistent with observed data
    data_are_consistent <- FALSE
    obs_dur <- length(observed_data)
    
    ## Solve up to the first paralysis case, and then keep solving until
    ## the paralysis incidence curve is no longer consistent with the data
    case_detect <- FALSE
    while((t <= first_paralysis | (t < (first_paralysis + obs_dur) & data_are_consistent)) & t <= tmax & sum(new_infections) < inf_max){

        ## Simulate new infections caused today by currently infected individuals
        ## Poisson draw of new infections with lambda = R0*g(t) where g(t) is infectiousness 
        ## on day t of infection
        min_t <- max(1, t-max_infectious_period)
        use_infections <- new_infections[min_t:t]
        cur_infectious <- sum(use_infections)
        
        ## If at least one infectious person
        if(cur_infectious >= 1){
            tmp_infectiousness <- R0*infectiousness((t - min_t:t)-1)
            inc <- sum(sapply(seq_along(min_t:t), 
                              function(x) sum(rpois(1, use_infections[x]*tmp_infectiousness[x]))))
            inc <- floor(susceptible[t-1]*(1-exp(-(inc/P))))
            susceptible[t] <- susceptible[t-1] - inc
            
            ## If we simulated an infection
            if(inc >= 1){
                ## Set infection states of new infections and record time of infection
                new_infections[t] <- inc
                
                ## Simulate paralysis cases from these new infections
                paralysis_cases <- rbinom(1, inc, prob_paralysis)
                total_paralysis_cases <- total_paralysis_cases + paralysis_cases
                
                ## If we generated a paralysis case, simulation incubation periods and
                ## update the paralysis incidence curve
                if(paralysis_cases > 0){
                    ## Randomly sample which infections will have paralysis and then simulate incubation periods
                  
                    incubation_periods <- incubation_period(paralysis_cases)
                    incu_table <- table(incubation_periods)
                    indices <- t + as.numeric(incu_table %>% names)
                    paralysis_incidence[indices] <- paralysis_incidence[indices] + as.numeric(incu_table)
                    first_paralysis <- min(first_paralysis, min(t + incubation_periods))
                }
                
            }
            
        }
        ## If we've had a paralysis case, start checking for consistency
        if(t >= first_paralysis & total_paralysis_cases > 0){        
            ## If this is the first time we see a paralysis case, then the data are temporarily consistent
            t_elapsed <- t - first_paralysis
            data_are_consistent <- isTRUE(all.equal(paralysis_incidence[first_paralysis:t],
                                                    observed_data[1:(1+t_elapsed)]))
        }
        
        ## Track how many infections we've generated
        cur_inf_number <- cur_inf_number + new_infections
        t <- t + 1
    }
    if(sum(new_infections) > inf_max) data_are_consistent <- FALSE
    if(data_are_consistent == TRUE & sum(paralysis_incidence[1:(t-1)]) > 1) browser()
    return(list(data.frame(t=1:(t-1), inc=new_infections[1:(t-1)]),
                n_paralysis=sum(paralysis_incidence[1:(t-1)]),
                data_are_consistent=data_are_consistent,
                paralysis_incidence=data.frame(t=1:(t-1), para=paralysis_incidence[1:(t-1)])))
}
rdgamma_by_mean <- function(n, mean1, var1){
    scale <- var1/mean1
    shape <- mean1/scale
    extraDistr::rdgamma(n, shape=shape, scale=scale)
}
rgamma_by_mean <- function(n, mean1, var1){
    scale <- var1/mean1
    shape <- mean1/scale
    rgamma(n, shape=shape, scale=scale)
}
get_beta_pars <- function(mu, var){
    a = mu*((mu*(1-mu) / var)  -1)
    b = a*(1-mu)/mu
    return(c(a,b))
}

random_simulation <- function( n=100,  
                               tmax=100,
                               max_infectious_period=50,
                               incu_mean_prior_mean=16,incu_mean_prior_var=5,
                               incu_var_prior_mean=10,incu_var_prior_var=3,
                               infect_mean_prior_mean=7,infect_mean_prior_var=5,
                               infect_var_prior_mean=5,infect_var_prior_var=3,
                               R0_prior_mean=1,R0_prior_var=0.5,
                               prob_paralysis_mean=1/2000,prob_paralysis_var=1/2000 / 1000,
                               prop_immune_mean = 0.5,prop_immune_var=0.001,
                               R0_min=0,R0_max=10,
                               ini_infs=5,
                               index_start=1
){
    ## Incubation period
    incu_mean_scale <- incu_mean_prior_var/incu_mean_prior_mean
    incu_mean_shape <- incu_mean_prior_mean/incu_mean_scale
    
    incu_var_scale <- incu_var_prior_var/incu_var_prior_mean
    incu_var_shape <- incu_var_prior_mean/incu_var_scale
    
    incu_mean <- rgamma(n, scale=incu_mean_scale, shape=incu_mean_shape)
    incu_var <- rgamma(n, scale=incu_var_scale, shape=incu_var_shape)
    
    incu_scale <- incu_var/incu_mean
    incu_shape <- incu_mean/incu_scale
    
    ## Infectious period
    infect_mean_scale <- infect_mean_prior_var/infect_mean_prior_mean
    infect_mean_shape <- infect_mean_prior_mean/infect_mean_scale
    
    infect_var_scale <- infect_var_prior_var/infect_var_prior_mean
    infect_var_shape <- infect_var_prior_mean/infect_var_scale
    
    infect_mean <- rgamma(n, scale=infect_mean_scale, shape=infect_mean_shape)
    infect_var <- rgamma(n, scale=infect_var_scale, shape=infect_var_shape)
    
    infect_scale <- infect_var/infect_mean
    infect_shape <- infect_mean/infect_scale
    
    R0_scale <- R0_prior_var/R0_prior_mean
    R0_shape <- R0_prior_mean/R0_scale
    
    R0 <- rgamma(n, scale=R0_scale,shape=R0_shape)
    R0 <- runif(n, R0_min, R0_max)
    
    para_pars <- get_beta_pars(prob_paralysis_mean, prob_paralysis_var)
    ##prob_para <- rbeta(n, para_pars[1],para_pars[2])
    
    prob_para <- rnorm(n, prob_paralysis_mean, sqrt(prob_paralysis_var))
    prob_para[prob_para < 0] <- 0
    prob_para[prob_para >1] <- 1
    
    immune_prop_pars <- get_beta_pars(prop_immune_mean, prop_immune_var)
    prop_immune <- rbeta(n, immune_prop_pars[1],immune_prop_pars[2])
    
    
    res <- NULL
    res_par <- NULL
    final_size <- numeric(n)
    incubation_periods <- NULL
    infectious_periods <- NULL
    n_paralysis <- numeric(n)
    data_consistent <- logical(n)
    
    for(i in 1:n){
        if(i %% 100 == 0) print(paste0("Sim number: ", i, " of ", n))
        tmp <- run_simulation(R0=R0[i], 
                              tmax=tmax,
                              infect_scale=infect_scale[i],
                              prob_paralysis=prob_para[i],
                              ini_infs=ini_infs,
                              S_ini = prop_immune[i],
                              infect_shape=infect_shape[i],
                              incu_scale=incu_scale[i],incu_shape=incu_shape[i])
        
        sim_no <- i + index_start - 1
        tmp[[1]]$sim <- sim_no
        tmp[[length(tmp)]]$sim <- sim_no
        res[[i]] <- tmp[[1]]
        res_par[[i]] <- tmp[[length(tmp)]]
        final_size[i] <- sum(tmp[[1]]$inc)
        
        n_paralysis[i] <- sum(tmp$n_paralysis)
        infectious_periods[[i]] <- data.frame(sim=sim_no, t=1:max_infectious_period - 1, 
                                              prob=extraDistr::ddgamma(1:max_infectious_period, 
                                                                       shape=infect_shape[i],
                                                                       scale=infect_scale[i]))
        
        
        incubation_periods[[i]] <- data.frame(sim=sim_no, t=1:50 - 1, 
                                              prob=extraDistr::ddgamma(1:50, 
                                                                       shape=incu_shape[i],
                                                                       scale=incu_scale[i]))
        data_consistent[i] <- tmp$data_are_consistent
        
    }
    par_table <- data.frame(sim=index_start:(index_start + n - 1),
                            R0=R0,infectious_period_mean=infect_mean,
                            prop_paralysis=prob_para,
                            prop_immune=prop_immune,
                            infectious_period_var=infect_var,
                            incubation_period_mean=incu_mean, incubation_period_var=incu_var)
    res <- do.call("bind_rows",res)
    res_par <- do.call("bind_rows",res_par)
    infectious_periods <- do.call("bind_rows",infectious_periods)
    incubation_periods <- do.call("bind_rows",incubation_periods)
    
    list(res=res, res_par=res_par,
         par_table=par_table,final_size=final_size,
         infectious_periods=infectious_periods,incubation_periods=incubation_periods,
         n_paralysis=data.frame(sim=index_start:(index_start+n-1), paralysis=n_paralysis),
         data_are_consistent=data_consistent)
}