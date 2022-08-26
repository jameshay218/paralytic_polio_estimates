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
is.scalar <- function(x) is.atomic(x) && length(x) == 1L

find_gamma_pars <- function(gamma_mean, gamma_var){
    gamma_scale <- gamma_var/gamma_mean
    gamma_shape <- gamma_mean/gamma_scale
    if(is.scalar(gamma_mean)){
        return(c(gamma_scale,gamma_shape))
    } else {
        return(list(gamma_scale,gamma_shape))
    }
}

## Run simulation with two types of immunity -- fully immune and partially immune
run_simulation_twoimmune <- function(
                            ## R0 and relative reduction in R0 for partially immune popn
                            R0=2, rel_R0=1, 
                            ## Vector of observed paralysis data
                            observed_data=c(1, rep(0,48)),
                            ## Generation interval distribution of fully susceptible population
                            latent_period=3,
                            infect_scale=0.71, infect_shape=9.8, 
                            ## Generation interval distribution of partially immune population
                            #infect_scale_ps=0.71, infect_shape_ps=9.8, 
                            infect_ps_par=0.231,
                            ## Paralysis onset delay distribution
                            incu_scale=0.625,incu_shape=25.6,
                            ## Run time, seed size, population size and maximum infectiousness period
                            tmax=100, ini_infs=5, max_infectious_period=50,max_incu_period=50,
                            P=10000000, 
                            ## Probability of paralysis for fully susceptible and partially immune
                            prob_paralysis_s=1/2000,prob_paralysis_ps=1/5000,
                            ## Initial distribution of susceptible and partially immune individuals
                            prop_immune_groups,
                            ## Continue the simulation after fitting the data
                            continue_run=FALSE, restart_simulation=FALSE, final_conditions=NULL, t_start=NULL, 
                            ## Vector of length 2: proportion of fully susceptible -> partially immune, 
                            ## and partially_immune -> fully immune
                            vaccinate_proportion=NULL
                            ){
    incubation_period <- function(n){pmax(extraDistr::rdgamma(n,scale=incu_scale,shape=incu_shape),max_incu_period)}
    infectiousness <- function(t){extraDistr::ddgamma(t-latent_period, scale=infect_scale, shape=infect_shape)/extraDistr::pdgamma(max_infectious_period, shape=infect_shape,scale=infect_scale)}
    #infectiousness_ps <- function(t){c(rep(0, latent_period), extraDistr::ddgamma(t, scale=infect_scale_ps, shape=infect_shape_ps))}
    infectiousness_ps <- function(t){dexp(t-latent_period, infect_ps_par)/pexp(max_infectious_period, infect_ps_par)}

    ## Initial number who are fully immune, partially immune and fully susceptible
    ini_pop <- rmultinom(1, P, prop_immune_groups)[,1]
if(restart_simulation){
        tmax1 <- tmax
        
        new_infections_s <- c(final_conditions$incidence_s, rep(0, (tmax - t_start)))
        new_infections_ps <- c(final_conditions$incidence_ps, rep(0, (tmax - t_start)))
        
        fully_susceptible <- c(final_conditions$fully_susceptible, rep(0, (tmax - t_start)))
        partially_susceptible <- c(final_conditions$partially_susceptible, rep(0, (tmax - t_start)))
        
        paralysis_incidence_s <- c(final_conditions$paralysis_s, rep(0, (tmax - t_start)))
        paralysis_incidence_ps <- c(final_conditions$paralysis_ps, rep(0, (tmax - t_start)))

        total_paralysis_cases <- sum(paralysis_incidence_s) + sum(paralysis_incidence_ps)
        
        t <- t_start 
        if(!is.null(vaccinate_proportion)){
            ## Vaccinate some fraction of fully susceptible
            fully_s_to_partially_immune <- rbinom(1, fully_susceptible[t-1],vaccinate_proportion[1])
            ## Vaccinate some fraction of partially immune people
            partially_immune_to_fully_immune <- rbinom(1, partially_susceptible[t-1], vaccinate_proportion[2])

            fully_susceptible[t-1] <- fully_susceptible[t-1] - fully_s_to_partially_immune
            partially_susceptible[t-1] <- partially_susceptible[t-1] + fully_s_to_partially_immune
            partially_susceptible[t-1] <- partially_susceptible[t-1] - partially_immune_to_fully_immune
            Rt <- final_conditions$Rt
        }
    } else {
        tmax1 <- tmax + max_incu_period
        ## Tracking all susceptible individuals, and partially/fully susceptible separately
        fully_susceptible <- rep(0, tmax + max_incu_period)
        partially_susceptible <- rep(0, tmax + max_incu_period)
        
        fully_susceptible[1] <- ini_pop[3]
        partially_susceptible[1] <- ini_pop[2]

        ## Initial infections in the fully susceptible group
        new_infections_s <- rep(0, tmax + max_incu_period)
        new_infections_s[1] <- ini_infs
        new_infections_ps <- rep(0, tmax + max_incu_period)
        
        fully_susceptible[2] <- fully_susceptible[1] - ini_infs
        
        paralysis_incidence_s <- rep(0, tmax + max_incu_period)
        paralysis_incidence_ps <- rep(0, tmax + max_incu_period)
        total_paralysis_cases <- 0
        t <- 2
        
        Rt <- rep(0, tmax1)
        
        Rt[1] <- R0 * prop_immune_groups[3]/sum(prop_immune_groups[2:3]) + R0*rel_R0*  prop_immune_groups[2]/sum(prop_immune_groups[2:3])
        Rt[1] <- Rt[1] * (fully_susceptible[1]+partially_susceptible[1])/P
    }
    
    ## Check if simulated paralysis incidence is consistent with observed data
    data_are_consistent <- FALSE
    obs_dur <- length(observed_data)
    first_paralysis <- tmax
    ## Solve up to the first paralysis case, and then keep solving until
    ## the paralysis incidence curve is no longer consistent with the data
    ## Or, keep solving after the observed data if the continue_run flag is TRUE
    while((t <= first_paralysis | (t < (first_paralysis + obs_dur) & data_are_consistent) | 
          (continue_run & t <=tmax)) & t <= tmax){
        ## Simulate new infections caused today by currently infected individuals
        ## Poisson draw of new infections with lambda = R0*g(t) where g(t) is infectiousness 
        ## on day t of infection
        min_t <- max(1, t-max_infectious_period)
        use_infections_s <- new_infections_s[min_t:t]
        use_infections_ps <- new_infections_ps[min_t:t]
        cur_infectious <- sum(use_infections_s) + sum(use_infections_ps)
        
        ## If at least one infectious person
        if(cur_infectious >= 1){
            ## Different infectiousness profile for susceptible and partially immune groups
            tmp_infectiousness_s <- R0*infectiousness((t - min_t:t))
            tmp_infectiousness_ps <- rel_R0*R0*infectiousness_ps((t-min_t:t))
            
            
            ## Number of new infections arising from previously susceptible and partially immune groups
            inc_s <- sum(sapply(seq_along(min_t:t), function(x) sum(rpois(1, use_infections_s[x]*tmp_infectiousness_s[x]))))
            inc_ps <- sum(sapply(seq_along(min_t:t), function(x) sum(rpois(1, use_infections_ps[x]*tmp_infectiousness_ps[x]))))
            
            
            ## Get number of contacts with susceptible/immune individuals
            inc <- inc_s + inc_ps
            inc_s <- rbinom(1,fully_susceptible[t-1],prob=(1-exp(-(inc/P))))
            inc_ps <- rbinom(1,partially_susceptible[t-1],prob=(1-exp(-(inc/P))))
            #if(inc > 25) browser()
            Rt[t] <- R0 * prop_immune_groups[3]/sum(prop_immune_groups[2:3]) + R0*rel_R0*  prop_immune_groups[2]/sum(prop_immune_groups[2:3])
            
            Rt[t] <- Rt[t] * (fully_susceptible[t-1]+partially_susceptible[t-1])/P
            
            ## Update susceptible pool
            fully_susceptible[t] <- fully_susceptible[t-1] - inc_s
            partially_susceptible[t] <- partially_susceptible[t-1] - inc_ps
            
            ## If we simulated an infection
            if(inc >= 1){
                ## Set infection states of new infections and record time of infection
                new_infections_s[t] <- inc_s
                new_infections_ps[t] <- inc_ps
                
                ## Simulate paralysis cases from these new infections
                paralysis_cases_s <- rbinom(1, inc_s, prob_paralysis_s)
                paralysis_cases_ps <- rbinom(1, inc_ps, (1.0-prob_paralysis_ps)*prob_paralysis_s)
                #paralysis_cases_ps <- rbinom(1, inc_ps, prob_paralysis_ps)
                #print(prob_paralysis_ps)
                total_paralysis_cases <- total_paralysis_cases + paralysis_cases_s + paralysis_cases_ps
                
                new_paralysis <- paralysis_cases_s + paralysis_cases_ps
                
                ## If we generated a paralysis case, simulate incubation periods and
                ## update the paralysis incidence curve
                if(new_paralysis > 0){
                    ## Randomly sample which infections will have paralysis and then simulate incubation periods
                    incu_periods <- c()
                    ## Fully susceptible individuals
                    if(paralysis_cases_s){
                        incubation_periods_s <- incubation_period(paralysis_cases_s)
                        incu_table_s <- table(incubation_periods_s)
                        indices <- t + as.numeric(incu_table_s %>% names)
                        paralysis_incidence_s[indices] <- paralysis_incidence_s[indices] + as.numeric(incu_table_s)
                        incu_periods <- c(incu_periods, incubation_periods_s)
                    }
                    
                    ## Partially immune individuals
                    if(paralysis_cases_ps){
                        incubation_periods_ps <- incubation_period(paralysis_cases_ps)
                        incu_table_ps <- table(incubation_periods_ps)
                        indices <- t + as.numeric(incu_table_ps %>% names)
                        paralysis_incidence_ps[indices] <- paralysis_incidence_ps[indices] + as.numeric(incu_table_ps)
                        incu_periods <- c(incu_periods, incubation_periods_ps)
                        
                    }
                    ## Keep track of when the first observed case of paralysis is
                    first_paralysis <- min(first_paralysis, min(t + incu_periods))
                }
            }
            
        }
        ## If we've had a paralysis case, start checking for consistency
        if(!restart_simulation & t >= first_paralysis & total_paralysis_cases > 0 & t < (first_paralysis + obs_dur)){    
            ## If this is the first time we see a paralysis case, then the data are temporarily consistent
            #print(paste0("First paralysis: ", first_paralysis))
            t_elapsed <- t - first_paralysis
            #print(paste0("Time elapsed: ", t_elapsed))
            #print("Sim cases: ")
            #print(paralysis_incidence_s[first_paralysis:t]+paralysis_incidence_ps[first_paralysis:t])
            #print("Obs cases: ")
            #print(observed_data[1:(1+t_elapsed)])
            data_are_consistent <- isTRUE(all.equal(paralysis_incidence_s[first_paralysis:t]+paralysis_incidence_ps[first_paralysis:t],
                                                    observed_data[1:(1+t_elapsed)]))

            #if(!data_are_consistent) print("BROKEN")
        }
        t <- t + 1
    }
    if(t >= tmax) data_are_consistent <- FALSE
    #if(data_are_consistent){
    #    print("SUCCESS")
    #    browser()
    #    break
    #}
    ## Inconsistent if generated more infections then there are people
    if((sum(new_infections_s) + sum(new_infections_ps)) > P) data_are_consistent <- FALSE
    inc_dat <- data.frame(t=1:(t-1), inc=new_infections_s[1:(t-1)] + new_infections_ps[1:(t-1)],
                          inc_s=new_infections_s[1:(t-1)], inc_ps = new_infections_ps[1:(t-1)],
                          para=paralysis_incidence_s[1:(t-1)]+paralysis_incidence_ps[1:(t-1)],
                          para_s=paralysis_incidence_s[1:(t-1)], para_ps=paralysis_incidence_ps[1:(t-1)],
                          Rt=Rt[1:(t-1)])   
    ## Store final conditions
    final_conditions <- data.frame(incidence_s=new_infections_s[1:tmax1], incidence_ps=new_infections_ps[1:tmax1],
                                   fully_susceptible=fully_susceptible[1:tmax1], partially_susceptible=partially_susceptible[1:tmax1],
                                   paralysis_s=paralysis_incidence_s[1:tmax1], paralysis_ps=paralysis_incidence_ps[1:tmax1],Rt=Rt[1:tmax1])
    return(list(dat=inc_dat,
                final_conditions=final_conditions,
                t_end=t,
                n_paralysis=sum(paralysis_incidence_s[1:(t-1)] + paralysis_incidence_ps[1:(t-1)]),
                data_are_consistent=data_are_consistent))
}


## Simulate from assumed priors
simulate_priors <- function(n=100,  
                           incu_fix=FALSE,
                           incu_mean_prior_mean=16,incu_mean_prior_var=5,
                           incu_var_prior_mean=10,incu_var_prior_var=3,
                           
                           ## Generation interval parameters
                           generation_fix=FALSE,
                           latent_period=3,
                           infect_mean_prior_mean=7,infect_mean_prior_var=5,
                           infect_var_prior_mean=5,infect_var_prior_var=3,
                           
                           infect_ps_par1_prior_mean=0.27,infect_ps_par1_prior_var=0.1,
                           #infect_ps_mean_prior_mean=infect_mean_prior_mean,infect_ps_mean_prior_var=infect_mean_prior_var,
                           #infect_ps_var_prior_mean=infect_var_prior_mean,infect_ps_var_prior_var=infect_var_prior_var,
                           
                           ## R0 draws, can choose distribution
                           R0_par1=0,R0_par2=10,R0_dist="uniform",
                           
                           ## Relative R0 draws
                           rel_R0_fixed=FALSE, rel_R0_mean=1,rel_R0_var=0.001,
                           
                           ## Prob paralysis parameters
                           prob_paralysis_fixed=FALSE,
                           prob_paralysis_mean=1/2000,prob_paralysis_var=1e-8,
                           prob_paralysis_ps_mean=prob_paralysis_mean,prob_paralysis_ps_var=prob_paralysis_var,
                           
                           ## Immunity groups
                           prop_immune_fixed=FALSE,
                           prop_immune_groups=c(0,0,1), prop_immune_pars=c(5,5,5)){
    ## Incubation period
    if(!incu_fix){
        incu_mean_scale <- find_gamma_pars(incu_mean_prior_mean,incu_mean_prior_var)[1]
        incu_mean_shape <- find_gamma_pars(incu_mean_prior_mean,incu_mean_prior_var)[2]
        
        incu_var_scale <- find_gamma_pars(incu_var_prior_mean,incu_var_prior_var)[1]
        incu_var_shape <- find_gamma_pars(incu_var_prior_mean,incu_var_prior_var)[2]
        
        ## Simulate means and variances of incubation period then convert to gamma parameters
        incu_mean <- rgamma(n, scale=incu_mean_scale, shape=incu_mean_shape)
        incu_var <- rgamma(n, scale=incu_var_scale, shape=incu_var_shape)
    } else {
        incu_mean <- rep(incu_mean_prior_mean, n)
        incu_var <- rep(incu_var_prior_mean, n)
    }
    incu_scale <- find_gamma_pars(incu_mean,incu_var)[[1]]
    incu_shape <- find_gamma_pars(incu_mean,incu_var)[[2]]
    
    ## Infectious period, fully susceptible
    if(!generation_fix){
        infect_mean_scale <- find_gamma_pars(infect_mean_prior_mean,infect_mean_prior_var)[1]
        infect_mean_shape <- find_gamma_pars(infect_mean_prior_mean,infect_mean_prior_var)[2]
        
        infect_var_scale <- find_gamma_pars(infect_var_prior_mean,infect_var_prior_var)[1]
        infect_var_shape <- find_gamma_pars(infect_var_prior_mean,infect_var_prior_var)[2]
        
        infect_mean <- rgamma(n, scale=infect_mean_scale, shape=infect_mean_shape)
        infect_var <- rgamma(n, scale=infect_var_scale, shape=infect_var_shape)
        
        ## Infectious period, partially immune
        infect_ps_par <- rlnorm(n, log(infect_ps_par1_prior_mean), infect_ps_par1_prior_var)
        
        #infect_ps_mean_scale <- find_gamma_pars(infect_ps_mean_prior_mean,infect_ps_mean_prior_var)[1]
        #infect_ps_mean_shape <- find_gamma_pars(infect_ps_mean_prior_mean,infect_ps_mean_prior_var)[2]
        
        #infect_ps_var_scale <- find_gamma_pars(infect_ps_var_prior_mean,infect_ps_var_prior_var)[1]
        #infect_ps_var_shape <- find_gamma_pars(infect_ps_var_prior_mean,infect_ps_var_prior_var)[2]
        
        #infect_ps_mean <- rgamma(n, scale=infect_ps_mean_scale, shape=infect_ps_mean_shape)
        #infect_ps_var <- rgamma(n, scale=infect_ps_var_scale, shape=infect_ps_var_shape)
    } else {
        infect_mean <- rep(infect_mean_prior_mean, n)
        infect_var <- rep(infect_var_prior_mean, n)
        infect_ps_par <- rep(infect_ps_par1_prior_mean, n)
        #infect_ps_mean <- rep(infect_ps_mean_prior_mean, n)
        #infect_ps_var <- rep(infect_ps_var_prior_mean, n)
    }
    
    infect_scale <- find_gamma_pars(infect_mean,infect_var)[[1]]
    infect_shape <- find_gamma_pars(infect_mean,infect_var)[[2]]
    
    ## If the priors are the same, assume we want identical generation intervals
    #if(infect_ps_mean_prior_mean == infect_mean_prior_mean &
    #    infect_ps_mean_prior_var == infect_mean_prior_var &
    #    infect_ps_var_prior_mean == infect_var_prior_mean &
    #    infect_ps_var_prior_var == infect_var_prior_var) {
    #    infect_ps_scale <- infect_scale
    #    infect_ps_shape <- infect_shape
    #} else {
    #    infect_ps_scale <- find_gamma_pars(infect_ps_var,infect_ps_mean)[[1]]
    #    infect_ps_shape <- find_gamma_pars(infect_ps_var,infect_ps_mean)[[2]]
    #}
    
    ## Simulate R0s
    if(R0_dist == "gamma"){
        R0_scale <- R0_par2/R0_par1
        R0_shape <- R0_par1/R0_scale
        R0 <- rgamma(n, scale=R0_scale,shape=R0_shape)
    } else if(R0_dist == "lognormal"){
        R0 <- rlnorm(n, R0_par1,R0_par2)
    } else if(R0_dist == "truncnorm"){
        R0 <- truncnorm::rtruncnorm(n, a=0,b=10,R0_par1,R0_par2)
    } else {
        R0 <- runif(n, R0_par1, R0_par2)
    }
    
    ## Simulate relative R0s
    if(!rel_R0_fixed) {
        rel_R0_pars <- get_beta_pars(rel_R0_mean, rel_R0_var)
        rel_R0s <- rbeta(n, rel_R0_pars[1], rel_R0_pars[2])
    } else {
        rel_R0s <- rep(rel_R0_mean, n)
    }
    
    ## Simulate paralysis probabilities
    if(!prob_paralysis_fixed){
        para_pars <- get_beta_pars(prob_paralysis_mean, prob_paralysis_var)
        para_pars_ps <- get_beta_pars(prob_paralysis_ps_mean, prob_paralysis_ps_var)
        prob_para <- rbeta(n, para_pars[1], para_pars[2])
        prob_para_ps <- rbeta(n, para_pars_ps[1], para_pars_ps[2])
    } else {
        prob_para <- rep(prob_paralysis_mean, n)
        prob_para_ps <- rep(prob_paralysis_ps_mean, n)
    }
    
    ## Simulate proportions in fully immune, partially immune and fully susceptible
    if(!prop_immune_fixed){
        prop_immune_groups <- LaplacesDemon::rdirichlet(n, prop_immune_pars)
    } else {
        prop_immune_groups <- matrix(rep(prop_immune_groups,n), ncol=3,byrow=TRUE)
    }
    
    main_pars <- data.frame(R0=R0, rel_R0s=rel_R0s,
                            latent_period=latent_period,
                            infect_scale=infect_scale,infect_shape=infect_shape,
                            infect_ps_par=infect_ps_par,
                            infect_mean=infect_mean,
                            infect_var=infect_var,
                            #infect_scale_ps=infect_ps_scale[i],infect_shape_ps=infect_ps_shape[i],
                            incu_scale=incu_scale,incu_shape=incu_shape,
                            prob_paralysis_s = prob_para, prob_paralysis_ps = prob_para_ps,
                            prop_immune_groups = prop_immune_groups)
    return(main_pars)
    
}


random_simulation_twoimmune <- function( n=100,  observed_data=c(1,rep(0,22)),
                               tmax=100,P=350000,continue_run=FALSE,
                               max_infectious_period=50,
                               
                               incu_fix=FALSE,
                               incu_mean_prior_mean=16,incu_mean_prior_var=5,
                               incu_var_prior_mean=10,incu_var_prior_var=3,
                               
                               ## Generation interval parameters
                               generation_fix=FALSE,
                               latent_period=3,
                               infect_mean_prior_mean=7,infect_mean_prior_var=5,
                               infect_var_prior_mean=5,infect_var_prior_var=3,
                               
                               infect_ps_par1_prior_mean=0.27,infect_ps_par1_prior_var=0.1,
                               #infect_ps_mean_prior_mean=infect_mean_prior_mean,infect_ps_mean_prior_var=infect_mean_prior_var,
                               #infect_ps_var_prior_mean=infect_var_prior_mean,infect_ps_var_prior_var=infect_var_prior_var,
                               
                               ## R0 draws, can choose distribution
                               R0_par1=0,R0_par2=10,R0_dist="uniform",
                               
                               ## Relative R0 draws
                               rel_R0_fixed=FALSE, rel_R0_mean=1,rel_R0_var=0.001,
                               
                               ## Prob paralysis parameters
                               prob_paralysis_fixed=FALSE,
                               prob_paralysis_mean=1/2000,prob_paralysis_var=1e-8,
                               prob_paralysis_ps_mean=prob_paralysis_mean,prob_paralysis_ps_var=prob_paralysis_var,
                               
                               ## Immunity groups
                               prop_immune_fixed=FALSE,
                               prop_immune_groups=c(0,0,1), prop_immune_pars=c(5,5,5),
                               
                               ## Seed size and simulation ID start
                               ini_infs=5, index_start=1
){
    
    pars <- simulate_priors(n=n,  
                            incu_fix=incu_fix,
                            incu_mean_prior_mean=incu_mean_prior_mean,incu_mean_prior_var=incu_mean_prior_var,
                            incu_var_prior_mean=incu_var_prior_mean,incu_var_prior_var=incu_var_prior_var,
                            
                            ## Generation interval parameters
                            generation_fix=generation_fix,
                            latent_period=latent_period,
                            infect_mean_prior_mean=infect_mean_prior_mean,infect_mean_prior_var=infect_mean_prior_var,
                            infect_var_prior_mean=infect_var_prior_mean,infect_var_prior_var=infect_var_prior_var,
                            
                            infect_ps_par1_prior_mean=infect_ps_par1_prior_mean,infect_ps_par1_prior_var=infect_ps_par1_prior_var,
                            #infect_ps_mean_prior_mean=infect_mean_prior_mean,infect_ps_mean_prior_var=infect_mean_prior_var,
                            #infect_ps_var_prior_mean=infect_var_prior_mean,infect_ps_var_prior_var=infect_var_prior_var,
                            
                            ## R0 draws, can choose distribution
                            R0_par1=R0_par1,R0_par2=R0_par2,R0_dist=R0_dist,
                            
                            ## Relative R0 draws
                            rel_R0_fixed=rel_R0_fixed, rel_R0_mean=rel_R0_mean,rel_R0_var=rel_R0_var,
                            
                            ## Prob paralysis parameters
                            prob_paralysis_fixed=prob_paralysis_fixed,
                            prob_paralysis_mean=prob_paralysis_mean,prob_paralysis_var=prob_paralysis_var,
                            prob_paralysis_ps_mean=prob_paralysis_ps_mean,prob_paralysis_ps_var=prob_paralysis_ps_var,
                            
                            ## Immunity groups
                            prop_immune_fixed=prop_immune_fixed,
                            prop_immune_groups=prop_immune_groups, prop_immune_pars=prop_immune_pars)
    ## Empty data structures to store outputs
    simulation_results <- NULL
    final_conditions <- NULL
    final_size <- numeric(n)
    
    incubation_periods <- NULL
    infectious_periods <- NULL
    infectious_ps_periods <- NULL

    tmax_vector <- numeric(n)
    n_paralysis <- numeric(n)
    data_consistent <- logical(n)
    
    for(i in 1:n){
        if(i %% 100 == 0) print(paste0("Sim number: ", i, " of ", n))
        tmp <- run_simulation_twoimmune(
            tmax=tmax,
            R0=pars$R0[i], rel_R0=pars$rel_R0s[i],P=P,
            observed_data=observed_data,
            latent_period=pars$latent_period[i],
            infect_scale=pars$infect_scale[i],infect_shape=pars$infect_shape[i],
            infect_ps_par=pars$infect_ps_par[i],
            #infect_scale_ps=infect_ps_scale[i],infect_shape_ps=infect_ps_shape[i],
            incu_scale=pars$incu_scale[i],incu_shape=pars$incu_shape[i],
            prob_paralysis_s = pars$prob_paralysis_s[i], prob_paralysis_ps = pars$prob_paralysis_ps[i],
            prop_immune_groups = as.numeric(pars[i,c("prop_immune_groups.1","prop_immune_groups.2","prop_immune_groups.3")]))
        
        sim_no <- i + index_start - 1
        
        dat <- tmp$dat
        dat$sim <- sim_no
        simulation_results[[i]] <- dat
        
        end <- tmp$final_conditions
        end$sim <- sim_no
        final_conditions[[i]] <- end
        
        final_size[i] <- sum(dat$inc)
        n_paralysis[i] <- sum(tmp$n_paralysis)
        
        infectious_periods[[i]] <- data.frame(sim=sim_no, t=1:max_infectious_period - 1,prob=extraDistr::ddgamma(1:max_infectious_period - pars$latent_period[i], shape=pars$infect_shape[i],scale=pars$infect_scale[i]))
        
        #infectious_ps_periods[[i]] <- data.frame(sim=sim_no, t=1:max_infectious_period - 1,prob=extraDistr::ddgamma(1:max_infectious_period, shape=infect_ps_shape[i],scale=infect_ps_scale[i]))
        infectious_ps_periods[[i]] <- data.frame(sim=sim_no, t=1:max_infectious_period - 1,prob=dexp(1:max_infectious_period - pars$latent_period[i], pars$infect_ps_par[i]))
        
        incubation_periods[[i]] <- data.frame(sim=sim_no, t=1:50 - 1, prob=extraDistr::ddgamma(1:50, shape=pars$incu_shape[i],scale=pars$incu_scale[i]))
        
        
        data_consistent[i] <- tmp$data_are_consistent
        
        ## Store final conditions for restart simulations
        tmax_vector[i] <- tmp$t_end
        
    }
    par_table <- data.frame(sim=index_start:(index_start + n - 1),
                            R0=pars$R0,rel_R0=pars$rel_R0s,
                            infectious_period_mean=pars$infect_mean,
                            infectious_period_var=pars$infect_var,
                            infect_shape=pars$infect_shape,infect_scale=pars$infect_scale,
                            infect_ps_par=pars$infect_ps_par,
                            #infectious_period_mean_ps=infect_ps_mean,
                            #infectious_period_var_ps=infect_ps_var,
                            #infect_shape_ps=infect_ps_shape,infect_scale_ps=infect_ps_scale,
                            latent_period=pars$latent_period,
                            incu_shape=pars$incu_shape,incu_scale=pars$incu_scale,
                            prob_paralysis_s = pars$prob_paralysis_s, 
                            prob_paralysis_ps = pars$prob_paralysis_ps,
                            prop_immune_groups = unname(pars[,c("prop_immune_groups.1","prop_immune_groups.2","prop_immune_groups.3")]))
    simulation_results <- do.call("bind_rows",simulation_results)
    infectious_periods <- do.call("bind_rows",infectious_periods)
    incubation_periods <- do.call("bind_rows",incubation_periods)
    infectious_ps_periods <- do.call("bind_rows",infectious_ps_periods)
    
    list(simulation_results=simulation_results, 
         final_conditions=final_conditions,
         pars=par_table,
         final_size=final_size,
         infectious_periods=infectious_periods,
         infectious_ps_periods=infectious_ps_periods,
         incubation_periods=incubation_periods,
         tmax_vector = tmax_vector,
         n_paralysis=data.frame(sim=index_start:(index_start+n-1), paralysis=n_paralysis),
         data_are_consistent=data_consistent)
}

restart_simulations_table_twoimmune <- function(use_sims,pars,
                                                final_conditions,
                                                t_starts,vaccinate_proportion=matrix(c(0,0),ncol=2),
                                                tmax,nruns=1,P=350000){
    
    
    match_sims <- match(use_sims,pars$sim)
    
    res_all <- NULL
    
    for(i in seq_along(match_sims)){
        res <- NULL
        
        index <- match_sims[i]
        print(paste0("Extending: ",index))
        x <- 1
        for(k in 1:nrow(vaccinate_proportion)){
            for(j in 1:nruns){
                tmp <- run_simulation_twoimmune(R0 = pars$R0[index], rel_R0 = pars$rel_R0[index], 
                                                incu_shape=pars$incu_shape[index],
                                                latent_period=pars$latent_period[index],
                                                infect_scale = pars$infect_scale[index], infect_shape = pars$infect_shape[index], 
                                                #infect_scale_ps = pars$infect_scale_ps[index],infect_shape_ps=pars$infect_shape_ps[index], 
                                                infect_ps_par=pars$infect_ps_par[index],
                                                tmax = tmax,
                                                P = P, 
                                                prob_paralysis_s = pars$prob_paralysis_s[index], 
                                                prob_paralysis_ps = pars$prob_paralysis_ps[index], 
                                                prop_immune_groups = as.numeric(pars[index,c("prop_immune_groups.1","prop_immune_groups.2","prop_immune_groups.3")]),
                                                       
                                                continue_run = TRUE, restart_simulation = TRUE, 
                                                final_conditions = final_conditions[[index]], t_start = t_starts[index], vaccinate_proportion = vaccinate_proportion[k,])
                
                dat <- tmp$dat
                dat$sim <- use_sims[i]
                dat$rep <- j
                dat$vacc_prop <- k
                res[[x]] <- dat
                
                x <- x + 1
            }
        }
        res_all[[i]] <- do.call("bind_rows",res)
    }
    
    res_all <- do.call("bind_rows",res_all)
    return(res_all)
}