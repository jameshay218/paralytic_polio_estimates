## We have a distribution of ages in the population
age_dist <- c(0.1,0.2,0.25,0.15,0.15,0.1,0.05)
ages <- seq_along(age_dist)

p_age <- ggplot(data.frame(prop=age_dist,Age=as.factor(ages))) +
    geom_bar(aes(x=Age,y=prop),stat="identity") +
    ylab("Proportion") +
    ggtitle("Age distribution of population")

## Vaccine groups:
## 1. OPV
## 2. IPV
## 3. Unvaccinated

## Assume that some proportion of each age group are in one of 3 vaccine groups:
## OPV, IPV and unvaccinated
prop_vacc_by_age <- matrix(c(
                        0.00,0.60,0.40,
                        0.00,0.80,0.20,
                        0.00,0.95,0.05,
                        0.00,0.95,0.05,
                        0.80,0.15,0.05,
                        0.80,0.10,0.10,
                        0.80,0.00,0.20
                        ),ncol=3,byrow=TRUE)

p_vacc_age <- prop_vacc_by_age %>% reshape2::melt() %>%
    rename(`Age group`=Var1, `Vaccine group`=Var2) %>%
    mutate(`Age group` = as.factor(`Age group`)) %>%
    ggplot() + geom_tile(aes(x=`Age group`,y=`Vaccine group`,fill=value))+
    scale_fill_viridis_c(name="Proportion") +
    ggtitle("Proportion of each age group in each vaccine group")

## Age represents time-since-vaccination, assuming that everyone who gets 
## vaccinated gets vaccinated in the first age class. 
## Therefore we can re-scale the proportion in each vaccine class by the proportion
## in each time-since-vaccination group
t_since_vacc_dist <- rev(age_dist)
prop_vacc_by_age_scaled <- age_dist * prop_vacc_by_age
prop_vacc_by_tsv_scaled <- prop_vacc_by_age_scaled[nrow(prop_vacc_by_age_scaled):1,]

p_vacc_age_scaled <- prop_vacc_by_age_scaled %>% as_tibble() %>% 
    mutate(age_group=1:n()) %>%
    pivot_longer(-age_group) %>%
    rename(vaccine=name) %>%
    ggplot() + geom_tile(aes(x=age_group,y=rev(vaccine),fill=value)) +
    scale_fill_viridis_c(name="Proportion") +
    scale_x_continuous(breaks=seq(0,length(ages))) +
    xlab("Age group") + ylab("Vaccine class") +
    ggtitle("Proportion of population in each vaccine\n class by age group")

## Assume that depending on time-since-vaccination, some proportion 
## end up in one of three immune classes: fully immune, partially immune, fully susceptible
## Define this separately for each of the 3 vaccine classes Note that 
vacc_class_immune_states_time <- 
    array(c(
        1.00,0.90,0.80,0.60,0.50,0.50,0.50,
        0.00,0.10,0.20,0.35,0.40,0.30,0.30,
        0.00,0.00,0.00,0.05,0.10,0.20,0.20,
        
        1.00,0.80,0.60,0.40,0.30,0.25,0.25,
        0.00,0.20,0.40,0.55,0.50,0.50,0.50, 
        0.00,0.00,0.00,0.05,0.20,0.25,0.25,
            
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,
        0.00,0.00,0.00,0.00,0.00,0.00,0.00,
        1.00,1.00,1.00,1.00,1.00,1.00,1.00
        ), dim=c(length(time_since_vacc),3,3)
           )

vacc_class_immune_states_time_melted <- reshape2::melt(vacc_class_immune_states_time)
vacc_class_immune_states_time_melted <- vacc_class_immune_states_time_melted %>% 
    rename(t_since_vacc=Var1,immune_class=Var2,vaccine_class=Var3)
p_vacc_immune_time <- ggplot(vacc_class_immune_states_time_melted) + 
    geom_tile(aes(x=t_since_vacc,y=immune_class,fill=value)) + 
    scale_x_continuous(breaks=seq(1,length(ages))) +
    facet_wrap(~vaccine_class,ncol=1) +
    scale_fill_viridis_c(name="Proportion") +
    xlab("Time since vaccination") +
    ylab("Immune class") +
    ggtitle("Proportion of individuals in each immunity class\n conditional on time since vaccination and vaccine type")

p_waning <- vacc_class_immune_states_time_melted %>% ggplot() + geom_bar(aes(x=t_since_vacc,y=value,fill=as.factor(immune_class)),stat="identity") + facet_wrap(~vaccine_class,ncol=1) + scale_fill_viridis_d(name="Immune class") + ylab("Proportion") + xlab("Time since vaccination") + scale_x_continuous(breaks=seq(1,length(ages)))

## We then re-scale everything by the proportion in each 
## vaccination and time-since-vaccination class
x <- prop_vacc_by_tsv_scaled[,3] %*% vacc_class_immune_states_time[,,3] +
    prop_vacc_by_tsv_scaled[,2] %*% vacc_class_immune_states_time[,,2] +
    prop_vacc_by_tsv_scaled[,1] %*% vacc_class_immune_states_time[,,1]
x

p_age
p_vacc_age
p_vacc_age_scaled
p_vacc_immune_time
p_waning
