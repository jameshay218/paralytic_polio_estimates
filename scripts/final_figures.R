library(ggsci)
load("outputs/fitted_trajectories.RData")

## Pad prematurely ended trajectories with zeros
N_tot <- length(unique(trajectories$sim))
## Pad up to this date
end_date <- max(trajectories$t)
## Create dummy rows from this date onward, unique to each sim
max_t1 <- trajectories %>% group_by(sim) %>% filter(t==max(t)) %>% select(sim, t, date_start) %>% mutate(t = t + 1) %>% rename(max_t=t) %>% mutate(max_date=end_date)

## Create the dummy rows
extra_rows <- max_t1 %>%
    filter(max_t < max_date) %>%
    group_by(sim) %>%
    distinct() %>%
    do(data.frame(sim=.$sim, t=seq(.$max_t,.$max_date,by="1 day")))
extra_rows <- extra_rows %>% mutate(inc=0,inc_s=0,inc_ps=0,para=0,para_s=0,para_ps=0, Rt=0)

## Add the dummy rows in
extra_rows <- extra_rows %>% left_join(trajectories %>% select(sim, date_start) %>% distinct())
trajectories <- bind_rows(trajectories,extra_rows) %>% arrange(sim, t)
trajectories <- trajectories %>% mutate(Rt = ifelse(is.na(Rt),0, Rt))

## For every day after 2022-08-20, if we've seen another 0, 1, 2 or 3 
## cases of paralysis, flag
max_date <- "2023-03-01"

flag_traj <- trajectories %>% group_by(sim) %>% 
    filter(t > "2022-08-20", t <= max_date) %>%
    group_by(sim) %>% 
    mutate(cumu_para=cumsum(para)) %>%
    ungroup() %>%
    ## Flag when we've seen 1, 2 or 3 total cases after 20th August
    mutate(
            no_new_para = cumu_para == 0, ## Have still no cases
            one_new_para=cumu_para == 1, ## I have only seen one case
           two_new_para=cumu_para == 2, ## I have only seen two cases
           three_new_para=cumu_para == 3 ## I have only seen three cases
           ) %>%
    select(sim, t, no_new_para, one_new_para,two_new_para, three_new_para)


## Find total paralysis cases by 31st December 2022
traj_cumu_para <- trajectories %>% 
    group_by(sim) %>% mutate(para=cumsum(para)) %>% 
    filter(t == max_date) %>%
    select(sim, para) %>%
    rename(total_para=para)
    
tmp <- flag_traj %>% left_join(traj_cumu_para)


omg1 <- tmp %>% group_by(sim) %>% filter(no_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_no_cases = t) %>%
    select(sim, last_date_no_cases,total_para)

omg2 <- tmp %>% group_by(sim) %>% filter(one_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_one_case = t) %>%
    select(sim, last_date_one_case,total_para)

omg3 <- tmp %>% group_by(sim) %>% filter(two_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_two_cases = t) %>%
    select(sim, last_date_two_cases,total_para)

omg4 <- tmp %>% group_by(sim) %>% filter(three_new_para == TRUE) %>%
    filter(t == max(t)) %>%
    mutate(last_date_three_cases = t) %>%
    select(sim, last_date_three_cases,total_para)

tmp_all <- omg1 %>% left_join(omg2) %>% left_join(omg3) %>% left_join(omg4) %>% ungroup()

dates <- seq(as.Date("2022-08-20"),as.Date(max_date),by="1 day")

ns <- numeric(length(dates))
subset_flags1 <- tmp_all
subset_flags2 <- tmp_all
subset_flags3 <- tmp_all
subset_flags4 <- tmp_all

all_no_para <- NULL
all_one_para <- NULL
all_two_para <- NULL
all_three_para <- NULL
quantile_lower <- 0.1
quantile_upper <- 0.9
for(i in seq_along(dates)){
    print(i)
    subset_flags1 <- subset_flags1 %>% filter(last_date_no_cases >= dates[i])
    subset_flags2 <- subset_flags2 %>% filter(last_date_one_case >= dates[i])
    subset_flags3 <- subset_flags3 %>% filter(last_date_two_cases >= dates[i])
    subset_flags4 <- subset_flags4 %>% filter(last_date_three_cases >= dates[i])
    ns[i] <- nrow(subset_flags)
    all_no_para[[i]] <- subset_flags1 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="0", t =dates[i])
    
    all_one_para[[i]] <- subset_flags2 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="1", t =dates[i])
    
    all_two_para[[i]] <- subset_flags3 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="2", t =dates[i])
    
    all_three_para[[i]] <- subset_flags4 %>% 
        dplyr::summarize(median_para=median(total_para),
                         mean_para=mean(total_para),
                         lower=quantile(total_para,quantile_lower),
                  upper=quantile(total_para,quantile_upper)) %>%
        mutate(model="3", t =dates[i])
    
}
all_no_para <- do.call("bind_rows",all_no_para)
all_one_para <- do.call("bind_rows",all_one_para)
all_two_para <- do.call("bind_rows",all_two_para)
all_three_para <- do.call("bind_rows",all_three_para)

tmp_comb <- bind_rows(all_no_para,all_one_para, all_two_para, all_three_para)
tmp_comb %>%
    mutate(model = as.factor(model)) %>%
    mutate(model=paste0("Further paralysis cases observed: ", model)) %>%
    ggplot() + 
    geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=model),alpha=0.05) +
    geom_line(aes(x=t,y=upper,col=model,linetype="90% quantiles"),size=0.25) +
    geom_line(aes(x=t,y=lower,col=model,linetype="90% quantiles"),size=0.25) +
    #geom_line(aes(x=t,y=upper,col=model),size=0.25) +
    #geom_line(aes(x=t,y=lower,col=model),size=0.25) +
    #geom_line(aes(x=t,y=lower,col=model),linetype="dotted") +
    #geom_line(aes(x=t,y=mean_para,col=model,linetype="Mean"),size=1) +
    
    geom_line(aes(x=t,y=median_para,col=model,linetype="Median"),size=0.8) +
    scale_y_continuous(expand=c(0,0),
                       limits=c(0,25)) +
    scale_fill_nejm(name="Additional paralysis\ncases observed") +
    scale_color_nejm(name="Additional paralysis\ncases observed") +
    scale_linetype_manual(name="",values=c("90% quantiles"="dashed",
                                    "Median"="solid")) +
    ylab(paste0("Projected incidence of paralysis by ",max_date)) +
    #scale_x_date(breaks="1 month") +
    #scale_x_date(breaks="1 month") +
    xlab("Date") +
    theme_classic() +
    theme(panel.grid.major = element_line(size=0.25,color="grey70"),
          panel.grid.minor = element_line(size=0.25,color="grey70"),
          legend.position="none") +
    facet_wrap(~model)

