library(ggsci)
library(Hmisc)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(lubridate)
library(tidyverse)
library(paletteer)
library(data.table)

setwd("~/Documents/GitHub/paralytic_polio_estimates/")
source("simulation_functions_twoimmune.R")
max_date <- "2023-04-01"
summarize <- dplyr::summarize
setwd("~/Documents/GitHub/paralytic_polio_estimates/")

scales::show_col(pal_nejm("default")(8))
nejm_palette <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF")

reload <- FALSE
# Read in simulation runs and combine -------------------------------------
if(reload){
    setwd("sims/")
    all_res <-  NULL
    for(i in 1:500){
        if(file.exists(paste0("simulation_",i,".RData"))){
            load(paste0("simulation_",i,".RData"))
            all_res[[i]] <- pars
        }
    }
    res <- as_tibble(do.call("bind_rows",all_res))
    res <- res %>% select(-c(Rt, inc_s,inc_ps,para,para_s,para_ps))
    
    
    setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_traj/")
    all_traj <-  NULL
    for(i in 1:500){
        if(file.exists(paste0("traj_",i,".RData"))){
            load(paste0("traj_",i,".RData"))
            all_traj[[i]] <- res2
        }
    }
    traj <- as_tibble(do.call("bind_rows",all_traj))
    
    setwd("~/Documents/GitHub/paralytic_polio_estimates/sims_nyc/")
    all_nyc <-  NULL
    for(i in 1:500){
        if(file.exists(paste0("nyc_",i,".RData"))){
            load(paste0("nyc_",i,".RData"))
            all_nyc[[i]] <- nyc
        }
    }
    traj_nyc <- as_tibble(do.call("bind_rows",all_nyc))
    
    
    ## Set day 0
    traj <- traj %>% left_join(res %>% select(sim, date_start)) %>%
        mutate(t = as.Date(date_start + t)) %>%
        drop_na() 
    
    ## Get cumulative incidence
    traj <- traj %>% group_by(sim) %>% mutate(cumu_inc = cumsum(inc),
                                              cumu_para = cumsum(para))
    #traj <- traj %>% ungroup() %>% mutate(month = round_date(t, "month"))
    
    ## Assume NYC is seeded 4 weeks later
    traj_nyc <- traj_nyc %>% 
        left_join(res %>% select(sim, date_start)) %>%
        mutate(t = as.Date(date_start + t) + 28) %>%
        drop_na() 
    
    traj_nyc <- traj_nyc %>% group_by(sim) %>% mutate(cumu_inc = cumsum(inc),
                                              cumu_para = cumsum(para))
    #traj_nyc <- traj_nyc %>% ungroup() %>% mutate(month = round_date(t, "month"))
    
    traj$vacc_prop <- as.factor(traj$vacc_prop)
    trajectories <- traj %>% filter(vacc_prop==1) %>% dplyr::select(-c(vacc_prop,rep))
    
    subset_sims <- sample(unique(traj$sim),min(10000,length(unique(traj$sim))))
    trajectories <- trajectories %>% filter(sim %in% subset_sims)
    traj_nyc <- traj_nyc %>% filter(sim %in% subset_sims)
    res <- res %>% filter(sim %in% subset_sims)
    
    
    save(trajectories,file="../outputs/fitted_trajectories.RData")
    save(traj_nyc,file="../outputs/nyc_trajectories.RData")
    save(res,file="../outputs/filtered_parameters.RData")
} else {
    load("outputs/fitted_trajectories.RData")
    load("outputs/nyc_trajectories.RData")
    load("outputs/filtered_parameters.RData")
}

# Clean trajectories ------------------------------------------------------
## Rockland County
pad_rows <- function(tmp_traj){
    ## Pad prematurely ended trajectories with zeros
    N_tot <- length(unique(tmp_traj$sim))
    ## Pad up to this date
    end_date <- max(tmp_traj$t)
    ## Create dummy rows from this date onward, unique to each sim
    max_t1 <- tmp_traj %>% group_by(sim) %>% filter(t==max(t)) %>% 
        select(sim, t, date_start) %>% mutate(t = t + 1) %>% 
        rename(max_t=t) %>% mutate(max_date=end_date)
    
    ## Create the dummy rows
    extra_rows <- max_t1 %>%
        filter(max_t < max_date) %>%
        group_by(sim) %>%
        distinct() %>%
        do(data.frame(sim=.$sim, t=seq(.$max_t,.$max_date,by="1 day")))
    extra_rows <- extra_rows %>% mutate(inc=0,inc_s=0,inc_ps=0,para=0,para_s=0,para_ps=0, Rt=0)
    
    ## Add the dummy rows in
    extra_rows <- extra_rows %>% left_join(tmp_traj %>% select(sim, date_start) %>% distinct())
    tmp_traj <- bind_rows(tmp_traj,extra_rows) %>% arrange(sim, t)
    tmp_traj <- tmp_traj %>% mutate(Rt = ifelse(is.na(Rt),0, Rt))
    tmp_traj <- tmp_traj %>% group_by(sim) %>% mutate(cumu_para=cumsum(para)) %>% ungroup()
    tmp_traj <- tmp_traj %>% group_by(sim) %>% 
        ## If cases aren't what they were 7 days ago, then the epidemic is still ongoing
        mutate(cumu_inc=cumsum(inc)) %>% 
        mutate(ongoing_7= cumu_inc != lag(cumu_inc,7,nafill(0))) %>%
        ungroup()
    tmp_traj
}
trajectories <- pad_rows(trajectories)
traj_nyc <- pad_rows(traj_nyc)

## Flag if further cases observed or not
tmp_flags <- trajectories %>% filter(t <= "2022-10-01") %>% group_by(sim) %>% dplyr::summarise(y=sum(para)) %>% mutate(no_further_cases=y==1)

## Summarize proportion of trajectories with ongoing incidence over time
traj_summary1 <- trajectories %>% filter(t <= max_date) %>% 
    mutate(y=ongoing_7 == TRUE) %>%
    group_by(t) %>%
    dplyr::summarize(
        median_para=median(cumu_para),lower1=quantile(cumu_para,0.1),upper1=quantile(cumu_para,0.9),
        y=sum(y),N=n(),prop=sum(y))%>% group_by(t) %>% 
    mutate(prop=y/N,lower=binconf(y,N)[2],upper=binconf(y,N)[3]) %>%
    mutate(model="Current data")

## Calculate extinction probability if no further cases observed
traj_summary2 <- trajectories %>% left_join(tmp_flags) %>%
    filter(no_further_cases==TRUE) %>%
    filter(t <= max_date) %>% 
    mutate(y=ongoing_7 == TRUE) %>%
    group_by(t) %>% 
    dplyr::summarize(
        median_para=median(cumu_para),lower1=quantile(cumu_para,0.1),upper1=quantile(cumu_para,0.9),
        y=sum(y),N=n())%>% group_by(t) %>% 
    mutate(prop=y/N,lower=binconf(y,N)[2],upper=binconf(y,N)[3])%>%
    mutate(model="No cases by October 1st")

## Calculate extinction probability if further cases observed
traj_summary3 <- trajectories %>% left_join(tmp_flags) %>%
    filter(no_further_cases==FALSE) %>%
    filter(t <= max_date) %>% 
    mutate(y=ongoing_7 == TRUE) %>%
    group_by(t) %>% 
    dplyr::summarize(
        median_para=median(cumu_para),lower1=quantile(cumu_para,0.1),upper1=quantile(cumu_para,0.9),
                     y=sum(y),N=n())%>% group_by(t) %>% 
    mutate(prop=y/N,lower=binconf(y,N)[2],upper=binconf(y,N)[3])%>%
    mutate(model="Further cases reported\n by October 1st")

## Combine
traj_summary <- bind_rows(traj_summary1,traj_summary2,traj_summary3)

## Plot some trajectories
set.seed(1234)
sub_sims1 <- trajectories %>% left_join(tmp_flags) %>% filter(no_further_cases==TRUE) %>% 
    select(sim) %>% distinct() %>% sample_n(15) %>% pull(sim)
sub_sims2 <- trajectories %>% left_join(tmp_flags) %>% filter(no_further_cases==FALSE) %>% 
    select(sim) %>% distinct() %>% sample_n(15) %>% pull(sim)
trajectories <- trajectories %>% 
    left_join(tmp_flags) %>%
    mutate(`Consistent with further\nparalytic polio\n cases by 2022-10-21` = ifelse(no_further_cases==TRUE,"Yes","No")) 

p_traj <- trajectories %>% 
    filter(sim %in% c(sub_sims1,sub_sims2)) %>% 
    filter(t <= as.Date("2022-10-21")) %>%
    ggplot() +
    geom_rect(data=data.frame(xmin=as.Date("2022-08-20"),xmax=as.Date("2022-11-01"),ymin=0,ymax=6000),
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="black",alpha=0.1) +
    
    geom_rect(data=data.frame(xmin=as.Date("2022-05-01"),xmax=as.Date("2022-08-20"),ymin=4300,ymax=5300),
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=nejm_palette[6],alpha=0.25) +
    
    geom_line(aes(x=t,y=inc,group=sim),col=nejm_palette[3],size=0.25) +
    geom_segment(data=data.frame(x=as.Date(c("2022-06-22","2022-06-22","2022-08-20")),
                                 xend=as.Date(c("2022-06-22","2022-06-22","2022-08-20")),y=c(0,5300,0),yend=c(4300,6000,6000)),
                 aes(x=x,xend=xend,y=y,yend=yend),
                 linetype="dashed") +    
    geom_text(data=data.frame(x=as.Date(c("2022-06-22","2022-08-20")),
                              y=6350,label=c("First case of\n paralysis","Last observation")),aes(x=x,y=y,label=label),
              size=1.75) +
    geom_text(data=data.frame(x=as.Date("2022-07-01"),y=4750,lab="Positive wastewater samples\n from Rockland County"),
              aes(x=x,y=y,label=lab),size=1.75) +
    ylab("Incidence of polio infections")+
    xlab("Date") +
    scale_x_date(limits=as.Date(c("2022-03-01","2022-11-01")),date_labels="%b",breaks="month") +
    scale_y_continuous(limits=c(0,6600),expand=c(0,0),breaks=seq(0,6000,by=1000)) +
    theme_classic() +
    scale_color_manual(values=c("Yes"=nejm_palette[1],"No"=nejm_palette[2])) +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.title=element_text(size=6),
          legend.position="none",
          panel.grid.major = element_line(size=0.1,color="grey70")) +
    labs(tag="A")
p_traj


p_final_size <- trajectories %>% group_by(sim) %>%
    mutate(y=cumsum(inc)) %>% 
    filter(t == "2022-08-20") %>% 
    mutate(model = "Current data") %>%
    ggplot() + geom_histogram(aes(x=y,fill=model),binwidth=10000,
                                  color="black",fill=nejm_palette[3],alpha=0.25) +
    scale_x_continuous(breaks=seq(0,250000,by=100000))+
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    xlab("Estimated infections by 2022-08-20") + ylab("Simulation count") +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=6),
          legend.position="none",
          legend.text=element_text(size=6),
          legend.title=element_text(size=6)) +
    scale_fill_nejm()+
    labs(tag="C")

p_start <- trajectories %>% 
    group_by(sim) %>%
    filter(t == min(t)) %>%
    mutate(model = "Current data") %>%
    ggplot() + geom_density(aes(x=date_start,fill=model),color="black",fill="grey70") +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    xlab("Seed date") + ylab("Density") +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=6),
          legend.position="none",
          legend.text=element_text(size=6),
          legend.title=element_text(size=6)) +
    scale_fill_nejm() +
    labs(tag="B")
p_re <- trajectories %>% 
    group_by(sim) %>%
    filter(t == min(t)) %>%
    #mutate(model = "Current data") %>%
    ggplot() + geom_density(aes(x=Rt),color="black",fill=nejm_palette[3],alpha=0.25) +
    geom_vline(xintercept=1,linetype="dashed") +
    scale_y_continuous(expand=c(0,0)) +
    theme_classic() +
    xlab("Initial effective reproductive number") + ylab("Density") +
    theme(axis.text=element_text(size=6),
          axis.title=element_text(size=6),
          legend.position="none",
          legend.text=element_text(size=6),
          legend.title=element_text(size=6)) +
    scale_fill_nejm() +
    labs(tag="B")

pX <- (p_re/p_final_size)
p_top <- p_traj + pX + plot_layout(ncol=2,widths=c(3,1))

p1 <- ggplot(traj_summary) + 
    geom_ribbon(aes(x=t,ymin=upper,ymax=lower,fill=model),alpha=0.25) +
    geom_line(aes(x=t,y=prop,col=model)) +
    geom_text(data=data.frame(x=as.Date(c("2022-06-22","2022-08-20")),
                              y=1.12,label=c("First case of\n paralysis","Last observation")),
              aes(x=x,y=y,label=label),size=1.75) +
    scale_y_continuous(limits=c(0,1.15),breaks=seq(0,1,by=0.2)) +
    geom_segment(data=data.frame(x=as.Date(c("2022-06-22","2022-08-20")),
                                 xend=as.Date(c("2022-06-22","2022-08-20")),
                                 y=0,yend=1),
                 aes(x=x,xend=xend,y=y,yend=yend),
               linetype="dashed") +    ylab("Probability epidemic is extinct\n (no cases in previous 7 days)")+
    theme_classic() +
    scale_fill_manual(values=c("Current data"=nejm_palette[3],"No cases by October 1st"=nejm_palette[2],
                               "Further cases reported\n by October 1st"=nejm_palette[1])) +
    scale_color_manual(values=c("Current data"=nejm_palette[3],"No cases by October 1st"=nejm_palette[2],
                               "Further cases reported\n by October 1st"=nejm_palette[1]))  +
    
    scale_x_date(date_label="%b",limits=as.Date(c("2022-06-01",max_date)),
                 breaks="1 month") +    
    theme(legend.position="none",axis.text.x=element_blank(),
          axis.line.x=element_blank(),axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),axis.text.y=element_text(size=6),
          axis.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.title=element_text(size=6),
          panel.grid.major = element_line(size=0.1,color="grey70")) +
    labs(tag="D")

p2 <- ggplot(traj_summary) + 
    geom_ribbon(aes(x=t,ymin=lower1,ymax=upper1,fill=model),alpha=0.25) +
    geom_line(aes(x=t,y=median_para,col=model)) +
    geom_vline(xintercept=as.Date(c("2022-06-22","2022-08-20")),
linetype="dashed") +
    ylab("Cumulative paralysis cases") +
    scale_x_date(date_label="%b",limits=as.Date(c("2022-06-01",max_date)),
                 breaks="1 month") +
    theme_classic() +
    xlab("Date") +
    scale_fill_manual(name="Data assumption",values=c("Current data"=nejm_palette[3],"No cases by October 1st"=nejm_palette[2],
                               "Further cases reported\n by October 1st"=nejm_palette[1])) +
    scale_color_manual(name="Data assumption",values=c("Current data"=nejm_palette[3],"No cases by October 1st"=nejm_palette[2],
                                "Further cases reported\n by October 1st"=nejm_palette[1])) +
    theme(legend.position="bottom",
          panel.grid.major = element_line(size=0.1,color="grey70"),
          axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6),
          axis.title=element_text(size=8),
          legend.text=element_text(size=6),
          legend.title=element_text(size=6))+
    labs(tag="E")

fig1 <- p_top / p1 / p2
fig1


# Hypothetical further observations ---------------------------------------
## For every day after 2022-08-20, if we've seen another 0, 1, 2 or 3 
## cases of paralysis, flag
get_flag_traj <- function(tmp_traj){
    flag_traj <- tmp_traj %>% group_by(sim) %>% 
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
    ## Find total paralysis cases by end date
    traj_cumu_para <- tmp_traj %>% 
        group_by(sim) %>% mutate(para=cumsum(para)) %>% 
        filter(t == max_date) %>%
        select(sim, para) %>%
        rename(total_para=para)
    tmp <- flag_traj %>% left_join(traj_cumu_para)
    
    ## Flag last date with no new paralysis
    omg1 <- tmp %>% group_by(sim) %>% filter(no_new_para == TRUE) %>%
        filter(t == max(t)) %>%
        mutate(last_date_no_cases = t) %>%
        select(sim, last_date_no_cases,total_para)
    
    ## Flag last date where we only had one cumulative case of paralysis
    omg2 <- tmp %>% group_by(sim) %>% filter(one_new_para == TRUE) %>%
        filter(t == max(t)) %>%
        mutate(last_date_one_case = t) %>%
        select(sim, last_date_one_case,total_para)
    ## Flag last date where we only had two cumulative cases of paralysis
    omg3 <- tmp %>% group_by(sim) %>% filter(two_new_para == TRUE) %>%
        filter(t == max(t)) %>%
        mutate(last_date_two_cases = t) %>%
        select(sim, last_date_two_cases,total_para)
    ## Flag last date where we only had three cumulative cases of paralysis
    omg4 <- tmp %>% group_by(sim) %>% filter(three_new_para == TRUE) %>%
        filter(t == max(t)) %>%
        mutate(last_date_three_cases = t) %>%
        select(sim, last_date_three_cases,total_para)
    
    tmp_all <- omg1 %>% left_join(omg2) %>% left_join(omg3) %>% left_join(omg4) %>% ungroup()
    tmp_all
}
tmp_all <- get_flag_traj(trajectories)
tmp_all_nyc <- get_flag_traj(traj_nyc)

filter_traj_by_cumu_para <- function(tmp_traj, quantile_lower=0.1,quantile_upper=0.9){
    dates <- seq(as.Date("2022-08-20"),as.Date(max_date),by="1 day")
    
    ns <- numeric(length(dates))
    subset_flags1 <- tmp_traj
    subset_flags2 <- tmp_traj
    subset_flags3 <- tmp_traj
    subset_flags4 <- tmp_traj
    
    all_no_para <- NULL
    all_one_para <- NULL
    all_two_para <- NULL
    all_three_para <- NULL
    
    ## Go through each date and filter only trajectories which 
    for(i in seq_along(dates)){
        print(i)
        ## Filter only trajectories whose last date of "no cases" is today or later
        subset_flags1 <- subset_flags1 %>% filter(
            !is.na(last_date_no_cases) &
                last_date_no_cases >= dates[i])
        ## Filter only trajectories whose last date of "only one case" is today or later, etc
        subset_flags2 <- subset_flags2 %>% filter(
            !is.na(last_date_one_case) &last_date_one_case >= dates[i])
        subset_flags3 <- subset_flags3 %>% filter(
            !is.na(last_date_two_cases) &last_date_two_cases >= dates[i])
        subset_flags4 <- subset_flags4 %>% filter(
            !is.na(last_date_three_cases) &last_date_three_cases >= dates[i])
        
        ns[i] <- nrow(subset_flags1)
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
    tmp_comb
}
tmp_comb <- filter_traj_by_cumu_para(tmp_all)
tmp_comb_nyc <- filter_traj_by_cumu_para(tmp_all_nyc)

plot_fig2 <- function(tmp_comb){
    ymax <- max(tmp_comb$upper) + 5
    fig2A <- tmp_comb %>%
        mutate(model = as.factor(model)) %>%
        mutate(model=paste0("Further paralysis cases observed: ", model)) %>%
        ggplot() + 
        geom_ribbon(aes(x=t,ymin=lower,ymax=upper,fill=model),alpha=0.05) +
        geom_line(aes(x=t,y=upper,col=model,linetype="90% quantiles"),size=0.25) +
        geom_line(aes(x=t,y=lower,col=model,linetype="90% quantiles"),size=0.25) +
        geom_line(aes(x=t,y=median_para,col=model,linetype="Median"),size=0.8) +
        scale_y_continuous(expand=c(0,0),
                           limits=c(0,ymax)) +
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
              legend.position="none",
              strip.background = element_blank(),
              axis.text=element_text(size=6),
              axis.title=element_text(size=8),
              strip.text=element_text(size=6)) +
        facet_wrap(~model,nrow=1)

    date_key <- c("2022-09-01"="Sep","2022-10-01"="Oct","2022-11-01"="Nov","2022-12-01"="Dec")
    tmp_comb$date1 <- date_key[as.character(tmp_comb$t)]
    
    tmp_comb1 <- tmp_comb %>% filter(t %in% as.Date(c("2022-09-01","2022-10-01","2022-11-01","2022-12-01"))) %>%
        mutate(date2=paste0("Cases observed by ", date1))
    tmp_comb1$date2 <- factor(tmp_comb1$date1,levels=c("Cases observed by Sep","Cases observed by Oct","Cases observed by Nov","Cases observed by Dec"))
    tmp_comb1$date1 <- factor(tmp_comb1$date1, levels=c("Sep","Oct","Nov","Dec"))
    fig2B <- tmp_comb1 %>%   
        mutate(model=paste0("k=",model)) %>%
        ggplot() +
        geom_bar(aes(x=model,y=median_para,fill=model),stat="identity",col="black") +
        geom_errorbar(aes(x=model,ymin=lower,ymax=upper),width=0.5) +
        facet_wrap(~date1,nrow=1,strip.position="bottom") +
        xlab("Month by which k additional cases observed") +
        ylab(paste0("Projected paralysis cases\n by ", max_date)) +
        theme_classic() +
        scale_fill_nejm() +
        theme(legend.position="none",strip.background = element_blank(),
              panel.grid.major=element_line(size=0.1,color="grey70"),
              strip.placement = "outside"
        ) +
        scale_y_continuous(expand=c(0,0),limits=c(0,ymax))
    return(list(fig2A,fig2B))
}
fig2_rockland <- plot_fig2(tmp_comb)
fig2_nyc <- plot_fig2(tmp_comb_nyc)


ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig1.pdf", fig1, height=8,width=7)
ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig1.png", fig1, height=8,width=7,units='in',dpi=300)

fig2_main <- (fig2_rockland[[2]]+labs(tag="A") + ggtitle("Rockland County") + theme(plot.title=element_text(size=10)))/(fig2_nyc[[2]] + labs(tag="B")+ ggtitle("New York City")+ theme(plot.title=element_text(size=10)))
fig2_alt <- (fig2_rockland[[1]]+labs(tag="A") + ggtitle("Rockland County") + theme(plot.title=element_text(size=10)))/(fig2_nyc[[1]] + labs(tag="B")+ ggtitle("New York City")+ theme(plot.title=element_text(size=10)))


ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig2_main.pdf", fig2_main, height=5,width=8)
ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig2_main.png", fig2_main, height=5,width=8,units="in",dpi=300)
ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig2_alt.pdf", fig2_alt, height=5,width=8)
ggsave(filename="~/Documents/GitHub/paralytic_polio_estimates/figures/fig2_alt.png", fig2_alt, height=5,width=8,units="in",dpi=300)


## Numbers for paper
## % showing continued spread by August 20th
trajectories %>% filter(t == "2022-08-20") %>%
    group_by(ongoing_7) %>% tally() %>%
    pivot_wider(values_from=n,names_from=ongoing_7) %>%
    mutate(prop=`TRUE`/(`TRUE` + `FALSE`)) %>%
    mutate(prop = (1-prop)*100)

trajectories %>% filter(t == "2022-08-20") %>%
    filter(`Consistent with further\nparalytic polio\n cases by 2022-10-21` == "Yes" ) %>%
    group_by(ongoing_7) %>% tally() %>%
    pivot_wider(values_from=n,names_from=ongoing_7) %>%
    mutate(prop=`TRUE`/(`TRUE` + `FALSE`)) %>%
    mutate(prop = (1-prop)*100)


## Effective reproduction number at seeding
trajectories %>% group_by(sim) %>%
    filter(t == min(t)) %>%
    ungroup() %>%
    dplyr::summarize(mean_re=mean(Rt),median_re=median(Rt),
              lower90=quantile(Rt,0.1),upper90=quantile(Rt,0.9))

## Generation interval estimate
res$infectious_period_mean %>% mean
quantile(res$infectious_period_mean,c(0.1,0.9))

## Doubling time

## 22nd June number of infections
trajectories %>% filter(t == as.Date("2022-06-22")) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## 18nd July number of infections
trajectories %>% filter(t == as.Date("2022-07-18")) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## 20th August number of infections
trajectories %>% filter(t == as.Date("2022-08-20")) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## 1st April 2023 number of infections
## 22nd June number of infections
trajectories %>% filter(t == as.Date(max_date)) %>% 
    ungroup() %>% summarize(mean_inf = mean(cumu_inc),
                            median_inf = median(cumu_inc),
                            lower90=quantile(cumu_inc,0.1),
                            upper90=quantile(cumu_inc,0.9))

## Number of paralytic polio by 1st April
trajectories %>% filter(t == as.Date(max_date)) %>% 
    ungroup() %>% summarize(mean_para = mean(cumu_para),
                            median_para = median(cumu_para),
                            lower90=quantile(cumu_para,0.1),
                            upper90=quantile(cumu_para,0.9))

## Number of paralytic polio by 20th August if no further cases by 1st Oct
traj_summary %>% filter(t == as.Date("2022-08-20")) %>%
    mutate(prop=1-prop,lower=1-lower,upper=1-upper)

## Number of paralytic polio by 1st April if 1 further case by 1st Oct
traj_summary %>% filter(t == max_date)

## Probability of not seeing any cases in NYC at all
traj_nyc %>% group_by(sim) %>% summarize(total_para = sum(para)) %>%
    mutate(any_para=total_para > 0) %>% group_by(any_para) %>%
    tally() %>% pivot_wider(names_from=any_para,values_from=n) %>%
    mutate(prop = `FALSE`/(`FALSE`+`TRUE`))

## Number of paralysis in NYC if we see 1 case on or after 1st October by 1st April 2023
tmp_comb_nyc %>% filter(t == "2022-10-01")


## Supplementary figures
priors <- read_csv("~/Documents/GitHub/paralytic_polio_estimates/pars/priors.csv")
scenarios <- c("rockland_high_coverage","rockland_low_coverage")
tmp_pars <- priors %>% filter(scenario == scenarios[1])
prior_draws <- simulate_priors(10000,incu_mean_prior_mean=14,
                               incu_mean_prior_var=3,
                               incu_var_prior_mean=15,
                               incu_var_prior_var=1,
                               gen_interval_susc_shape_par1=tmp_pars %>% 
                                   filter(`model parameter` == "susceptible_generation_interval_shape") %>% pull(par1),
                               gen_interval_susc_shape_par2=tmp_pars %>% 
                                   filter(`model parameter` == "susceptible_generation_interval_shape") %>% pull(par2),
                               gen_interval_susc_rate_par1=tmp_pars %>% 
                                   filter(`model parameter` == "susceptible_generation_interval_rate") %>% pull(par1),
                               gen_interval_susc_rate_par2=tmp_pars %>% 
                                   filter(`model parameter` == "susceptible_generation_interval_rate") %>% pull(par2),
                               gen_interval_partial_shape_par1=tmp_pars %>% 
                                   filter(`model parameter` == "partial_generation_interval_shape") %>% pull(par1),
                               gen_interval_partial_shape_par2=tmp_pars %>% 
                                   filter(`model parameter` == "partial_generation_interval_shape") %>% pull(par2),
                               gen_interval_partial_rate_par1=tmp_pars %>% 
                                   filter(`model parameter` == "partial_generation_interval_rate") %>% pull(par1),
                               gen_interval_partial_rate_par2=tmp_pars %>% 
                                   filter(`model parameter` == "partial_generation_interval_rate") %>% pull(par2),
                               prob_paralysis_mean=0.0005,
                               prob_paralysis_ps_par1 = -4,
                               prob_paralysis_ps_par2 = -1,
                               prob_paralysis_var = 2e-8,
                               R0_dist="truncnorm",
                               R0_par1=4.9,R0_par2=2,
                               rel_R0_par1=tmp_pars %>% filter(`model parameter` == "relative_infectiousness") %>% pull(par1),
                               rel_R0_par2=tmp_pars %>% filter(`model parameter` == "relative_infectiousness") %>% pull(par2),
                               prop_susceptible_par1 = tmp_pars %>% filter(`model parameter` == "prop_susceptible") %>% pull(par1),
                               prop_susceptible_par2 = tmp_pars %>% filter(`model parameter` == "prop_susceptible") %>% pull(par2),
                               prop_refractory_par1 = tmp_pars %>% filter(`model parameter` == "prop_refractory") %>% pull(par1),
                               prop_refractory_par2 = tmp_pars %>% filter(`model parameter` == "prop_refractory") %>% pull(par2)
)
colnames(prior_draws)


prior_key <- c("R0"="R[0]", 
               "Re"="R[e]",
               "rel_R0s"="φ", 
               "rel_R0"="φ", 
               "infect_mean"="g[S]^μ", 
               "infect_var"="g[S]^σ", 
               "infectious_period_var"="g[S]^σ",
               "infectious_period_mean"="g[S]^μ",
               "infectious_period_mean_ps"="f[PS]^μ",
               "infectious_period_var_ps"="f[PS]^σ",
               "infect_partial_mean"="f[PS]^μ", 
               "infect_partial_var"="f[PS]^σ", 
               "incu_mean"="γ^μ",
               "incu_var"="γ^σ",
               "prob_paralysis_s"="α[S]", 
               "prob_paralysis_ps"="α[PS]", 
               "prop_immune_groups.1"="π[R]", 
               "prop_immune_groups.2"="π[PS]", 
               "prop_immune_groups.3"="π[S]")

prior_draws <- prior_draws %>% mutate(Re = (prop_immune_groups.3 + prop_immune_groups.2) * (R0 * (prop_immune_groups.3/(prop_immune_groups.3+prop_immune_groups.2)) + R0*rel_R0s*(prop_immune_groups.3/(prop_immune_groups.3+prop_immune_groups.2))))
prior_draws <- prior_draws %>% 
    mutate(prob_paralysis_ps = (1-prob_paralysis_ps)*prob_paralysis_s)


prior_draws<- prior_draws%>%mutate("incu_mean"=incu_scale*incu_shape,"incu_var"=incu_shape*incu_scale*incu_scale)


prior_draws_long <-prior_draws %>% mutate(i=1:n()) %>% pivot_longer(-i) %>%
    filter(value<= 100)

prior_draws_long <- prior_draws_long %>% filter(name %in% names(prior_key))
prior_draws_long$Parameter = prior_key[prior_draws_long$name]



res <-res%>% mutate(Re = (prop_immune_groups.3 + prop_immune_groups.2) * (R0 * (prop_immune_groups.3/(prop_immune_groups.3+prop_immune_groups.2)) + R0*rel_R0*(prop_immune_groups.3/(prop_immune_groups.3+prop_immune_groups.2))))
res <- res %>% mutate(incu_mean=incu_shape*incu_scale,incu_var=incu_shape*incu_scale*incu_scale)

res1 <- res[,colnames(res) %in% c("sim",names(prior_key))]
res1 <- res1 %>% mutate(prob_paralysis_ps = (1-prob_paralysis_ps)*prob_paralysis_s)
res1 <- res1 %>% pivot_longer(-sim)
res1 <- res1 %>% filter(value < 100)
res1$Parameter <- prior_key[res1$name]

all_draws <- bind_rows(prior_draws_long %>% mutate(Distribution="Prior"), 
                       res1 %>% mutate(Distribution="Post-filter"))

all_draws <- all_draws %>% filter(name != "prob_paralysis_ps" | (name == "prob_paralysis_ps" & value < 1e3))

p_prior <- ggplot(all_draws) + 
    geom_density(aes(x=value,fill=Distribution),alpha=0.25) +
    geom_vline(data=data.frame(Parameter=c("R[0]","R[e]"),xintercept=1),aes(xintercept=xintercept),linetype="dashed",col="red")+
    facet_wrap(~Parameter,scales="free",labeller=label_parsed,ncol=3) +
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position=c(0.8,0.1)) +
    scale_fill_manual(values=c("Post-filter"=nejm_palette[3],"Prior"="grey40"))+
    xlab("Value") +
    ylab("Density")

ggsave(filename = "figures/prior_distributions.pdf",p_prior,height=8,width=8)
ggsave(filename = "figures/prior_distributions.png",p_prior,height=8,width=8,units="in",dpi=300)



## Get post-filter and prior draws for generation interval dist and incubation period dist
## Generation interval plots
incubation_period <- function(t, incu_scale, incu_shape, max_incu_period=50){
    extraDistr::ddgamma(t,scale=incu_scale,shape=incu_shape)/extraDistr::pdgamma(max_incu_period,scale=incu_scale,shape=incu_shape)
}
infectiousness <- function(t,infect_rate,infect_shape,max_infectious_period=50){
    extraDistr::ddgamma(t, rate=infect_rate, shape=infect_shape)/extraDistr::pdgamma(max_infectious_period, shape=infect_shape,rate=infect_rate)
}


ts <- 0:25
tmp_samp <- sample(1:nrow(prior_draws),1000)
res1 <- prior_draws[tmp_samp,]
res2 <- res[tmp_samp,]
incu_periods_prior <- matrix(0,nrow=nrow(res1),ncol=length(ts))
infectiousness_prior <- matrix(0,nrow=nrow(res1),ncol=length(ts))
infectiousness_ps_prior <- matrix(0,nrow=nrow(res1),ncol=length(ts))
incu_periods_dat <- matrix(0,nrow=nrow(res1),ncol=length(ts))
infectiousness_dat <- matrix(0,nrow=nrow(res1),ncol=length(ts))
infectiousness_ps_dat <- matrix(0,nrow=nrow(res1),ncol=length(ts))
for(x in 1:nrow(res1)){
    ## Priors
    incu_periods_prior[x,] <- incubation_period(ts, res1$incu_scale[x],res1$incu_shape[x])
    infectiousness_prior[x,] <- infectiousness(ts, infect_rate=res1$infect_rate[x],infect_shape=res1$infect_shape[x])
    infectiousness_ps_prior[x,] <- infectiousness(ts, infect_rate=res1$infect_partial_rate[x],infect_shape=res1$infect_partial_shape[x])
    
    ## Post filter
    incu_periods_dat[x,] <- incubation_period(ts, res2$incu_scale[x],res2$incu_shape[x])
    infectiousness_dat[x,] <- infectiousness(ts, infect_rate=res2$infect_rate[x],infect_shape=res2$infect_shape[x])
    infectiousness_ps_dat[x,] <- infectiousness(ts, infect_rate=res2$infect_partial_rate[x],infect_shape=res2$infect_partial_shape[x])
}
incu_periods_dat <- reshape2::melt(incu_periods_dat)
colnames(incu_periods_dat) <- c("sim","time_since_infection","probability")
incu_periods_dat$model <- "Paralysis incubation period"
infectiousness_dat <- reshape2::melt(infectiousness_dat)
colnames(infectiousness_dat) <- c("sim","time_since_infection","probability")
infectiousness_dat$model <- "Fully susceptible\ngeneration interval"
infectiousness_ps_dat <- reshape2::melt(infectiousness_ps_dat)
colnames(infectiousness_ps_dat) <- c("sim","time_since_infection","probability")
infectiousness_ps_dat$model <- "Partially susceptible\ngeneration interval"

dists <- bind_rows(incu_periods_dat,infectiousness_dat,infectiousness_ps_dat)

incu_periods_prior <- reshape2::melt(incu_periods_prior)
colnames(incu_periods_prior) <- c("sim","time_since_infection","probability")
incu_periods_prior$model <- "Paralysis incubation period"
infectiousness_prior <- reshape2::melt(infectiousness_prior)
colnames(infectiousness_prior) <- c("sim","time_since_infection","probability")
infectiousness_prior$model <- "Fully susceptible\ngeneration interval"
infectiousness_ps_prior <- reshape2::melt(infectiousness_ps_prior)
colnames(infectiousness_ps_prior) <- c("sim","time_since_infection","probability")
infectiousness_ps_prior$model <- "Partially susceptible\ngeneration interval"

dists_prior <- bind_rows(incu_periods_prior,infectiousness_prior,infectiousness_ps_prior)

dists_all<- bind_rows(dists%>%mutate(Distribution="Post-filter"),
                      dists_prior %>% mutate(Distribution="Prior"))
dists_summary <- dists_all %>% group_by(Distribution,model, time_since_infection) %>%
    summarize(mean_p=mean(probability),lower_90=quantile(probability,0.1),
              upper_90=quantile(probability,0.9))


p_prior_intervals <-ggplot(dists_summary) + 
    #geom_line(aes(x=time_since_infection,y=probability,group=sim),alpha=0.1) +
    geom_ribbon(aes(x=time_since_infection,ymin=lower_90,ymax=upper_90,
                    fill=Distribution),alpha=0.25) +
    geom_line(aes(x=time_since_infection,y=mean_p,color=Distribution)) +
    geom_point(aes(x=time_since_infection,y=mean_p,color=Distribution)) +
    facet_wrap(~model,scales="free_y",ncol=3) +
    scale_x_continuous(breaks=seq(0,25,by=5))+
    theme_classic() +
    theme(strip.background = element_blank(),
          legend.position="bottom") +
    scale_fill_manual(values=c("Post-filter"=nejm_palette[3],"Prior"="grey40"))+
    scale_color_manual(values=c("Post-filter"=nejm_palette[3],"Prior"="grey40"))+
    xlab("Days since infection") +
    ylab("Probability mass")

ggsave(filename = "figures/prior_intervals.pdf",p_prior_intervals,height=3,width=8)
ggsave(filename = "figures/prior_intervals.png",p_prior_intervals,height=3,width=8,units="in",dpi=300)

get_beta_meanvar <- function(a,b){
    c(a/(a+b), (a*b)/((a+b)^2 * (a+b+1)))
}
priors1 <- cbind(priors, t(apply(priors, 1, function(x){
    x1 <- as.numeric(x)
    if(x[3] == "beta"){
        print(signif(get_beta_meanvar(x1[4], x1[5])),3)
    } else if(x[3] == "gamma") {
        print(signif(get_gamma_meanvar(x1[4], x1[5])))
    } else {
        print(signif(c(x1[4], x1[5])))
    }
})))
colnames(priors1)[6:7] <- c("mean","var")
priors1 %>% filter(scenario == "rockland_high_coverage") %>% select(`model parameter`, mean, var)
priors1 %>% filter(scenario == "NYC") %>% select(`model parameter`, mean, var)
