# Comparative analysis of stopping criteria for multi-objective evolutionary algorithms
# Manon Perrignon

source("script/01_functions.R")

########################################
########## Industrial case #############
########################################
# Evaluation of the stopping criteria on a cheese-making process optimization problem : Sect 4.3
## Import random forest models
modY1 <- read_rds("data/industrial_case/modY1.rds")
modY2 <- read_rds("data/industrial_case/modY2.rds")
modY3 <- read_rds("data/industrial_case/modY3.rds")
modY4 <- read_rds("data/industrial_case/modY4.rds")

## Import bounds
Mini <- read_rds("data/industrial_case/bornes_min_optim.rds")

Maxi <- read_rds("data/industrial_case/bornes_max_optim.rds")

### Import fixed variables
Milk <- read_rds("data/industrial_case/Var_fixe_optim.rds")

## Constraint
constraint_industrial <- function(x) {
  contrainte <- numeric(11)
  contrainte[1] <- (sum(x[c(1:4,6,7)])) -  0.98 
  contrainte[2] <- 1.02- (sum(x[c(1:4,6,7)])) 
  contrainte[3] <- 0.035- x[7] 
  contrainte[4] <- 41.6 - x[24] 
  contrainte[5] <- x[24] - 41.4 
  contrainte[6] <- 39.4 - x[23] 
  contrainte[7] <- x[23] - 39.2 
  contrainte[8] <- (x[23]+0.2)-((x[1]*44.55)+(x[2]*0.03)+(x[3]*442.2)+(x[6]*0.9)+
                                  (x[7]*14)+(x[4]*116.7)) # bilan mat MG
  contrainte[9] <- ((x[1]*44.55)+(x[2]*0.03)+(x[3]*442.2)+(x[6]*0.9)+(x[7]*14)+
                      (x[4]*116.7))- (x[23]-0.2) # bilan mat MG
  contrainte[10] <- (x[24]+0.2)-((x[1]*34.49)+(x[2]*35.5)+(x[3]*17.2)+(x[6]*96.3)+
                                   (x[7]*31)+(x[4]*4.6)) # bilan mat MP
  contrainte[11] <- ((x[1]*34.49)+(x[2]*35.5)+(x[3]*17.2)+(x[6]*96.3)+(x[7]*31)+
                       (x[4]*4.6))-(x[24]-0.2) # bilan mat MP
  return(contrainte)
}

## Objective function
Funct_milk <- function(x){
  y <- Milk
  y[,38] <- 0
  y[,39] <- 0
  x <- data.frame(t(x))
  x <- cbind(x,y)
  names(x) <- modY1[["forest"]][["independent.variable.names"]]
  
  predY1 <<- predict(modY1,data=x)
  predY1 <- abs(60.8-predY1$predictions)
  predY2 <<- predict(modY2,data=x)
  predY2 <- predY2$predictions
  predY3 <- predict(modY3,data=x)
  predY3 <- predY3$predictions
  predY4 <- predict(modY4,data=x)
  predY4 <- predY4$predictions
  return(c(predY1,predY2,predY3,predY4))
}

## Function to obtain 100 repetitions and results
optim_rep_real <- function(probleme,var,obj,low,upp,rep){
  result_tot <- data.frame("critere"=NA, "HV"=NA, "Spread"=NA,"entropy" = NA, "MDR"=NA , "gen"=NA, "time" = NA, "time_tot" = NA ,"rep" = NA, "nb_obj" = NA)
  tot_crit <- c()
  
  for (i in 1:rep){
    optim <- mco::nsga2(probleme, var, obj, generations = c(1:5000), popsize = 100,
                        lower.bounds = Mini, upper.bounds = Maxi, constraints = constraint_industrial,cdim = 11) 
    
    
    front <- extract_pareto_front(optim)
    
    ref_point <- c(1,330,310,7)
    ideal_point <- c(0,320,300,3)
    
    # MPF
    tic()
    stop_mp_crit_hv <- MP_crit_HV(front, ref_point)
    time_mp_hv <- toc()
    time_mp_crit_hv <- (time_mp_hv$toc[["elapsed"]]-time_mp_hv$tic[["elapsed"]])/stop_mp_crit_hv
    
    time_tot_mp_crit_hv <- (time_mp_hv$toc[["elapsed"]]-time_mp_hv$tic[["elapsed"]])
    
    if (is.na(stop_mp_crit_hv)) {
      hv_mp_crit_hv <- NA
      spd_mp_crit_hv <- NA
      mdr_mp_crit_hv <- NA
      ent_mp_crit_hv <- NA
    } else {
      hv_mp_crit_hv <- dominatedHypervolume(front[[stop_mp_crit_hv]], ref_point)
      
      if (nrow(front[[stop_mp_crit_hv]]) <= 2) {
        spd_mp_crit_hv <- NA
      } else {
        spd_mp_crit_hv <- generalizedSpread(front[[stop_mp_crit_hv]], front[[stop_mp_crit_hv+1]])}
      
      mdr_mp_crit_hv <- mutual_dominance_rate(front[[stop_mp_crit_hv]], front[[stop_mp_crit_hv+1]])
      
      ent_mp_crit_hv <- entropy(front[[stop_mp_crit_hv]], front[[stop_mp_crit_hv+1]])
    }
    
    
    ## critere mgbm
    tic()
    stop_MGBM_crit <- MGBM_crit(front)
    time_MGBM <- toc()
    time_MGBM_crit <- (time_MGBM$toc[["elapsed"]]-time_MGBM$tic[["elapsed"]])/stop_MGBM_crit
    time_tot_MGBM_crit <- (time_MGBM$toc[["elapsed"]]-time_MGBM$tic[["elapsed"]])
    
    if (is.na(stop_MGBM_crit)) {
      hv_MGBM_crit <- NA
      spd_MGBM_crit <- NA
      mdr_MGBM_crit <- NA
      ent_MGBM_crit <- NA
    } else {
      hv_MGBM_crit <- dominatedHypervolume(front[[stop_MGBM_crit]], ref_point)
      
      if (nrow(front[[stop_MGBM_crit]]) <= 2) {
        spd_MGBM_crit <- NA
      } else {
        spd_MGBM_crit <- generalizedSpread(front[[stop_MGBM_crit]], front[[stop_MGBM_crit+1]])}
      
      mdr_MGBM_crit <- mutual_dominance_rate(front[[stop_MGBM_crit]], front[[stop_MGBM_crit+1]])
      
      ent_MGBM_crit <- entropy(front[[stop_MGBM_crit]], front[[stop_MGBM_crit+1]])
    }
    
    
    ## critere OCD_HV
    tic()
    stop_OCD_crit <- OCD_HV(front, ref_point)
    time_OCD <- toc()
    time_OCD_crit <- (time_OCD$toc[["elapsed"]]-time_OCD$tic[["elapsed"]])/stop_OCD_crit
    time_tot_OCD_crit <- (time_OCD$toc[["elapsed"]]-time_OCD$tic[["elapsed"]])
    
    if (is.na(stop_OCD_crit)) {
      hv_OCD_crit <- NA
      spd_OCD_crit <- NA
      mdr_OCD_crit <- NA
      ent_OCD_crit <- NA
    } else {
      hv_OCD_crit <- dominatedHypervolume(front[[stop_OCD_crit]], ref_point)
      
      if (nrow(front[[stop_OCD_crit]]) <= 2) {
        spd_OCD_crit <- NA
      } else {
        spd_OCD_crit <- generalizedSpread(front[[stop_OCD_crit]], front[[stop_OCD_crit+1]])}
      
      mdr_OCD_crit <- mutual_dominance_rate(front[[stop_OCD_crit]], front[[stop_OCD_crit+1]])
      
      ent_OCD_crit <- entropy(front[[stop_OCD_crit]], front[[stop_OCD_crit+1]])
    }
    
    ## critere LSSC
    tic()
    stop_LSSC_crit <- LSSC(front, ref_point)
    time_LSSC <- toc()
    time_LSSC_crit <- (time_LSSC$toc[["elapsed"]]-time_LSSC$tic[["elapsed"]])/stop_LSSC_crit
    time_tot_LSSC_crit <- (time_LSSC$toc[["elapsed"]]-time_LSSC$tic[["elapsed"]])
    
    if (is.na(stop_LSSC_crit)) {
      hv_LSSC_crit <- NA
      spd_LSSC_crit <- NA
      mdr_LSSC_crit <- NA
      ent_LSSC_crit <- NA
    } else {
      hv_LSSC_crit <- dominatedHypervolume(front[[stop_LSSC_crit]], ref_point)
      
      if (nrow(front[[stop_LSSC_crit]]) <= 2) {
        spd_LSSC_crit <- NA
      } else {
        spd_LSSC_crit <- generalizedSpread(front[[stop_LSSC_crit]], front[[stop_LSSC_crit+1]])}
      
      mdr_LSSC_crit <- mutual_dominance_rate(front[[stop_LSSC_crit]], front[[stop_LSSC_crit+1]])
      
      ent_LSSC_crit <- entropy(front[[stop_LSSC_crit]], front[[stop_LSSC_crit+1]])
    }
    
    ## critere ertropy
    tic()
    stop_Ent_crit <- Entropy_crit(front)
    time_Ent <- toc()
    time_Ent_crit <- (time_Ent$toc[["elapsed"]]-time_Ent$tic[["elapsed"]])/stop_Ent_crit
    time_tot_Ent_crit <- (time_Ent$toc[["elapsed"]]-time_Ent$tic[["elapsed"]])
    
    if (is.na(stop_Ent_crit)) {
      hv_Ent_crit <- NA
      spd_Ent_crit <- NA
      mdr_Ent_crit <- NA
      ent_Ent_crit <- NA
    } else {
      hv_Ent_crit <- dominatedHypervolume(front[[stop_Ent_crit]], ref_point)
      
      if (nrow(front[[stop_Ent_crit]]) <= 2) {
        spd_Ent_crit <- NA
      } else {
        spd_Ent_crit <- generalizedSpread(front[[stop_Ent_crit]], front[[stop_Ent_crit+1]])}
      
      mdr_Ent_crit <- mutual_dominance_rate(front[[stop_Ent_crit]], front[[stop_Ent_crit+1]])
      
      ent_Ent_crit <- entropy(front[[stop_Ent_crit]], front[[stop_Ent_crit+1]])
    }
    
    ## Recuperation des resultats
    tot_crit$critere <- c("MPF","MGBM", "OCD_HV", "LSSC", "Entropy")
    tot_crit$HV <- c(hv_mp_crit_hv,hv_MGBM_crit,hv_OCD_crit,hv_LSSC_crit,hv_Ent_crit)
    tot_crit$Spread <- c(spd_mp_crit_hv,spd_MGBM_crit,spd_OCD_crit,spd_LSSC_crit,spd_Ent_crit)
    tot_crit$entropy <- c(ent_mp_crit_hv,ent_MGBM_crit,ent_OCD_crit,ent_LSSC_crit,ent_Ent_crit)
    tot_crit$R2 <- c(r2_mp_crit_hv,r2_MGBM_crit,r2_OCD_crit,r2_LSSC_crit,r2_Ent_crit)
    tot_crit$MDR<- c(mdr_mp_crit_hv,mdr_MGBM_crit,mdr_OCD_crit,mdr_LSSC_crit,mdr_Ent_crit)
    tot_crit$gen <- c(stop_mp_crit_hv,stop_MGBM_crit,stop_OCD_crit,stop_LSSC_crit,stop_Ent_crit)
    tot_crit$time <- c(time_mp_crit_hv,time_MGBM_crit,time_OCD_crit,time_LSSC_crit,time_Ent_crit)
    tot_crit$time_tot <- c(time_tot_mp_crit_hv,time_tot_MGBM_crit,time_tot_OCD_crit,time_tot_LSSC_crit,time_tot_Ent_crit)
    tot_crit$rep <- i
    tot_crit$nb_obj <- obj
    tot_crit <- as.data.frame(tot_crit)
    
    result_tot <- rbind(result_tot,tot_crit)
    print(i)
  }
  return(result_tot[-1,])
}

set.seed(123)
results_real_case <- optim_rep_real(Funct_milk,36,4,Mini,Maxi,100)

write_rds(results_real_case,"data/industrial_case/results_industrial_case.rds")