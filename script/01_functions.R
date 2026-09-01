# Comparative analysis of stopping criteria for multi-objective evolutionary algorithms
# Manon Perrignon

########################## OPTIMIZATION ##################################

# Load packages 

install.packages("mco")
install.packages("readr")
install.packages("tidyverse")
install.packages("smoof")
install.packages("tictoc")
install.packages("rPref")
install.packages("ranger")

library(mco)
library(readr)
library(tidyverse)
library(smoof)
library(tictoc)
library(rPref)
library(ranger)

## Creating the necessary functions
### Evaluation metrics 
#### Entropy
MultidimensionalHistogram <- function(P, Q, nb) {
  # get the cell ID
  GetCellId <- function(solution, nb) {
    # generate a cell ID
    return(floor(solution * nb))  # assumes that the objctives are normalised between 0 and 1
  }
  
  C <- list()
  Cq <- list()
  Pc <- c()
  Qc <- c()
  Qcq <- c()
  
  # scan population P
  for (i in 1:nrow(P)) {
    s <- P[i, ]
    c <- paste(GetCellId(s, nb), collapse = "-")
    
    if (c %in% C) {
      k <- match(c, C)
      Pc[k] <- Pc[k] + 1
    } else {
      C <- c(C, c)
      Pc <- c(Pc, 1)
      Qc <- c(Qc, 0)
    }
  }
  
  # scan population Q
  for (i in 1:nrow(Q)) {
    s <- Q[i, ]
    c <- paste(GetCellId(s, nb), collapse = "-")
    
    if (c %in% C) {
      k <- match(c, C)
      Qc[k] <- Qc[k] + 1
    } else if (c %in% Cq) {
      k <- match(c, Cq)
      Qcq[k] <- Qcq[k] + 1
    } else {
      Cq <- c(Cq, c)
      Qcq <- c(Qcq, 1)
    }
  }
  
  return(list(C = C, Cq = Cq, Pc = Pc, Qc = Qc, Qcq = Qcq))
}

entropy <- function(P, Q, nb = 10) {
  # Multidimensional histogram
  histo <- MultidimensionalHistogram(P, Q, nb)
  C <- histo$C
  Cq <- histo$Cq
  Pc <- histo$Pc
  Qc <- histo$Qc
  Qcq <- histo$Qcq
  
  Dt <- 0
  
  # Cellules communes
  for (i in seq_along(C)) {
    p <- Pc[i] / nrow(P)
    q <- Qc[i] / nrow(Q)
    
    if (q > 0 && p > 0) {
      Dt <- Dt - ((q / 2) * log2(q / p) + (p / 2) * log2(p / q))
    } else if (p > 0) {
      Dt <- Dt - p * log2(p)
    }
  }
  
  # Cells specific to Q
  for (i in seq_along(Cq)) {
    q <- Qcq[i] / nrow(Q)
    Dt <- Dt - q * log2(q)
  }
  
  return(Dt)
}

#### Mutual dominance rate
dominates <- function(a, b) {
  all(a <= b) && any(a < b)
}

mutual_dominance_rate <- function(front_n, front_n_1) {
  n <- nrow(front_n)
  m <- nrow(front_n_1)
  
  
  count_dominance_n <- sum(sapply(1:n, function(i) {
    any(sapply(1:m, function(j) dominates(front_n[i, ], front_n_1[j, ])))
  }))
  
  count_dominance_n_1 <- sum(sapply(1:m, function(j) {
    any(sapply(1:n, function(i) dominates(front_n_1[j, ], front_n[i, ])))
  }))
  
  MDR <-  count_dominance_n / n - count_dominance_n_1 / m 
  return(MDR)
}

#### For Hypervolume and Spread, functions are included in the package: mco

### chi2 test
chi2_test <- function(indic, varLimit) {
  N <- length(indic) - 1
  observed_var <- var(indic)
  chi_stat <- observed_var * N / varLimit
  p_value <- pchisq(chi_stat,N,lower.tail = T)
  return(p_value)
}

### Extraction of the Pareto front
extract_pareto_front <- function(MOO) {
  lapply(seq_along(MOO), function(i) {
    MOO[[i]]$value[MOO[[i]]$pareto.optimal, , drop = FALSE]
  })
}

### Sliding window
sliding_window <- function(x, window_size) {
  windows <- list()
  
  for (i in 1:(length(x) - window_size + 1)) {
    windows[[i]] <- x[i:(i + window_size - 1)]
  }
  
  return(windows)
}

## Implementation of criteria 

### MGBM
MGBM_crit <- function(fronts, R = 0.1, Q = 0.0001, seuil = 0.05) {
  num_generations <- length(fronts)
  if (num_generations < 2) stop("Pas assez de générations")
  
  # Initialisation
  A <- 1
  P <- 0.1
  mdr_values <- numeric(num_generations - 1)
  filtered_values <- numeric(num_generations - 1)
  
  seuil_atteint_index <- NA
  
  gen_init <- 2
  while (gen_init <= num_generations && (is.null(fronts[[gen_init]]) || is.null(fronts[[gen_init - 1]]) ||
                                         nrow(fronts[[gen_init]]) == 0 || nrow(fronts[[gen_init - 1]]) == 0)) {
    gen_init <- gen_init + 1
  }
  if (gen_init > num_generations) return(NA)
  
  mdr_values[gen_init - 1] <- mutual_dominance_rate(fronts[[gen_init]], fronts[[gen_init - 1]])
  x_hat <- mdr_values[gen_init - 1]
  filtered_values[gen_init - 1] <- x_hat
  
  if (x_hat <= seuil) {
    seuil_atteint_index <- gen_init
    return(as.numeric(seuil_atteint_index))
  }
  
  
  
  for (gen in (gen_init + 1):num_generations) {
    mdr <- mutual_dominance_rate(fronts[[gen]], fronts[[gen - 1]])
    mdr_values[gen - 1] <- mdr
    
    # kalman filter
    x_hat_pred <- A * x_hat
    P_pred <- A * P * A + Q
    
    # updated
    K <- P_pred / (P_pred + R)
    x_hat <- x_hat_pred + K * (mdr - x_hat_pred)
    P <- (1 - K) * P_pred
    
    filtered_values[gen - 1] <- x_hat
    if (is.na(seuil_atteint_index) && x_hat <= seuil) {
      seuil_atteint_index <- gen
      break
    }
  }
  
  return(as.numeric(seuil_atteint_index))
}

## OCD_HV 
OCD_HV <- function(front, ref_point, varLimit = 1e-3, window = 10, threshold = 0.05) {
  
  n_gen   <- length(front)
  m       <- ncol(front[[1]])    
  start   <- window                    
  
  # check the front outside the for loop
  valid <- vapply(front, function(f) !is.null(f) && nrow(f) > 0, logical(1))
  
  # pre-calculation of HV
  HVs_all <- vapply(seq_len(n_gen), function(i) {
    if (valid[i]) dominatedHypervolume(front[[i]], ref_point) else NA
  }, numeric(1))
  
  for (gen in start:n_gen) {
    
    idx <- (gen - window + 1):gen
    
    if (!all(valid[idx])) next
    
    # Direct extraction from the pre-calculated HV window
    if (chi2_test(HVs_all[idx], varLimit) <= threshold) {
      return(as.numeric(gen))
    }
  }
  
  return(as.numeric(n_gen))
}

## LSSC 
LSSC <- function(front, ref_point, window_size = 10, slope_min = 0.002) {
  n <- length(front)
  if (n < window_size) return(FALSE)
  
  # pre-calculation (removed from the for loop)
  time      <- seq_len(window_size)
  mean_res  <- 1 - (2 / window_size)
  var_res   <- (2 / window_size) - (4 / window_size^2)
  thres     <- mean_res + 3 * sqrt(var_res)
  
  # pre-calcultation
  time_mean <- mean(time)
  time_var  <- sum((time - time_mean)^2)
  
  # pre-calculation of hv
  hv_values <- vapply(front, function(f) {
    if (is.null(f) || nrow(f) == 0) NA
    else dominatedHypervolume(f, ref_point)
  }, numeric(1))
  
  # Index of HVs for reconstructing the window
  valid_idx <- which(!is.na(hv_values))
  hv_valid  <- hv_values[valid_idx]
  nv        <- length(hv_valid)
  
  if (nv < window_size) return(as.numeric(n))
  
  for (k in window_size:nv) {
    window <- hv_valid[(k - window_size + 1):k]
    
    # manual linear regression (faster than lm)
    w_mean <- mean(window)
    b      <- sum((time - time_mean) * (window - w_mean)) / time_var
    
    if (abs(b) < slope_min) {          
      a        <- w_mean - b * time_mean
      residuals <- window - (a + b * time)
      res_norm  <- sum(residuals^2) / window_size
      
      if (res_norm < thres) {
        return(as.numeric(valid_idx[k]))
      }
    }
  }
  
  return(as.numeric(n))
}

## Entropy 
Entropy_crit <- function(fronts, window=10, np=2, nb=10) {
  # initialisation
  nGen <- length(fronts)
  history_Mt <- c()
  history_St <- c()
  Dt_list <- c()
  c1 <- FALSE
  c2 <- FALSE
  
  for (t in 1:(nGen-1)) {
    P <- fronts[[t]]
    Q <- fronts[[t+1]]
    
    if (is.null(P) || is.null(Q) || nrow(P) == 0 || nrow(Q) == 0) {
      next
    }
    
    # multidimensional histogram 
    histo <- MultidimensionalHistogram(P, Q, nb)
    C <- histo$C
    Cq <- histo$Cq
    Pc <- histo$Pc
    Qc <- histo$Qc
    Qcq <- histo$Qcq
    
    Dt <- 0
    
    # search for common cells
    for (i in seq_along(C)) {
      p <- Pc[i] / nrow(P)
      q <- Qc[i] / nrow(Q)
      # application of the Kullback–Leibler divergence formula
      if (q > 0) {
        Dt <- Dt - ((q / 2) * log2(q / p) + (p / 2) * log2(p / q))
      } else {
        Dt <- Dt - p * log2(p)
      }
    }
    
    # cells specific to Q
    for (i in seq_along(Cq)) {
      q <- Qcq[i] / nrow(Q)
      Dt <- Dt - q * log2(q)
    }
    
    Dt_list <- c(Dt_list, Dt)
    
    
    history_Mt[t] <- round(mean(Dt_list), np)
    history_St[t] <- round(sd(Dt_list), np)
    
    if (t > window) {
      if (
        all(!is.na(history_Mt[(t - (window - 1)):t])) &&
        all(!is.na(history_St[(t - (window - 1)):t]))
      ) {
        if (all(history_Mt[(t - (window-1)):t] == history_Mt[t])) c1 <- TRUE
        if (all(history_St[(t - (window-1)):t] == history_St[t])) c2 <- TRUE
      }
      
      if (c1 && c2) {
        return(as.numeric(t))
      }
    }
    
  }
  return(as.numeric(t))
}

## MPF
MP_crit_HV <- function(fronts, ref_point,window_size = 10, threshold = 1e-3, dec=2) {
  
  hv_values <- c()
  ent_values <- c()
  previous_crit_hv <- NA
  previous_crit_ent <- NA
  
  for (i in 2:length(fronts)) {
    front_i   <- fronts[[i]]
    front_prev <- fronts[[i - 1]]
    
    if (is.null(front_i) || is.null(front_prev) || nrow(front_i) == 0 || nrow(front_prev) == 0) {
      next
    }
    
    # Calculation of HV and Entropy
    hv_i  <- dominatedHypervolume(front_i, ref_point)
    ent_i <- entropy(front_prev, front_i)
    
    hv_values <- c(hv_values, hv_i)
    ent_values <- c(ent_values, ent_i)
    
    if (i >= window_size) {
      window_hv <- tail(hv_values, window_size)
      window_ent <- tail(ent_values, window_size)
      
      crit_stop_hv <- abs(mean(window_hv))
      
      crit_stop_ent <- abs(mean(window_ent))
      
      if (!is.na(previous_crit_hv) && !is.na(previous_crit_ent)) {
        delta_hv <- abs((crit_stop_hv - previous_crit_hv)/previous_crit_hv)*100
        
        
        if (delta_hv < threshold && round(previous_crit_ent,dec)==round(crit_stop_ent,dec)) {
          return(as.numeric(i))
        }
        
      }
      # Updates the criteria for the next comparison
      previous_crit_hv <- crit_stop_hv
      previous_crit_ent <- crit_stop_ent
      
    }
    
    
  }
  return(as.numeric(i))
}


## Function to obtain results after 100 repetitions
optim_rep <- function(probleme,var,obj,low,upp,rep,ref){
  result_tot <- data.frame("critere"=NA, "HV"=NA, "Spread"=NA,"entropy" = NA, "MDR"=NA , "gen"=NA, "time" = NA, "time_tot" = NA ,"rep" = NA, "nb_obj" = NA)
  tot_crit <- c()
  
  for (i in 1:rep){
    optim <- mco::nsga2(probleme, var, obj, generations = c(1:5000), popsize = 100,
                        lower.bounds = rep(low, var), upper.bounds = upp) 
    front <- extract_pareto_front(optim)
    
    ref_point <- rep(ref,obj)
    
    
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
      
      spd_mp_crit_hv <- generalizedSpread(front[[stop_mp_crit_hv-1]], front[[stop_mp_crit_hv]])
      
      mdr_mp_crit_hv <- mutual_dominance_rate(front[[stop_mp_crit_hv-1]], front[[stop_mp_crit_hv]])
      
      ent_mp_crit_hv <- entropy(front[[stop_mp_crit_hv-1]], front[[stop_mp_crit_hv]])
    }
    
    
    ## MGBM
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
      
      spd_MGBM_crit <- generalizedSpread(front[[stop_MGBM_crit-1]], front[[stop_MGBM_crit]])
      
      mdr_MGBM_crit <- mutual_dominance_rate(front[[stop_MGBM_crit-1]], front[[stop_MGBM_crit]])
      
      ent_MGBM_crit <- entropy(front[[stop_MGBM_crit-1]], front[[stop_MGBM_crit]])
    }
    
    
    ## OCD_HV
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
      
      spd_OCD_crit <- generalizedSpread(front[[stop_OCD_crit-1]], front[[stop_OCD_crit]])
      
      mdr_OCD_crit <- mutual_dominance_rate(front[[stop_OCD_crit-1]], front[[stop_OCD_crit]])
      
      ent_OCD_crit <- entropy(front[[stop_OCD_crit-1]], front[[stop_OCD_crit]])
    }
    
    ## LSSC
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
      
      spd_LSSC_crit <- generalizedSpread(front[[stop_LSSC_crit-1]], front[[stop_LSSC_crit]])
      
      mdr_LSSC_crit <- mutual_dominance_rate(front[[stop_LSSC_crit-1]], front[[stop_LSSC_crit]])
      
      ent_LSSC_crit <- entropy(front[[stop_LSSC_crit-1]], front[[stop_LSSC_crit]])
    }
    
    ## Entropy
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
      
      spd_Ent_crit <- generalizedSpread(front[[stop_Ent_crit-1]], front[[stop_Ent_crit]])
      
      mdr_Ent_crit <- mutual_dominance_rate(front[[stop_Ent_crit-1]], front[[stop_Ent_crit]])
      
      ent_Ent_crit <- entropy(front[[stop_Ent_crit-1]], front[[stop_Ent_crit]])
    }
    
    ## Results
    tot_crit$critere <- c("MPF","MGBM", "OCD_HV", "LSSC", "Entropy")
    tot_crit$HV <- c(hv_mp_crit_hv,hv_MGBM_crit,hv_OCD_crit,hv_LSSC_crit,hv_Ent_crit)
    tot_crit$Spread <- c(spd_mp_crit_hv,spd_MGBM_crit,spd_OCD_crit,spd_LSSC_crit,spd_Ent_crit)
    tot_crit$entropy <- c(ent_mp_crit_hv,ent_MGBM_crit,ent_OCD_crit,ent_LSSC_crit,ent_Ent_crit)
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
