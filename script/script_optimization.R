# Comparative analysis of stopping criteria for multi-objective evolutionary algorithms
# Manon Perrignon

########################## OPTIMIZATION ##################################

# Load packages 
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
  # obtenir l'identifiant de cellule
  GetCellId <- function(solution, nb) {
    # générer un identifiant de cellule
    return(floor(solution * nb))  # suppose que les objectifs sont normalisés entre 0 et 1
  }
  
  C <- list()
  Cq <- list()
  Pc <- c()
  Qc <- c()
  Qcq <- c()
  
  # on parcourt la population p
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
  
  # parcourir la population Q
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
  # Histogramme multidimensionnel
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
  
  # Cellules propres à Q
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
    
    # kalman 
    x_hat_pred <- A * x_hat
    P_pred <- A * P * A + Q
    
    # maj
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
  m <- ncol(front[[1]]) #nb obj
  
  for (gen in 1:length(front)) {
    
    if (gen >= window){
      window_front <- front[(gen - window + 1):gen]
      
      if (any(sapply(window_front, function(f) is.null(f) || nrow(f) == 0))) {
        next
      }
      
      # calcul hv
      HVs <- sapply(window_front, function(window_front) dominatedHypervolume(window_front, ref_point))
      # test chi²
      p_chi2 <- chi2_test(HVs, varLimit)
      if (p_chi2 <= threshold) {
        return(as.numeric(gen))
      }
    }
    
  }
  return(as.numeric(gen))
}

## LSSC 
LSSC <- function(front, ref_point, window_size = 10, slope_min = 0.002) {
  n <- length(front)
  if (n < window_size) return(FALSE)
  
  hv_values <- c()
  
  for (i in 1:length(front)) {
    front_i <- front[[i]]
    if (is.null(front_i) || nrow(front_i) == 0) {
      next
    }
    hv_i <- dominatedHypervolume(front_i, ref_point)
    
    hv_values <- c(hv_values, hv_i)
    
    if (length(hv_values) >= window_size) {
      window <- tail(hv_values, window_size)
      time <- 1:window_size
      
      model <- lm(window~ time)
      a <- coef(model)[1]
      b <- coef(model)[2]
      
      residuals <- window - (a + b * time)
      res_norm <- sum(residuals^2) / window_size
      
      mean_res <- 1 - (2/window_size)
      var_res <- (2 / (window_size)) - (4/window_size^2)
      thres <- mean_res + 3 * sqrt(var_res)
      
      if (abs(b) < slope_min && res_norm < thres) {
        return(as.numeric(i))
      }
    }
  }
  return(as.numeric(i))
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
    
    # histogramme multidimensionnel (voir article)
    histo <- MultidimensionalHistogram(P, Q, nb)
    C <- histo$C
    Cq <- histo$Cq
    Pc <- histo$Pc
    Qc <- histo$Qc
    Qcq <- histo$Qcq
    
    Dt <- 0
    
    # on cherche les cellules communes
    for (i in seq_along(C)) {
      p <- Pc[i] / nrow(P)
      q <- Qc[i] / nrow(Q)
      # application de la formule Kullback–Leibler divergence
      if (q > 0) {
        Dt <- Dt - ((q / 2) * log2(q / p) + (p / 2) * log2(p / q))
      } else {
        Dt <- Dt - p * log2(p)
      }
    }
    
    # cellules propres à Q
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
    
    # Calcul R2 et Entropy
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
      # Met à jour les critères pour la prochaine comparaison
      previous_crit_hv <- crit_stop_hv
      previous_crit_ent <- crit_stop_ent
      
    }
    
    
  }
  return(as.numeric(i))
}


## Function to obtain results after 100 repetitions
optim_rep <- function(probleme,var,obj,low,upp,rep){
  result_tot <- data.frame("critere"=NA, "HV"=NA, "Spread"=NA,"entropy" = NA, "MDR"=NA , "gen"=NA, "time" = NA, "time_tot" = NA ,"rep" = NA, "nb_obj" = NA)
  tot_crit <- c()
  
  for (i in 1:rep){
    optim <- mco::nsga2(probleme, var, obj, generations = c(1:5000), popsize = 100,
                        lower.bounds = rep(low, var), upper.bounds = rep(upp, var)) 
    front <- extract_pareto_front(optim)
    
    ref_point <- rep(10,obj)
    ideal_point <- rep(0,obj)
    
    
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

########################################
####### Simulated Pareto Front #########
########################################
# Comparison with the simulated Pareto Front : Sect 4.1
## Creating the grid
set.seed(555)
grid_size <- 50 
x1_seq <- seq(-10, 10, length.out = grid_size)
x2_seq <- seq(-10, 10, length.out = grid_size)

X1_grid <- rep(x1_seq, each = grid_size)
X2_grid <- rep(x2_seq, times = grid_size)

## Gaussian function
gaussian_function <- function(X1, X2, centers, sigmas, heights) {
  result <- rep(0, length(X1))
  for (i in seq_along(centers)) {
    center <- centers[[i]]
    sigma <- sigmas[i]
    height <- heights[i]
    distance_squared <- (X1 - center[1])^2 + (X2 - center[2])^2
    result <- result + height * exp(-distance_squared / (2 * sigma^2)) 
  }
  return(result)
}

## Initialization of Gaussian parameters
### Y1
centers_Y1 <- list(c(-7, -5), c(1, 4))  
sigmas_Y1 <- c(2, 1)
heights_Y1 <- c(0.7, 0.5)

### Y2
centers_Y2 <- list(c(-5,-3), c(6, -4))  
sigmas_Y2 <- c(2, 1.5)  
heights_Y2 <- c(0.6, 0.4)

## Simulation of Y1 and Y2
set.seed(456)
Y1_grid <- gaussian_function(X1_grid, X2_grid, centers_Y1, sigmas_Y1,heights_Y1)
set.seed(456)
Y2_grid <- gaussian_function(X1_grid, X2_grid, centers_Y2, sigmas_Y2,heights_Y2)

## Data recovery
Y1_matrix <- matrix(Y1_grid, nrow = grid_size, ncol = grid_size)
Y2_matrix <- matrix(Y2_grid, nrow = grid_size, ncol = grid_size)

dataXY <- data.frame(
  X1 = X1_grid,
  X2 = X2_grid,
  Y1 = Y1_grid,
  Y2 = Y2_grid
)

#write_rds(dataXY,"data/simulated_front/data_simulated_front.rds")

## Calculating the pareto front using the notion of point dominance
### Indicate that we want to maximise Y1 and Y2
p_test <- high(dataXY$Y1) * high(dataXY$Y2)

pareto_test <- psel(dataXY, p_test)

## Creation of the objective function
obj_function <- function(x) {
  
  val_opt <- data.frame(t(x))
  
  predictions <- matrix(NA, nrow = nrow(val_opt), ncol = 2)
  
  predY1 <- rep(0, length(val_opt$X1))
  center_Y1 <- list(c(-7, -5), c(1, 4))
  sigmas <-  c(2, 1)
  heights <- c(0.7, 0.5)
  
  for (i in seq_along(center_Y1)) {
    center <- center_Y1[[i]]
    sigma <- sigmas[i]
    height <- heights[i]
    distance_squared <- (val_opt$X1 - center[1])^2 + (val_opt$X2 - center[2])^2
    predY1 <- predY1 + height * exp(-distance_squared / (2 * sigma^2)) # gaussienne avec hauteur
  }
  
  predY2 <- rep(0, length(val_opt$X1))
  center_Y2 <- list(c(-5,-3), c(6, -4)) 
  sigmas <-  c(2, 1.5)
  heights <- c(0.6, 0.4)
  
  for (i in seq_along(center_Y2)) {
    center <- center_Y2[[i]]
    sigma <- sigmas[i]
    height <- heights[i]
    distance_squared <- (val_opt$X1 - center[1])^2 + (val_opt$X2 - center[2])^2
    predY2 <- predY2 + height * exp(-distance_squared / (2 * sigma^2)) # gaussienne avec hauteur
  }
  
  predictions <- as.matrix(data.frame(Y1 = -predY1,Y2=-predY2))
  return(predictions = predictions)
}

## Results
results_simu <- optim_rep(obj_function,2,2,-10,10,100)

results_simu$critere <- fct_relevel(results_simu$critere, c("OCD_HV", "LSSC", "MGBM","Entropy", "MPF"))

#write_rds(results_simu,"data/simulated_front/results_simulated_front.rds")

table_simu <- results_simu %>%
  group_by(critere) %>%
  summarise(
    "Mean HV" = mean(HV, na.rm = TRUE),
    "Std HV" = sd(HV, na.rm = TRUE),
    "Mean Spread" = mean(Spread, na.rm = TRUE),
    "Std Spread" = sd(Spread, na.rm = TRUE),
    "Mean Gen" = mean(gen, na.rm = TRUE),
    "Std Gen" = sd(gen, na.rm = TRUE),
    "Mean Time" = mean(time_tot, na.rm = TRUE),
    "Std Time" = sd(time_tot, na.rm = TRUE),
    "Mean Criterion Time" = mean(time, na.rm = TRUE),
    "Std Criterion Time" = sd(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c( `Mean Gen`, `Std Gen`), ~ round(., 0))) %>%
  mutate(across(c(`Mean HV`, `Std HV`,`Mean Spread`, `Std Spread`,`Mean Time`, `Std Time`,`Mean Criterion Time`, `Std Criterion Time`), ~ round(., 3))) 

#write_rds(table_simu, "tables/table_simu.rds")

########################################
########## Benchmark problem ###########
########################################
### Results for two-objective problems : Sect 4.2.1
set.seed(2503)
results_zdt1 <- optim_rep(zdt1,30,2,0,1,100)

set.seed(2503)
wfg2_2Y <- makeWFG2Function(2,6,14)
results_wfg2_2Y <- optim_rep(wfg2_2Y,20,2,0,1,100)

set.seed(2503)
wfg3_2Y <- makeWFG3Function(2,6,14)
results_wfg3_2Y <- optim_rep(wfg3_2Y,20,2,0,1,100)

set.seed(2503)
wfg4_2Y <- makeWFG4Function(2,6,14)
results_wfg4_2Y <- optim_rep(wfg4_2Y,20,2,0,1,100)

### Results for four-objective problems : Sect 4.2.2
set.seed(123)
wfg2_4Y <- makeWFG2Function(4,6,14)
results_wfg2_4Y <- optim_rep(wfg2_4Y,20,4,0,1,100)

set.seed(2503)
wfg3_4Y <- makeWFG3Function(4,6,14)
results_wfg3_4Y <- optim_rep(wfg3_4Y,20,4,0,1,100)

set.seed(2503)
wfg4_4Y <- makeWFG4Function(4,6,14)
results_wfg4_4Y <- optim_rep(wfg4_4Y,20,4,0,1,100)

results_tot <- bind_rows(results_zdt1,results_wfg2_2Y,results_wfg3_2Y,
                         results_wfg4_2Y,results_wfg2_4Y,results_wfg3_4Y,
                         results_wfg4_4Y,.id = "Problem")

results_tot$Problem <- results_tot$Problem %>%  fct_recode("ZDT1" = "1",  "WFG2" = "2",
                                                           "WFG3" = "3", "WFG4" = "4","WFG2" = "5", "WFG3" = "6",
                                                           "WFG4" = "7")

results_tot$critere <- fct_relevel(results_tot$critere, c("OCD_HV", "LSSC", "MGBM","Entropy","MPF"))

#write_rds(results_tot,
          "data/benchmark/results_benchmark_problems.rds")

### Table for 2 objective
table_2Y <- results_tot %>%
  filter(nb_obj == 2) %>%
  group_by(Problem, critere) %>%
  summarise(
    "Mean HV" = mean(HV, na.rm = TRUE),
    "Std HV" = sd(HV, na.rm = TRUE),
    "Mean Spread" = mean(Spread, na.rm = TRUE),
    "Std Spread" = sd(Spread, na.rm = TRUE),
    "Mean Gen" = mean(gen, na.rm = TRUE),
    "Std Gen" = sd(gen, na.rm = TRUE),
    "Mean Time" = mean(time_tot, na.rm = TRUE),
    "Std Time" = sd(time_tot, na.rm = TRUE),
    "Mean Criterion Time" = mean(time, na.rm = TRUE),
    "Std Criterion Time" = sd(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c( `Mean Gen`, `Std Gen`), ~ round(., 0))) %>%
  mutate(across(c(`Mean HV`, `Std HV`,`Mean Spread`, `Std Spread`,`Mean Time`, `Std Time`,`Mean Criterion Time`, `Std Criterion Time`), ~ round(., 3))) 

#write_rds(table_2Y, "tables/table_2Y.rds")

### Table for 4 objective
table_4Y <- results_tot %>%
  filter(nb_obj == 4) %>%
  group_by(Problem, critere) %>%
  summarise(
    "Mean HV" = mean(HV, na.rm = TRUE),
    "Std HV" = sd(HV, na.rm = TRUE),
    "Mean Spread" = mean(Spread, na.rm = TRUE),
    "Std Spread" = sd(Spread, na.rm = TRUE),
    "Mean Gen" = mean(gen, na.rm = TRUE),
    "Std Gen" = sd(gen, na.rm = TRUE),
    "Mean Time" = mean(time_tot, na.rm = TRUE),
    "Std Time" = sd(time_tot, na.rm = TRUE),
    "Mean Criterion Time" = mean(time, na.rm = TRUE),
    "Std Criterion Time" = sd(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c( `Mean Gen`, `Std Gen`), ~ round(., 0))) %>%
  mutate(across(c(`Mean HV`, `Std HV`,`Mean Spread`, `Std Spread`,`Mean Time`, `Std Time`), ~ round(., 2))) %>% 
  mutate(across(c(`Mean Criterion Time`, `Std Criterion Time`), ~ round(., 3))) 

#write_rds(table_4Y, "tables/table_4Y.rds")


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

#write_rds(results_real_case,"data/industrial_case/results_industrial_case.rds")

## Table for industrial case
table_real_case <- results_real_case %>%
  group_by(critere) %>%
  summarise(
    "Mean HV" = mean(HV, na.rm = TRUE),
    "Std HV" = sd(HV, na.rm = TRUE),
    "Mean Spread" = mean(Spread, na.rm = TRUE),
    "Std Spread" = sd(Spread, na.rm = TRUE),
    "Mean Gen" = mean(gen, na.rm = TRUE),
    "Std Gen" = sd(gen, na.rm = TRUE),
    "Mean Time" = mean(time_tot, na.rm = TRUE),
    "Std Time" = sd(time_tot, na.rm = TRUE),
    "Mean Criterion Time" = mean(time, na.rm = TRUE),
    "Std Criterion Time" = sd(time, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(c(`Mean Gen`, `Std Gen`), ~ round(., 0))) %>%
  mutate(across(c(`Mean HV`, `Std HV`, `Mean Spread`, `Std Spread`, 
                  `Mean Time`, `Std Time`, `Mean Criterion Time`, `Std Criterion Time`), 
                ~ round(., 3)))

#write_rds(table_real_case,"tables/table_industrial_case.rds")

########################## VISUALIZATION ##################################

# Load packages 
library(tidyverse)
library(ggplot2)
library(readr)
library(rPref)
library(patchwork)
library(ggh4x)
library(mco)
library(ggpubr)

########################################
####### Simulated Pareto Front #########
########################################
# Comparison with the simulated Pareto Front : Sect 4.1
# Data recovery
dataXY <- read_rds("data/simulated_front/data_simulated_front.rds")

results_simu <- read_rds("data/simulated_front/results_simulated_front.rds")

## Calculating the pareto front using the notion of point dominance
### Indicate that we want to maximise Y1 and Y2
p_test <- high(dataXY$Y1) * high(dataXY$Y2)

pareto_test <- psel(dataXY, p_test)

# Figure 1
plot_Y1 <- ggplot(dataXY, aes(x = X1, y = X2, z = Y1)) +
  geom_contour(colour = "dodgerblue3") +
  labs(x = "X1", y = "X2") +
  xlim(-10,10) + ylim(-10,10) +
  labs(title = "Y1 based on X1 and X2") +
  theme_bw() + 
  theme(
    plot.title = element_text(face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 


plot_Y2 <- ggplot(dataXY, aes(x = X1, y = X2, z = Y2)) +
  geom_contour(colour = "dodgerblue3") +
  labs(x = "X1", y = "X2") +
  xlim(-10,10) + ylim(-10,10) +
  theme_bw() + labs(title = "Y2 based on X1 and X2") +
  theme(
    plot.title = element_text(face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 
ParetFront_theo <- ggplot(pareto_test) + geom_line(aes(x=Y1, y=Y2)) +
  theme_bw() + labs(title = "Approximate Pareto Front")  +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 


(plot_Y1 | plot_Y2) / ParetFront_theo 

# Figure 2
results_simu$time <- log(results_simu$time)
results_simu$time_tot <- log(results_simu$time_tot)
results_simu_long <- results_simu %>% select(HV,critere,gen,time,time_tot,rep,Spread) %>% pivot_longer(!c(rep,critere), names_to = "Variables", values_to = "Valeurs")

supp.labs <- c("Stopping Generation", "Hypervolume", "Spread" ,"Criterion Time", "Overall Time")
names(supp.labs) <- c("gen", "HV", "Spread", "time","time_tot")

results_simu_long %>% 
  ggplot(aes(x = critere, y = Valeurs,fill=Variables)) + 
  geom_boxplot( alpha = 0.3) + 
  facet_grid2(~Variables,scale="free",independent="all", labeller = labeller(Variables = supp.labs)) +  
  labs(
    title = "",
    x = "Stopping criterion",
    y = "Values"
  ) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Dark2") +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

# Figure 3
## Recovery of Pareto fronts for each stop proposed by the criteria for the first repetition
set.seed(2503)
front_MPF <- nsga2(obj_function, 2, 2, generations = results_simu$gen[1], popsize = 100,
                   lower.bounds = rep(-10, 2), upper.bounds = rep(10, 2))

set.seed(2503)
front_MGBM <- nsga2(obj_function, 2, 2, generations = results_simu$gen[2], popsize = 100,
                    lower.bounds = rep(-10, 2), upper.bounds = rep(10, 2))
set.seed(2503)
front_LSSC <- nsga2(obj_function, 2, 2, generations = results_simu$gen[4], popsize = 100,
                    lower.bounds = rep(-10, 2), upper.bounds = rep(10, 2))
set.seed(2503)
front_entropy <- nsga2(obj_function, 2, 2, generations = results_simu$gen[5], popsize = 100,
                       lower.bounds = rep(-10, 2), upper.bounds = rep(10, 2))
set.seed(2503)
front_ocd <- nsga2(obj_function, 2, 2, generations = results_simu$gen[3], popsize = 100,
                   lower.bounds = rep(-10, 2), upper.bounds = rep(10, 2))

MPF_plot <- ggplot() + 
  geom_line(data = pareto_test, aes(x=Y1,Y2),size=1) +
  geom_point(aes(x=-front_MPF$value[front_MPF$pareto.optimal, , drop = FALSE][,1],
                 y=-front_MPF$value[front_MPF$pareto.optimal,,drop = FALSE][,2]), col = "firebrick", size = 2, alpha= .5, shape= 16) +
  labs(
    title = "",
    x = "Y1",
    y = "Y2"
  ) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )  + 
  annotate("text", x = .3, y = .25, label = paste("Gen = ", results_simu$gen[1]), 
           hjust = -0.1, vjust = 1.1, size = 4)

MGBM_plot <- ggplot() + 
  geom_line(data = pareto_test, aes(x=Y1,Y2),size=1) +
  geom_point(aes(x=-front_MGBM$value[front_MGBM$pareto.optimal,, drop = FALSE][,1], 
                 y=-front_MGBM$value[front_MGBM$pareto.optimal, , drop = FALSE][,2]), col = "#F194B4", size = 2, alpha= .5, shape= 16) +
  labs(
    title = "",
    x = "Y1",
    y = "Y2"
  ) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("text", x = .3, y = .25, label = paste("Gen = ", results_simu$gen[2]), 
           hjust = -0.1, vjust = 1.1, size = 4)

LSSC_plot <- ggplot() + 
  geom_line(data = pareto_test, aes(x=Y1,Y2),size=1) +
  geom_point(aes(x=-front_LSSC$value[front_LSSC$pareto.optimal,, drop = FALSE][,1], 
                 y=-front_LSSC$value[front_LSSC$pareto.optimal, , drop = FALSE][,2]), col = "darkcyan", size = 2, alpha= .5, shape= 16) +
  labs(
    title = "",
    x = "Y1",
    y = "Y2"
  ) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("text",x = .3, y = .25, label = paste("Gen = ", results_simu$gen[4]), 
           hjust = -0.1, vjust = 1.1, size = 4)

entropy_plot <- ggplot() + 
  geom_line(data = pareto_test, aes(x=Y1,Y2),size=1) +
  geom_point(aes(x=-front_entropy$value[front_entropy$pareto.optimal, , drop = FALSE][,1], 
                 y=-front_entropy$value[front_entropy$pareto.optimal, , drop = FALSE][,2]), col = "goldenrod2", size = 2, alpha= .5, shape= 16) +
  labs(
    title = "",
    x = "Y1",
    y = "Y2"
  ) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("text", x = .3, y = .25, label = paste("Gen = ", results_simu$gen[5]), 
           hjust = -0.1, vjust = 1.1, size = 4)

ocd_plot <- ggplot() + 
  geom_line(data = pareto_test, aes(x=Y1,Y2),size=1) +
  geom_point(aes(x=-front_ocd$value[front_ocd$pareto.optimal,, drop = FALSE][,1], 
                 y=-front_ocd$value[front_ocd$pareto.optimal, , drop = FALSE][,2]), col = "darkorchid", size = 2, alpha= .5, shape= 16) +  
  labs(
    title = "",
    x = "Y1",
    y = "Y2"
  ) + 
  theme_bw() + 
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13),
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("text", x = .3, y = .25, label = paste("Gen = ", results_simu$gen[3]), 
           hjust = -0.1, vjust = 1.1, size = 4)

figure <- ggarrange(ocd_plot, LSSC_plot, MGBM_plot, entropy_plot, MPF_plot,
                    labels = c("OCD_HV", "LSSC", "MGBM", "Entropy", "MPF"),
                    ncol = 3, nrow = 2, hjust = c(-.8,-1.1,-.9,-.8,-1.5))

########################################
########## Benchmark problem ###########
########################################

results_tot <- read_rds("data/benchmark/results_benchmark_problems.rds")

# Results for two-objective problems : Sect 4.2.1
## Figure 4
pb_2Y <-  results_tot %>% select(HV,critere,gen,time,time_tot,rep,Problem,Spread,nb_obj) %>% filter(nb_obj == 2) 

pb_2Y$time <- log(pb_2Y$time)
pb_2Y$time_tot <- log(pb_2Y$time_tot)
pb_2Y <- pb_2Y %>%  filter(gen > 2)
pb_2Y_long <- pb_2Y %>% select(HV,critere,gen,time,time_tot,rep,Problem,Spread) %>% pivot_longer(!c(rep,critere,Problem), names_to = "Variables", values_to = "Valeurs")

supp.labs <- c("Stopping Generation", "Hypervolume", "Spread" ,"Criterion Time", "Overall Time")
names(supp.labs) <- c("gen", "HV", "Spread", "time","time_tot")

pb_2Y_long %>% 
  ggplot(aes(x = critere, y = Valeurs,fill=Variables)) + 
  geom_boxplot( alpha = 0.3) + 
  facet_grid2(Problem~Variables,scale="free",independent="all", labeller = labeller(Variables = supp.labs)) +  
  labs(
    title = "",
    x = "Stopping criterion",
    y = "Values"
  ) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Dark2") +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

## Figure 5
### ZDT1
set.seed(123)
optim_zdt1 <- mco::nsga2(zdt1, 30,2, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 30), upper.bounds = rep(1,30))
pareto_fronts_zdt <- extract_pareto_front(optim_zdt1)

hv_values_zdt1 <- sapply(pareto_fronts_zdt, function(pareto_fronts_zdt) dominatedHypervolume(pareto_fronts_zdt, rep(10,2)))

set.seed(123)
zdt1_plot <- optim_rep(zdt1,30,2,0,1,1)

p1 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_zdt1[1:500])) +
  theme_bw() + labs(title = "ZDT1, Y = 2", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = zdt1_plot$gen[3], xend = zdt1_plot$gen[3], y = min(hv_values_zdt1)+2.5, yend = max(hv_values_zdt1), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[4], xend = zdt1_plot$gen[4], y = min(hv_values_zdt1)+2.5, yend = max(hv_values_zdt1), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[2], xend = zdt1_plot$gen[2], y = min(hv_values_zdt1)+2.5, yend = max(hv_values_zdt1), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[5], xend = zdt1_plot$gen[5], y = min(hv_values_zdt1)+2.5, yend = max(hv_values_zdt1), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[1], xend = zdt1_plot$gen[1], y = min(hv_values_zdt1)+2.5, yend = max(hv_values_zdt1), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = zdt1_plot$gen[3], y = min(hv_values_zdt1) +0.5, label = "OCD HV", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[4], y = min(hv_values_zdt1)+0.5 , label = "LSSC", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[2], y = min(hv_values_zdt1)+0.5 , label = "MGBM", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[5], y = min(hv_values_zdt1) +0.5, label = "Entropy", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[1], y = min(hv_values_zdt1) +0.5, label = "MPF", size = 4, angle = 90, vjust = 0.2)

### WFG2
set.seed(123)
optim_wfg2 <- mco::nsga2(wfg2_2Y, 20,2, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = rep(1,20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg2)

hv_values_wfg2 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(10,2)))

set.seed(123)
wfg2_plot <- optim_rep(wfg2_2Y,20,2,0,1,1)

p2 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_wfg2[1:500])) +
  theme_bw() + labs(title = "WFG2, Y = 2", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment",x = wfg2_plot$gen[3], xend = wfg2_plot$gen[3], y = min(hv_values_wfg2)+.3, yend = max(hv_values_wfg2), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[4], xend = wfg2_plot$gen[4], y = min(hv_values_wfg2)+.3, yend = max(hv_values_wfg2), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[2], xend = wfg2_plot$gen[2], y = min(hv_values_wfg2)+.3, yend = max(hv_values_wfg2), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[5], xend = wfg2_plot$gen[5], y = min(hv_values_wfg2)+.3, yend = max(hv_values_wfg2), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[1], xend = wfg2_plot$gen[1], y = min(hv_values_wfg2)+.3, yend = max(hv_values_wfg2), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg2_plot$gen[3], y = min(hv_values_wfg2)+.05 , label = "OCD HV", size = 4, angle = 90, vjust = -0.2) +
  annotate("text", x = wfg2_plot$gen[4], y = min(hv_values_wfg2)+.05 , label = "LSSC", size = 4, angle = 90, vjust = 0.5) +
  annotate("text", x = wfg2_plot$gen[2], y = min(hv_values_wfg2)+.05 , label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text",  x = wfg2_plot$gen[5], y = min(hv_values_wfg2)+.05 , label = "Entropy", size = 4, angle = 90, vjust = 0.5) +
  annotate("text", x = wfg2_plot$gen[1], y = min(hv_values_wfg2) +.05, label = "MPF", size = 4, angle = 90, vjust = 0.2)

### WFG3
set.seed(123)
optim_wfg3 <- mco::nsga2(wfg3_2Y, 20,2, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = rep(1,20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg3)

hv_values_wfg3 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(10,2)))

set.seed(123)
wfg3_plot <- optim_rep(wfg3_2Y,20,2,0,1,1)

p3 <- ggplot() + geom_point(aes(x=c(1:750),y=hv_values_wfg3[1:750])) +
  theme_bw() + labs(title = "WFG3, Y = 2", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = wfg3_plot$gen[3], xend = wfg3_plot$gen[3], y = min(hv_values_wfg3)+.25, yend = max(hv_values_wfg3), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[4], xend = wfg3_plot$gen[4], y = min(hv_values_wfg3)+.25, yend = max(hv_values_wfg3), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[2], xend = wfg3_plot$gen[2], y = min(hv_values_wfg3)+.25, yend = max(hv_values_wfg3), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[5], xend = wfg3_plot$gen[5], y = min(hv_values_wfg3)+.25, yend = max(hv_values_wfg3), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[1], xend = wfg3_plot$gen[1], y = min(hv_values_wfg3)+.25, yend = max(hv_values_wfg3), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg3_plot$gen[3], y = min(hv_values_wfg3)+.05 , label = "OCD HV", size = 4, angle = 90, vjust = 0.1) +
  annotate("text", x = wfg3_plot$gen[4], y = min(hv_values_wfg3)+.05 , label = "LSSC", size = 4, angle = 90, vjust = 1) +
  annotate("text", x = wfg3_plot$gen[2], y = min(hv_values_wfg3)+.05 , label = "MGBM", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = wfg3_plot$gen[5], y = min(hv_values_wfg3)+.05 , label = "Entropy", size = 4, angle = 90, vjust = -0.7) +
  annotate("text", x = wfg3_plot$gen[1], y = min(hv_values_wfg3) +.05, label = "MPF", size = 4, angle = 90, vjust = 0.3)

### WFG4
set.seed(123)
optim_wfg4 <- mco::nsga2(wfg4_2Y, 20,2, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = rep(1,20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg4)

hv_values_wfg4 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(10,2)))

set.seed(123)
wfg4_plot <- optim_rep(wfg4_2Y,20,2,0,1,1)

p4 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_wfg4[1:500])) +
  theme_bw() + labs(title = "WFG4, Y = 2", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = wfg4_plot$gen[3], xend = wfg4_plot$gen[3], y = min(hv_values_wfg4)+1.1, yend = max(hv_values_wfg4), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[4], xend = wfg4_plot$gen[4], y = min(hv_values_wfg4)+1.1, yend = max(hv_values_wfg4), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[2], xend = wfg4_plot$gen[2], y = min(hv_values_wfg4)+1.1, yend = max(hv_values_wfg4), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[5], xend = wfg4_plot$gen[5], y = min(hv_values_wfg4)+1.1, yend = max(hv_values_wfg4), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[1], xend = wfg4_plot$gen[1], y = min(hv_values_wfg4)+1.1, yend = max(hv_values_wfg4), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg4_plot$gen[3], y = min(hv_values_wfg4)+.2 , label = "OCD HV", size = 3.6, angle = 90, vjust = -0.5) +
  annotate("text", x = wfg4_plot$gen[4], y = min(hv_values_wfg4)+.2 , label = "LSSC", size = 3.6, angle = 90, vjust = 0.8) +
  annotate("text", x = wfg4_plot$gen[2], y = min(hv_values_wfg4)+.2 , label = "MGBM", size = 3.6, angle = 90, vjust = 0.5) +
  annotate("text", x = wfg4_plot$gen[5], y = min(hv_values_wfg4)+.2 , label = "Entropy", size = 3.6, angle = 90, vjust = -0.2) +
  annotate("text", x = wfg4_plot$gen[1], y = min(hv_values_wfg4) +.2, label = "MPF", size = 4, angle = 90, vjust = 0.3)


(p1 | p2) / (p3 |p4) 

# Results for four-objective problems : Sect 4.2.2
## Figure 6
pb_4Y <-  results_tot %>% select(HV,critere,gen,time,time_tot,rep,Problem,Spread,nb_obj) %>% filter(nb_obj == 4)
pb_4Y$time <- log(pb_4Y$time)
pb_4Y$time_tot <- log(pb_4Y$time_tot)

pb_4Y_long <- pb_4Y %>% select(HV,critere,gen,time,time_tot,rep,Problem,Spread) %>% pivot_longer(!c(rep,critere,Problem), names_to = "Variables", values_to = "Valeurs")

supp.labs <- c("Stopping Generation", "Hypervolume", "Spread" ,"Criterion Time", "Overall Time")
names(supp.labs) <- c("gen", "HV", "Spread", "time","time_tot")

pb_4Y_long %>% 
  ggplot(aes(x = critere, y = Valeurs,fill=Variables)) + 
  geom_boxplot( alpha = 0.3) + 
  facet_grid2(Problem~Variables,scale="free",independent="all", labeller = labeller(Variables = supp.labs)) +  
  labs(
    title = "",
    x = "Stopping criterion",
    y = "Values"
  ) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Dark2") +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 


## Figure 7
### WFG2
set.seed(2503)
optim_wfg2_4Y <- mco::nsga2(wfg2_4Y, 20,4, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = rep(1,20))

pareto_fronts_wfg2_4Y <- extract_pareto_front(optim_wfg2_4Y)

hv_values_wfg2_4Y <- sapply(pareto_fronts_wfg2_4Y, function(pareto_fronts_wfg2_4Y) dominatedHypervolume(pareto_fronts_wfg2_4Y, rep(10,4)))

set.seed(2503)
wfg2_plot_4Y <- optim_rep(wfg2_4Y,20,4,0,1,1)

p2 <- ggplot() + geom_point(aes(x=c(1:800),y=hv_values_wfg2_4Y[1:800])) +
  theme_bw() + labs(title = "WFG2, Y = 4", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = wfg2_plot_4Y$gen[4], xend = wfg2_plot_4Y$gen[4], y = min(hv_values_wfg2_4Y)+20, yend = max(hv_values_wfg2_4Y), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot_4Y$gen[2], xend = wfg2_plot_4Y$gen[2], y = min(hv_values_wfg2_4Y)+20, yend = max(hv_values_wfg2_4Y), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot_4Y$gen[5], xend = wfg2_plot_4Y$gen[5], y = min(hv_values_wfg2_4Y)+20, yend = max(hv_values_wfg2_4Y), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot_4Y$gen[1], xend = wfg2_plot_4Y$gen[1], y = min(hv_values_wfg2_4Y) + 20, , yend = max(hv_values_wfg2_4Y), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg2_plot_4Y$gen[4], y = min(hv_values_wfg2_4Y) + 5, label = "LSSC", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = wfg2_plot_4Y$gen[2], y = min(hv_values_wfg2_4Y) + 5, label = "MGBM", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = wfg2_plot_4Y$gen[5], y = min(hv_values_wfg2_4Y) + 5, label = "Entropy", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = wfg2_plot_4Y$gen[1], y = min(hv_values_wfg2_4Y) +5, label = "MPF", size = 4, angle = 90, vjust = 0.2)


### WFG3
set.seed(123)
optim_wfg3 <- mco::nsga2(wfg3_4Y, 20,4, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = rep(1,20))

pareto_fronts_wfg3_4Y <- extract_pareto_front(optim_wfg3)

hv_values_wfg3_4Y <- sapply(pareto_fronts_wfg3_4Y, function(pareto_fronts_wfg3_4Y) dominatedHypervolume(pareto_fronts_wfg3_4Y, rep(10,4)))

set.seed(123)
wfg3_plot_4Y <- optim_rep(wfg3_4Y,20,4,0,1,1)

p3 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_wfg3_4Y[1:500])) +
  theme_bw() + labs(title = "WFG3, Y = 4", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = wfg3_plot_4Y$gen[4], xend = wfg3_plot_4Y$gen[4], y = min(hv_values_wfg3_4Y)+20, yend = max(hv_values_wfg3_4Y), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot_4Y$gen[2], xend = wfg3_plot_4Y$gen[2], y = min(hv_values_wfg3_4Y)+20, yend = max(hv_values_wfg3_4Y), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot_4Y$gen[5], xend = wfg3_plot_4Y$gen[5], y = min(hv_values_wfg3_4Y)+20, yend = max(hv_values_wfg3_4Y), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot_4Y$gen[1], xend = wfg3_plot_4Y$gen[1], y = min(hv_values_wfg3_4Y)+20, yend = max(hv_values_wfg3_4Y), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg3_plot_4Y$gen[4], y = min(hv_values_wfg3_4Y) + 5, label = "LSSC", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot_4Y$gen[2], y = min(hv_values_wfg3_4Y) + 5, label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot_4Y$gen[5], y = min(hv_values_wfg3_4Y) + 5, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot_4Y$gen[1], y = min(hv_values_wfg3_4Y) +5, label = "MPF", size = 4, angle = 90, vjust = 0.3)



### WFG4
set.seed(2503)
optim_wfg4_4Y <- mco::nsga2(wfg4_4Y, 20,4, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = rep(1,20))

pareto_fronts_wfg4_4Y <- extract_pareto_front(optim_wfg4_4Y)

hv_values_wfg4_4Y <- sapply(pareto_fronts_wfg4_4Y, function(pareto_fronts_wfg4_4Y) dominatedHypervolume(pareto_fronts_wfg4_4Y, rep(10,4)))

set.seed(2503)
wfg4_plot_4Y <- optim_rep(wfg4_4Y,20,4,0,1,1)

p4 <- ggplot() + geom_point(aes(x=c(1:5000),y=hv_values_wfg4_4Y[1:5000])) +
  theme_bw() + labs(title = "WFG4, Y = 4", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = wfg4_plot_4Y$gen[2], xend = wfg4_plot_4Y$gen[2], y = min(hv_values_wfg4_4Y)+200, yend = max(hv_values_wfg4_4Y), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot_4Y$gen[5], xend = wfg4_plot_4Y$gen[5], y = min(hv_values_wfg4_4Y)+200, yend = max(hv_values_wfg4_4Y), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot_4Y$gen[1], xend = wfg4_plot_4Y$gen[1], y = min(hv_values_wfg4_4Y)+200, yend = max(hv_values_wfg4_4Y), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg4_plot_4Y$gen[2], y = min(hv_values_wfg4_4Y) + 50, label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot_4Y$gen[5], y = min(hv_values_wfg4_4Y) + 50, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot_4Y$gen[1], y = min(hv_values_wfg4_4Y) + 50, label = "MPF", size = 4, angle = 90, vjust = 0.3)

p5 <- ggplot() + 
  theme_void() + 
  theme(panel.background = element_rect(fill = "white", color = NA))

(p2 | p3) / (p4 |p5)

########################################
########## Industrial case #############
########################################
# Evaluation of the stopping criteria on a cheese-making process optimization problem : Sect 4.3
results_real_case <- read_rds("data/industrial_case/results_industrial_case.rds")

## Figure 8
results_real_case$time <- log(results_real_case$time)
results_real_case$time_tot <- log(results_real_case$time_tot)

results_real_case_long <- results_real_case %>% select(HV,critere,gen,time,time_tot,rep,Spread) %>% pivot_longer(!c(rep,critere), names_to = "Variables", values_to = "Valeurs")

supp.labs <- c("Stopping Generation", "Hypervolume", "Spread" ,"Criterion Time", "Overall Time")
names(supp.labs) <- c("gen", "HV", "Spread", "time","time_tot")

results_real_case_long %>% 
  ggplot(aes(x = critere, y = Valeurs,fill=Variables)) + 
  geom_boxplot( alpha = 0.3) + 
  facet_grid2(~Variables,scale="free",independent="all", labeller = labeller(Variables = supp.labs)) +  
  labs(
    title = "",
    x = "Stopping criterion",
    y = "Values"
  ) + 
  theme_bw() + 
  scale_fill_brewer(palette = "Dark2") +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

## Figure 9
Front_real <- read_rds("data/industrial_case/results_HV.rds")

Front_real <- extract_pareto_front(Front_real)

### Delate the first iterations without solutions
Front_real <- Front_real[18:1500]

hv_values_real <- sapply(Front_real, function(Front_real) dominatedHypervolume(Front_real, c(1,330,310,7)))


ggplot() + geom_point(aes(x=c(1:1000),y=hv_values_real[1:1000])) +
  theme_bw() + labs(title = "Industrial case, Y = 4", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = results_real_case$gen[498]-17, xend = results_real_case$gen[498]-17, y = min(hv_values)+6, yend = max(hv_values), colour = "darkorchid", size=1.5, alpha=0.6) + # Don't forget to remove iterations without solutions
  annotate("segment", x = results_real_case$gen[499]-17, xend = results_real_case$gen[499]-17, y = min(hv_values)+6, yend = max(hv_values), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = results_real_case$gen[497]-17, xend = results_real_case$gen[497]-17, y = min(hv_values)+6, yend = max(hv_values), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = results_real_case$gen[500]-17, xend = results_real_case$gen[500]-17, y = min(hv_values)+6, yend = max(hv_values), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = results_real_case$gen[496]-17, xend = results_real_case$gen[496]-17, y = min(hv_values)+6, yend = max(hv_values), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = results_real_case$gen[498]-17, y = min(hv_values) + 1, label = "OCD HV", size = 4, angle = 90, vjust = 0.9) +
  annotate("text", x = results_real_case$gen[499]-17, y = min(hv_values) + 1, label = "LSSC", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = results_real_case$gen[497]-17, y = min(hv_values) + 1, label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = results_real_case$gen[500]-17, y = min(hv_values) + 1, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = results_real_case$gen[496]-17, y = min(hv_values) + 1, label = "MPF", size = 4, angle = 90, vjust = -0.1)

########################################
############ Future Work ###############
########################################

set.seed(123)
wfg1_2Y <- makeWFG1Function(2,6,14)

set.seed(123)
optim_wfg1 <- mco::nsga2(wfg1_2Y, 20,2, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = rep(1,20))

pareto_fronts_wfg1 <- extract_pareto_front(optim_wfg1)

hv_values_wfg1 <- sapply(pareto_fronts_wfg1, function(pareto_fronts_wfg1) dominatedHypervolume(pareto_fronts_wfg1, rep(10,2)))

set.seed(123)
wfg1_plot_2Y <- optim_rep(wfg1_2Y,20,2,0,1,1)

ggplot() + geom_point(aes(x=c(1:1000),y=hv_values_wfg1[1:1000])) +
  theme_bw() + labs(title = "WFG1, Y = 2", x = "Generations", y="Hypervolume" )   +
  theme(
    plot.title = element_text( face = "bold", size = 16),  
    plot.subtitle = element_text(hjust = 0.5, size = 12),  
    axis.title = element_text(size = 15),  
    axis.text = element_text(size = 13),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  annotate("segment", x = wfg1_plot_2Y$gen[3], xend = wfg1_plot_2Y$gen[3], y = min(hv_values_wfg1)+1.2, yend = max(hv_values_wfg1), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[4], xend = wfg1_plot_2Y$gen[4], y = min(hv_values_wfg1)+1.2, yend = max(hv_values_wfg1), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[2], xend = wfg1_plot_2Y$gen[2], y = min(hv_values_wfg1)+1.2, yend = max(hv_values_wfg1), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[5], xend = wfg1_plot_2Y$gen[5], y = min(hv_values_wfg1)+1.2, yend = max(hv_values_wfg1), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[1], xend = wfg1_plot_2Y$gen[1], y = min(hv_values_wfg1)+1.2, yend = max(hv_values_wfg1), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg1_plot_2Y$gen[3], y = min(hv_values_wfg1) +0.2, label = "OCD HV", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = wfg1_plot_2Y$gen[4], y = min(hv_values_wfg1)+0.2 , label = "LSSC", size = 4, angle = 90, vjust = 1.5) +
  annotate("text", x = wfg1_plot_2Y$gen[2], y = min(hv_values_wfg1)+0.2 , label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg1_plot_2Y$gen[5], y = min(hv_values_wfg1) +0.2, label = "Entropy", size = 4, angle = 90, vjust = -0.3) +
  annotate("text", x = wfg1_plot_2Y$gen[1], y = min(hv_values_wfg1) +0.2, label = "MPF", size = 4, angle = 90, vjust = 0.3)
