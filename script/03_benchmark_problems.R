# Comparative analysis of stopping criteria for multi-objective evolutionary algorithms
# Manon Perrignon

source("script/01_functions.R")

########################################
########## Benchmark problem ###########
########################################
### Results for two-objective problems : Sect 4.2.1
set.seed(2503)
results_zdt1 <- optim_rep(zdt1,30,2,0,rep(1,30),100,10)

set.seed(2503)
wfg2_2Y <- makeWFG2Function(2,6,14)
results_wfg2_2Y <- optim_rep(wfg2_2Y,20,2,0,2*seq_len(20),100,10)

set.seed(2503)
wfg3_2Y <- makeWFG3Function(2,6,14)
results_wfg3_2Y <- optim_rep(wfg3_2Y,20,2,0,2*seq_len(20),100,10)

set.seed(2503)
wfg4_2Y <- makeWFG4Function(2,6,14)
results_wfg4_2Y <- optim_rep(wfg4_2Y,20,2,0,2*seq_len(20),100,10)

### Results for four-objective problems : Sect 4.2.2
set.seed(123)
wfg2_4Y <- makeWFG2Function(4,6,14)
results_wfg2_4Y <- optim_rep(wfg2_4Y,20,4,0,2*seq_len(20),100,10)

set.seed(2503)
wfg3_4Y <- makeWFG3Function(4,6,14)
results_wfg3_4Y <- optim_rep(wfg3_4Y,20,4,0,2*seq_len(20),100,10)

set.seed(2503)
wfg4_4Y <- makeWFG4Function(4,6,14)
results_wfg4_4Y <- optim_rep(wfg4_4Y,20,4,0,2*seq_len(20),100,10)

### Results for six-objective problems : not included
set.seed(123)
wfg2_6Y <- makeWFG2Function(6,10,10)
results_wfg2_6Y <- optim_rep(wfg2_6Y,20,6,0,2*seq_len(20),100,25)

set.seed(2503)
wfg3_6Y <- makeWFG3Function(6,10,10)
results_wfg3_6Y <- optim_rep(wfg3_6Y,20,6,0,2*seq_len(20),100,25)

set.seed(2503)
wfg4_6Y <- makeWFG4Function(6,10,10)
results_wfg4_6Y <- optim_rep(wfg4_6Y,20,6,0,2*seq_len(20),100,25)

### Results for eight-objective problems : Sect 4.2.3
set.seed(123)
wfg2_8Y <- makeWFG2Function(8,14,6)
results_wfg2_8Y <- optim_rep(wfg2_8Y,20,8,0,2*seq_len(20),100,25)

set.seed(2503)
wfg3_8Y <- makeWFG3Function(8,14,6)
results_wfg3_8Y <- optim_rep(wfg3_8Y,20,8,0,2*seq_len(20),100,25)

set.seed(2503)
wfg4_8Y <- makeWFG4Function(8,14,6)
results_wfg4_8Y <- optim_rep(wfg4_8Y,20,8,0,2*seq_len(20),100,25)


results_tot <- bind_rows(results_zdt1,results_wfg2_2Y,results_wfg3_2Y,results_wfg4_2Y,
                         results_wfg2_4Y,results_wfg3_4Y,results_wfg4_4Y,
                         results_wfg2_8Y,results_wfg3_8Y,results_wfg4_8Y,
                         .id = "Problem")


results_tot$Problem <- results_tot$Problem %>%  fct_recode("ZDT1" = "1",  "WFG2" = "2",
                                                           "WFG3" = "3", "WFG4" = "4","WFG2" = "5", "WFG3" = "6",
                                                           "WFG4" = "7","WFG2" = "8", "WFG3" = "9",
                                                           "WFG4" = "10")


results_tot$critere <- fct_relevel(results_tot$critere, c("OCD_HV", "LSSC", "MGBM","Entropy","MPF"))

results_tot <- results_tot %>%
  mutate(
    HV     = ifelse(gen == 5000, NA, HV),
    Spread = ifelse(gen == 5000, NA, Spread),
    gen    = ifelse(gen == 5000, NA, gen)
  )


write_rds(results_tot,"data/benchmark/results_benchmark_problems.rds")