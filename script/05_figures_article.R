# Comparative analysis of stopping criteria for multi-objective evolutionary algorithms
# Manon Perrignon

source("script/01_functions.R")

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
  xlim(-10,10) + ylim(-10,10) +
  labs(
    x     = expression(X[1]),
    y     = expression(X[2]),
    title = expression(bold(Y[1]~"based on"~X[1]~"and"~X[2]))
  ) +
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
  xlim(-10,10) + ylim(-10,10) +
  theme_bw() + labs(
    x     = expression(X[1]),
    y     = expression(X[2]),
    title = expression(bold(Y[2]~"based on"~X[1]~"and"~X[2]))
  )  +
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
  theme_bw() +
  labs(x = expression(Y[1]),
       y = expression(Y[2]),
       title = "Approximate Pareto Front")  +
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
    axis.title = element_text(size = 17),  
    axis.text.x = element_text(size = 13, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 13),
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
    x = expression(Y[1]),
    y = expression(Y[2])
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
    x = expression(Y[1]),
    y = expression(Y[2])
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
    x = expression(Y[1]),
    y = expression(Y[2])
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
    x = expression(Y[1]),
    y = expression(Y[2])
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
    x = expression(Y[1]),
    y = expression(Y[2])
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
    axis.text.x = element_text(size = 13, face = "bold"),  
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
set.seed(2503)
optim_zdt1 <- mco::nsga2(zdt1, 30,2, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 30), upper.bounds = rep(1,30))
pareto_fronts_zdt <- extract_pareto_front(optim_zdt1)

hv_values_zdt1 <- sapply(pareto_fronts_zdt, function(pareto_fronts_zdt) dominatedHypervolume(pareto_fronts_zdt, rep(10,2)))

set.seed(2503)
zdt1_plot <- optim_rep(zdt1,30,2,0,rep(1,30),1)

p1 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_zdt1[1:500])) +
  theme_bw() + labs(title = "ZDT1, k = 2", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = zdt1_plot$gen[3], xend = zdt1_plot$gen[3], y = min(hv_values_zdt1)+3, yend = max(hv_values_zdt1), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[4], xend = zdt1_plot$gen[4], y = min(hv_values_zdt1)+3, yend = max(hv_values_zdt1), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[2], xend = zdt1_plot$gen[2], y = min(hv_values_zdt1)+3, yend = max(hv_values_zdt1), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[5], xend = zdt1_plot$gen[5], y = min(hv_values_zdt1)+3, yend = max(hv_values_zdt1), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = zdt1_plot$gen[1], xend = zdt1_plot$gen[1], y = min(hv_values_zdt1)+3, yend = max(hv_values_zdt1), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = zdt1_plot$gen[3], y = min(hv_values_zdt1) +0.8, label = "OCD HV", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[4], y = min(hv_values_zdt1)+0.8 , label = "LSSC", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[2], y = min(hv_values_zdt1)+0.8 , label = "MGBM", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[5], y = min(hv_values_zdt1) +0.8, label = "Entropy", size = 4, angle = 90, vjust = 0.2) +
  annotate("text", x = zdt1_plot$gen[1], y = min(hv_values_zdt1) +0.8, label = "MPF", size = 4, angle = 90, vjust = 0.2)

### WFG2
set.seed(123)
optim_wfg2 <- mco::nsga2(wfg2_2Y, 20,2, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg2)

hv_values_wfg2 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(10,2)))

set.seed(123)
wfg2_plot <- optim_rep(wfg2_2Y,20,2,0,2*seq_len(20),1)

p2 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_wfg2[1:500])) +
  theme_bw() + labs(title = "WFG2, k = 2", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment",x = wfg2_plot$gen[3], xend = wfg2_plot$gen[3], y = min(hv_values_wfg2)+1.2, yend = max(hv_values_wfg2), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[4], xend = wfg2_plot$gen[4], y = min(hv_values_wfg2)+1.2, yend = max(hv_values_wfg2), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[2], xend = wfg2_plot$gen[2], y = min(hv_values_wfg2)+1.2, yend = max(hv_values_wfg2), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[5], xend = wfg2_plot$gen[5], y = min(hv_values_wfg2)+1.2, yend = max(hv_values_wfg2), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[1], xend = wfg2_plot$gen[1], y = min(hv_values_wfg2)+1.2, yend = max(hv_values_wfg2), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg2_plot$gen[3], y = min(hv_values_wfg2)+.3 , label = "OCD HV", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg2_plot$gen[4], y = min(hv_values_wfg2)+.3 , label = "LSSC", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg2_plot$gen[2], y = min(hv_values_wfg2)+.3 , label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text",  x = wfg2_plot$gen[5], y = min(hv_values_wfg2)+.3 , label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg2_plot$gen[1], y = min(hv_values_wfg2) +.3, label = "MPF", size = 4, angle = 90, vjust = 0.3)

### WFG3
set.seed(2503)
optim_wfg3 <- mco::nsga2(wfg3_2Y, 20,2, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg3)

hv_values_wfg3 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(10,2)))

set.seed(2503)
wfg3_plot <- optim_rep(wfg3_2Y,20,2,0,2*seq_len(20),1)

p3 <- ggplot() + geom_point(aes(x=c(1:750),y=hv_values_wfg3[1:750])) +
  theme_bw() + labs(title = "WFG3, k = 2", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg3_plot$gen[3], xend = wfg3_plot$gen[3], y = min(hv_values_wfg3)+1.4, yend = max(hv_values_wfg3), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[4], xend = wfg3_plot$gen[4], y = min(hv_values_wfg3)+1.4, yend = max(hv_values_wfg3), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[2], xend = wfg3_plot$gen[2], y = min(hv_values_wfg3)+1.4, yend = max(hv_values_wfg3), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[5], xend = wfg3_plot$gen[5], y = min(hv_values_wfg3)+1.4, yend = max(hv_values_wfg3), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[1], xend = wfg3_plot$gen[1], y = min(hv_values_wfg3)+1.4, yend = max(hv_values_wfg3), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg3_plot$gen[3], y = min(hv_values_wfg3)+.4 , label = "OCD HV", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot$gen[4], y = min(hv_values_wfg3)+.4 , label = "LSSC", size = 4, angle = 90, vjust = -0.2) +
  annotate("text", x = wfg3_plot$gen[2], y = min(hv_values_wfg3)+.4 , label = "MGBM", size = 4, angle = 90, vjust = 0.5) +
  annotate("text", x = wfg3_plot$gen[5], y = min(hv_values_wfg3)+.4 , label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot$gen[1], y = min(hv_values_wfg3) +.4, label = "MPF", size = 4, angle = 90, vjust = 0.3)

### WFG4
set.seed(2503)
optim_wfg4 <- mco::nsga2(wfg4_2Y, 20,2, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg4)

hv_values_wfg4 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(10,2)))

set.seed(2503)
wfg4_plot <- optim_rep(wfg4_2Y,20,2,0,2*seq_len(20),1)

p4 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_wfg4[1:500])) +
  theme_bw() + labs(title = "WFG4, k = 2", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg4_plot$gen[3], xend = wfg4_plot$gen[3], y = min(hv_values_wfg4)+1.7, yend = max(hv_values_wfg4), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[4], xend = wfg4_plot$gen[4], y = min(hv_values_wfg4)+1.7, yend = max(hv_values_wfg4), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[2], xend = wfg4_plot$gen[2], y = min(hv_values_wfg4)+1.7, yend = max(hv_values_wfg4), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[5], xend = wfg4_plot$gen[5], y = min(hv_values_wfg4)+1.7, yend = max(hv_values_wfg4), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[1], xend = wfg4_plot$gen[1], y = min(hv_values_wfg4)+1.7, yend = max(hv_values_wfg4), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg4_plot$gen[3], y = min(hv_values_wfg4)+.4 , label = "OCD HV", size = 3.6, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot$gen[4], y = min(hv_values_wfg4)+.4 , label = "LSSC", size = 3.6, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot$gen[2], y = min(hv_values_wfg4)+.4 , label = "MGBM", size = 3.6, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot$gen[5], y = min(hv_values_wfg4)+.4 , label = "Entropy", size = 3.6, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot$gen[1], y = min(hv_values_wfg4) +.4, label = "MPF", size = 4, angle = 90, vjust = 0.3)


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
    axis.text.x = element_text(size = 13, face = "bold"),  
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
set.seed(123)
optim_wfg2_4Y <- mco::nsga2(wfg2_4Y, 20,4, generations = c(1:5000), popsize = 100,
                            lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg2_4Y <- extract_pareto_front(optim_wfg2_4Y)

hv_values_wfg2_4Y <- sapply(pareto_fronts_wfg2_4Y, function(pareto_fronts_wfg2_4Y) dominatedHypervolume(pareto_fronts_wfg2_4Y, rep(10,4)))

set.seed(2503)
wfg2_plot_4Y <- optim_rep(wfg2_4Y,20,4,0,2*seq_len(20),1)

p2 <- ggplot() + geom_point(aes(x=c(1:800),y=hv_values_wfg2_4Y[1:800])) +
  theme_bw() + labs(title = "WFG2, k = 4", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg2_plot_4Y$gen[2], xend = wfg2_plot_4Y$gen[2], y = min(hv_values_wfg2_4Y)+250, yend = max(hv_values_wfg2_4Y), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot_4Y$gen[5], xend = wfg2_plot_4Y$gen[5], y = min(hv_values_wfg2_4Y)+250, yend = max(hv_values_wfg2_4Y), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot_4Y$gen[1], xend = wfg2_plot_4Y$gen[1], y = min(hv_values_wfg2_4Y) + 250, , yend = max(hv_values_wfg2_4Y), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg2_plot_4Y$gen[2], y = min(hv_values_wfg2_4Y) + 50, label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg2_plot_4Y$gen[5], y = min(hv_values_wfg2_4Y) + 50, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg2_plot_4Y$gen[1], y = min(hv_values_wfg2_4Y) +50, label = "MPF", size = 4, angle = 90, vjust = 0.3)


### WFG3
set.seed(2503)
optim_wfg3 <- mco::nsga2(wfg3_4Y, 20,4, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg3_4Y <- extract_pareto_front(optim_wfg3)

hv_values_wfg3_4Y <- sapply(pareto_fronts_wfg3_4Y, function(pareto_fronts_wfg3_4Y) dominatedHypervolume(pareto_fronts_wfg3_4Y, rep(10,4)))

set.seed(2503)
wfg3_plot_4Y <- optim_rep(wfg3_4Y,20,4,0,2*seq_len(20),1)

p3 <- ggplot() + geom_point(aes(x=c(1:500),y=hv_values_wfg3_4Y[1:500])) +
  theme_bw() + labs(title = "WFG3, k = 4", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg3_plot_4Y$gen[2], xend = wfg3_plot_4Y$gen[2], y = min(hv_values_wfg3_4Y)+150, yend = max(hv_values_wfg3_4Y), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot_4Y$gen[5], xend = wfg3_plot_4Y$gen[5], y = min(hv_values_wfg3_4Y)+150, yend = max(hv_values_wfg3_4Y), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot_4Y$gen[1], xend = wfg3_plot_4Y$gen[1], y = min(hv_values_wfg3_4Y)+150, yend = max(hv_values_wfg3_4Y), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg3_plot_4Y$gen[2], y = min(hv_values_wfg3_4Y) + 50, label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot_4Y$gen[5], y = min(hv_values_wfg3_4Y) + 50, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot_4Y$gen[1], y = min(hv_values_wfg3_4Y) +50, label = "MPF", size = 4, angle = 90, vjust = 0.3)



### WFG4
set.seed(123)
optim_wfg4_4Y <- mco::nsga2(wfg4_4Y, 20,4, generations = c(1:5000), popsize = 100,
                            lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg4_4Y <- extract_pareto_front(optim_wfg4_4Y)

hv_values_wfg4_4Y <- sapply(pareto_fronts_wfg4_4Y, function(pareto_fronts_wfg4_4Y) dominatedHypervolume(pareto_fronts_wfg4_4Y, rep(10,4)))

set.seed(123)
wfg4_plot_4Y <- optim_rep(wfg4_4Y,20,4,0,2*seq_len(20),1)

p4 <- ggplot() + geom_point(aes(x=c(1:5000),y=hv_values_wfg4_4Y[1:5000])) +
  theme_bw() + labs(title = "WFG4, k = 4", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg4_plot_4Y$gen[2], xend = wfg4_plot_4Y$gen[2], y = min(hv_values_wfg4_4Y)+250, yend = max(hv_values_wfg4_4Y), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot_4Y$gen[5], xend = wfg4_plot_4Y$gen[5], y = min(hv_values_wfg4_4Y)+250, yend = max(hv_values_wfg4_4Y), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot_4Y$gen[1], xend = wfg4_plot_4Y$gen[1], y = min(hv_values_wfg4_4Y)+250, yend = max(hv_values_wfg4_4Y), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg4_plot_4Y$gen[2], y = min(hv_values_wfg4_4Y) + 50, label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot_4Y$gen[5], y = min(hv_values_wfg4_4Y) + 50, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot_4Y$gen[1], y = min(hv_values_wfg4_4Y) + 50, label = "MPF", size = 4, angle = 90, vjust = 0.3)

p5 <- ggplot() + 
  theme_void() + 
  theme(panel.background = element_rect(fill = "white", color = NA))

(p2 | p3) / (p4 |p5)

# Results for 8-objective problems : Sect 4.2.3
## Figure 8
pb_8Y <-  results_tot %>% select(HV,critere,gen,time,time_tot,rep,Problem,Spread,nb_obj) %>% filter(nb_obj == 8)
pb_8Y$time <- log(pb_8Y$time)
pb_8Y$time_tot <- log(pb_8Y$time_tot)

pb_8Y <- pb_8Y %>%
  mutate(
    HV     = ifelse(gen == 5000, NA, HV),
    Spread = ifelse(gen == 5000, NA, Spread),
    gen    = ifelse(gen == 5000, NA, gen)
  )


pb_8Y_long <- pb_8Y %>% pivot_longer(!c(rep,critere,Problem,nb_obj), names_to = "Variables", values_to = "Valeurs")

supp.labs <- c("Stopping Generation", "Hypervolume", "Spread" ,"Criterion Time", "Overall Time")
names(supp.labs) <- c("gen", "HV", "Spread", "time","time_tot")

pb_8Y_long %>% 
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
    axis.text.x = element_text(size = 13, face = "bold"),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

#Figure 9
### WFG2
set.seed(123)
optim_wfg2 <- mco::nsga2(wfg2_8Y, 20,8, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg2)

hv_values_wfg2 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(25,8)))

set.seed(123)
wfg2_plot <- optim_rep(wfg2_8Y,20,8,0,2*seq_len(20),1)

p2 <- ggplot() + geom_point(aes(x=c(1:1000),y=hv_values_wfg2[1:1000])) +
  theme_bw() + labs(title = "WFG2, k = 8", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg2_plot$gen[2], xend = wfg2_plot$gen[2], y = min(hv_values_wfg2)+0.05e+11, yend = max(hv_values_wfg2), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[5], xend = wfg2_plot$gen[5], y = min(hv_values_wfg2)+0.05e+11, yend = max(hv_values_wfg2), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg2_plot$gen[1], xend = wfg2_plot$gen[1], y = min(hv_values_wfg2)+0.05e+11, yend = max(hv_values_wfg2), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg2_plot$gen[2], y = min(hv_values_wfg2)+0.01e+11 , label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text",  x = wfg2_plot$gen[5], y = min(hv_values_wfg2)+0.01e+11 , label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg2_plot$gen[1], y = min(hv_values_wfg2)+0.01e+11, label = "MPF", size = 4, angle = 90, vjust = 0.3)
p2
### WFG3
set.seed(2503)
optim_wfg3 <- mco::nsga2(wfg3_8Y, 20,8, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg3)

hv_values_wfg3 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(25,8)))

set.seed(2503)
wfg3_plot <- optim_rep(wfg3_8Y,20,8,0,2*seq_len(20),1)

p3 <- ggplot() + geom_point(aes(x=c(1:1000),y=hv_values_wfg3[1:1000])) +
  theme_bw() + labs(title = "WFG3, k = 8", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg3_plot$gen[2], xend = wfg3_plot$gen[2], y = min(hv_values_wfg3)+0.02e+11, yend = max(hv_values_wfg3), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[5], xend = wfg3_plot$gen[5], y = min(hv_values_wfg3)+0.02e+11, yend = max(hv_values_wfg3), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg3_plot$gen[1], xend = wfg3_plot$gen[1], y = min(hv_values_wfg3)+0.02e+11, yend = max(hv_values_wfg3), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg3_plot$gen[2], y = min(hv_values_wfg3)+0.005e+11 , label = "MGBM", size = 4, angle = 90, vjust = .8) +
  annotate("text", x = wfg3_plot$gen[5], y = min(hv_values_wfg3)+0.005e+11 , label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg3_plot$gen[1], y = min(hv_values_wfg3) +0.005e+11, label = "MPF", size = 4, angle = 90, vjust = 0.3)
p3


### WFG4
set.seed(2503)
optim_wfg4 <- mco::nsga2(wfg4_8Y, 20,8, generations = c(1:5000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg <- extract_pareto_front(optim_wfg4)

hv_values_wfg4 <- sapply(pareto_fronts_wfg, function(pareto_fronts_wfg) dominatedHypervolume(pareto_fronts_wfg, rep(25,8)))

set.seed(2503)
wfg4_plot <- optim_rep(wfg4_8Y,20,8,0,2*seq_len(20),1)

p4 <- ggplot() + geom_point(aes(x=c(1:5000),y=hv_values_wfg4[1:5000])) +
  theme_bw() + labs(title = "WFG4, k = 8", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg4_plot$gen[2], xend = wfg4_plot$gen[2], y = min(hv_values_wfg4)+0.05e+11, yend = max(hv_values_wfg4), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[5], xend = wfg4_plot$gen[5], y = min(hv_values_wfg4)+0.05e+11, yend = max(hv_values_wfg4), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg4_plot$gen[1], xend = wfg4_plot$gen[1], y = min(hv_values_wfg4)+0.05e+11, yend = max(hv_values_wfg4), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg4_plot$gen[2], y = min(hv_values_wfg4)+0.02e+11 , label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot$gen[5], y = min(hv_values_wfg4)+0.02e+11 , label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg4_plot$gen[1], y = min(hv_values_wfg4) +0.02e11, label = "MPF", size = 4, angle = 90, vjust = 0.3)
p4

p5 <- ggplot() + 
  theme_void() + 
  theme(panel.background = element_rect(fill = "white", color = NA))

(p2 | p3) / (p4 |p5)


########################################
########## Industrial case #############
########################################
# Evaluation of the stopping criteria on a cheese-making process optimization problem : Sect 4.3
results_real_case <- read_rds("data/industrial_case/results_industrial_case.rds")
## Figure 10
results_real_case$time <- log(results_real_case$time)
results_real_case$time_tot <- log(results_real_case$time_tot)
results_real_case$critere <- fct_relevel(results_real_case$critere, c("OCD_HV", "LSSC", "MGBM","Entropy","MPF"))
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
    axis.text.x = element_text(size = 13, face = "bold"),  
    legend.title = element_blank(),
    legend.text = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 15, face = "bold"),
    strip.background = element_rect(
      color="black", fill="white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) 

## Figure 11
Front_real <- read_rds("data/industrial_case/results_HV.rds")

Front_real <- extract_pareto_front(Front_real)

### Delate the first iterations without solutions
Front_real <- Front_real[18:1500]

hv_values_real <- sapply(Front_real, function(Front_real) dominatedHypervolume(Front_real, c(1,330,310,7)))


ggplot() + geom_point(aes(x=c(1:1000),y=hv_values_real[1:1000])) +
  theme_bw() + labs(title = "Industrial case, k = 4", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = results_real_case$gen[498]-17, xend = results_real_case$gen[498]-17, y = min(hv_values_real)+8.5, yend = max(hv_values_real), colour = "darkorchid", size=1.5, alpha=0.6) + # Don't forget to remove iterations without solutions
  annotate("segment", x = results_real_case$gen[499]-17, xend = results_real_case$gen[499]-17, y = min(hv_values_real)+8.5, yend = max(hv_values_real), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = results_real_case$gen[497]-17, xend = results_real_case$gen[497]-17, y = min(hv_values_real)+8.5, yend = max(hv_values_real), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = results_real_case$gen[500]-17, xend = results_real_case$gen[500]-17, y = min(hv_values_real)+8.5, yend = max(hv_values_real), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = results_real_case$gen[496]-17, xend = results_real_case$gen[496]-17, y = min(hv_values_real)+8.5, yend = max(hv_values_real), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = results_real_case$gen[498]-17, y = min(hv_values_real) + 3, label = "OCD HV", size = 4, angle = 90, vjust = 0.9) +
  annotate("text", x = results_real_case$gen[499]-17, y = min(hv_values_real) + 3, label = "LSSC", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = results_real_case$gen[497]-17, y = min(hv_values_real) + 3, label = "MGBM", size = 4, angle = 90, vjust = -0.3) +
  annotate("text", x = results_real_case$gen[500]-17, y = min(hv_values_real) + 3, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = results_real_case$gen[496]-17, y = min(hv_values_real) + 3, label = "MPF", size = 4, angle = 90, vjust = -0.1)

########################################
############ Future Work ###############
########################################

set.seed(123)
wfg1_2Y <- makeWFG1Function(2,6,14)

set.seed(123)
optim_wfg1 <- mco::nsga2(wfg1_2Y, 20,2, generations = c(1:1000), popsize = 100,
                         lower.bounds = rep(0, 20), upper.bounds = 2*seq_len(20))

pareto_fronts_wfg1 <- extract_pareto_front(optim_wfg1)

hv_values_wfg1 <- sapply(pareto_fronts_wfg1, function(pareto_fronts_wfg1) dominatedHypervolume(pareto_fronts_wfg1, rep(10,2)))

set.seed(123)
wfg1_plot_2Y <- optim_rep(wfg1_2Y,20,2,0,2*seq_len(20),1)

ggplot() + geom_point(aes(x=c(1:1000),y=hv_values_wfg1[1:1000])) +
  theme_bw() + labs(title = "WFG1, k = 2", x = "Generations", y="Hypervolume" )   +
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
  annotate("segment", x = wfg1_plot_2Y$gen[3], xend = wfg1_plot_2Y$gen[3], y = min(hv_values_wfg1)+1.6, yend = max(hv_values_wfg1), colour = "darkorchid", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[4], xend = wfg1_plot_2Y$gen[4], y = min(hv_values_wfg1)+1.6, yend = max(hv_values_wfg1), colour = "darkcyan", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[2], xend = wfg1_plot_2Y$gen[2], y = min(hv_values_wfg1)+1.6, yend = max(hv_values_wfg1), colour = "#F194B4", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[5], xend = wfg1_plot_2Y$gen[5], y = min(hv_values_wfg1)+1.6, yend = max(hv_values_wfg1), colour = "goldenrod2", size=1.5, alpha=0.6) + 
  annotate("segment", x = wfg1_plot_2Y$gen[1], xend = wfg1_plot_2Y$gen[1], y = min(hv_values_wfg1)+1.6, yend = max(hv_values_wfg1), colour = "firebrick", size=1.5, alpha=0.6) +
  annotate("text", x = wfg1_plot_2Y$gen[3], y = min(hv_values_wfg1) +0.4, label = "OCD HV", size = 4, angle = 90, vjust = 0) +
  annotate("text", x = wfg1_plot_2Y$gen[4], y = min(hv_values_wfg1)+0.4 , label = "LSSC", size = 4, angle = 90, vjust = 0.5) +
  annotate("text", x = wfg1_plot_2Y$gen[2], y = min(hv_values_wfg1)+0.4 , label = "MGBM", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg1_plot_2Y$gen[5], y = min(hv_values_wfg1) +0.4, label = "Entropy", size = 4, angle = 90, vjust = 0.3) +
  annotate("text", x = wfg1_plot_2Y$gen[1], y = min(hv_values_wfg1) +0.4, label = "MPF", size = 4, angle = 90, vjust = 0.3)
