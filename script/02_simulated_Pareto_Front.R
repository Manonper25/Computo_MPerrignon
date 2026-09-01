# Comparative analysis of stopping criteria for multi-objective evolutionary algorithms
# Manon Perrignon

source("script/01_functions.R")

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

write_rds(dataXY,"data/simulated_front/data_simulated_front.rds")

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
    predY1 <- predY1 + height * exp(-distance_squared / (2 * sigma^2)) # Gaussian with height
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
    predY2 <- predY2 + height * exp(-distance_squared / (2 * sigma^2)) # Gaussian with height
  }
  
  predictions <- as.matrix(data.frame(Y1 = -predY1,Y2=-predY2))
  return(predictions = predictions)
}

## Results
results_simu <- optim_rep(obj_function,2,2,-10,rep(10,2),1,10)

results_simu$critere <- fct_relevel(results_simu$critere, c("OCD_HV", "LSSC", "MGBM","Entropy", "MPF"))

write_rds(results_simu,"data/simulated_front/results_simulated_front.rds")