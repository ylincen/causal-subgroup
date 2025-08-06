# This script collects the simulations we are gonna use
simulate1 = function(n, iter){
  set.seed(iter)
  X1 = rnorm(n, 0, 1)
  X2 = rnorm(n, 0, 1)
  Treatment = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  
  Y = rep(0, n)
  Y[(X1 > 1) & (Treatment == 1)] = rbinom(sum(X1 > 1 & Treatment == 1), 1, 0.8)
  Y[(X1 < -1) & (Treatment == 0)] = rbinom(sum(X1 < -1 & Treatment == 0), 1, 0.75)
  otherwise_bool = !(((X1 > 1) & (Treatment == 1)) | ((X1 < -1) & (Treatment == 0)))
  Y[otherwise_bool] = rbinom(sum(otherwise_bool), 1, 0.2)
  
  stopifnot(sum((X1 > 1) & (Treatment == 1)) + 
              sum((X1 < -1) & (Treatment == 0)) + 
              sum(otherwise_bool) == n)
  
  d = data.frame(
    X1 = X1,
    X2 = X2,
    T = Treatment,
    Y = Y
  )
  return(d)
}

simulate_imbalance_treatment = function(n, iter){
  set.seed(iter)
  X1 = rnorm(n, 0, 1)
  X2 = rnorm(n, 0, 1)
  
  Treatment = rep(0, n)
  Treatment[X1 >= 0] = sample(c(0, 1), sum(X1 > 0), replace = TRUE, prob = c(0.8, 0.2))
  Treatment[X1 < 0] = sample(c(0, 1), sum(X1 < 0), replace = TRUE, prob = c(0.2, 0.8))
  
  Y = rep(0, n)
  Y[(X1 > 1) & (Treatment == 1)] = rbinom(sum(X1 > 1 & Treatment == 1), 1, 0.8)
  Y[(X1 < -1) & (Treatment == 0)] = rbinom(sum(X1 < -1 & Treatment == 0), 1, 0.75)
  otherwise_bool = !(((X1 > 1) & (Treatment == 1)) | ((X1 < -1) & (Treatment == 0)))
  Y[otherwise_bool] = rbinom(sum(otherwise_bool), 1, 0.2)
  
  stopifnot(sum((X1 > 1) & (Treatment == 1)) + 
              sum((X1 < -1) & (Treatment == 0)) + 
              sum(otherwise_bool) == n)
  
  d = data.frame(
    X1 = X1,
    X2 = X2,
    T = Treatment,
    Y = Y
  )
  return(d)
}

simulate_long_rule = function(n, iter){
  set.seed(iter)
  X1 = rnorm(n, 0, 1)
  X2 = rnorm(n, 0, 1)
  X3 = rnorm(n, 0, 1)
  X4 = rnorm(n, 0, 1)
  X5 = rnorm(n, 0, 1)
  Treatment = sample(c(0, 1), n, replace = TRUE, prob = c(0.5, 0.5))
  
  
  rule1 = (X1 > -1) & (X2 > -1) & (X3 > -1)
  rule2 = (X1 > -1) & (X2 > -1) & (!rule1)
  rule3 = (X1 > -1) & (!rule1) & (!rule2)
  
  Y = rep(0, n)
  Y[rule3 & (Treatment == 1)] = rbinom(sum(rule3 & Treatment == 1), 1, 0.4)
  Y[rule2 & (Treatment == 1)] = rbinom(sum(rule2 & Treatment == 1), 1, 0.6)
  Y[rule1 & (Treatment == 1)] = rbinom(sum(rule1 & Treatment == 1), 1, 0.8)
  
  otherwise_bool = !((rule1 & (Treatment == 1)) | 
                       (rule2 & (Treatment == 1)) | 
                       (rule3 & (Treatment == 1)))
  Y[otherwise_bool] = rbinom(sum(otherwise_bool), 1, 0.2)
  
  stopifnot(sum(rule1 & (Treatment == 1)) + 
              sum(rule2 & (Treatment == 1)) + 
              sum(rule3 & (Treatment == 1)) +
              sum(otherwise_bool) == n)
  
  d = data.frame(
    X1 = X1,
    X2 = X2,
    X3 = X3,
    X4 = X4,
    X5 = X5,
    T = Treatment,
    Y = Y
  )
  return(d)
}


ns = seq(500,5000,500)
iters = 1:50
simulators = list(
  simulate1 = simulate1,
  simulate_imbalance_treatment = simulate_imbalance_treatment,
  simulate_long_rule = simulate_long_rule
)
save_path = "../datasets/synthetic/"
for(simulator_name in names(simulators)) {
  simulator = simulators[[simulator_name]]
  save_path_ = paste0(save_path, simulator_name, "/")
  # create dir
  if(!dir.exists(save_path_)) {
    dir.create(save_path_, recursive = TRUE)
  }
  for(n in ns) {
    for(iter in iters) {
      cat("Processing simulator:", simulator_name, 
          "n:", n, "iter:", iter, "\n")
      d = simulator(n, iter)
      write.csv(d, 
                file = paste0(save_path_, "n_", n, "_iter_", iter, ".csv"), 
                row.names = FALSE, 
                quote = FALSE)
    }
  }
}