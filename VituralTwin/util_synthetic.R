get_gt_bool = function(dd, simulator_name){
  # simulator_name is the name of the simulator
  if(simulator_name == "simulate1"){
    return(dd$X1 > 1)
  } else if(simulator_name == "simulate_imbalance_treatment"){
    return(dd$X1 > 1)
  } else if(simulator_name == "simulate_long_rule"){
    return(dd$X1 > -1 & dd$X2 > -1 & dd$X3 > -1)
  } else {
    stop("Unknown simulator name")
  }
}

gt_te_per_sample <- function(d, simulator_name){
  if (simulator_name == "simulate1" || simulator_name == "simulate_imbalance_treatment") {
    # Same outcome mechanism; only treatment assignment differs.
    mu1 <- ifelse(d$X1 > 1, 0.8, 0.2)       # E[Y | T=1, X]
    mu0 <- ifelse(d$X1 < -1, 0.75, 0.2)     # E[Y | T=0, X]
    tau <- mu1 - mu0
    return(tau)
  } else if (simulator_name == "simulate_long_rule") {
    # Recreate the rules exactly as in the simulator
    rule1 <- (d$X1 > -1) & (d$X2 > -1) & (d$X3 > -1)
    rule2 <- (d$X1 > -1) & (d$X2 > -1) & (!rule1)
    rule3 <- (d$X1 > -1) & (!rule1) & (!rule2)
    
    # E[Y | T=1, X] by rule; baseline is 0.2 when no rule fires
    mu1 <- ifelse(rule1, 0.8,
                  ifelse(rule2, 0.6,
                         ifelse(rule3, 0.4, 0.2)))
    
    # Under T=0 everything is baseline 0.2 in this simulator
    mu0 <- rep(0.2, nrow(d))
    
    tau <- mu1 - mu0
    return(tau)
  } else {
    stop("Unknown simulator_name: ", simulator_name)
  }
}

gt_rule = list(
  simulator1 = "X1 > 1",
  simulator_imbalance_treatment = "X1 > 1",
  simulator_long_rule = "X1 > -1 & X2 > -1 & X3 > -1 & X4 > -1"
)