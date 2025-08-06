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

gt_rule = list(
  simulator1 = "X1 > 1",
  simulator_imbalance_treatment = "X1 > 1",
  simulator_long_rule = "X1 > -1 & X2 > -1 & X3 > -1 & X4 > -1"
)