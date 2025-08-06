read_and_process_semi_data = function(path){
  d_original = read.csv(path)
  gt_te = rep(NA, nrow(d_original))
  gt_te[d_original$treatment == 1] = 
    d_original$y_factual[d_original$treatment == 1] - 
    d_original$y_cfactual[d_original$treatment == 1]
  gt_te[d_original$treatment == 0] =
    d_original$y_cfactual[d_original$treatment == 0] - 
    d_original$y_factual[d_original$treatment == 0]
  
  d <- d_original[ , !(colnames(d_original) %in% c("y_factual", "y_cfactual")) ]
  d$Y = d_original$y_factual
  return(list(d=d, gt_te=gt_te))
}

