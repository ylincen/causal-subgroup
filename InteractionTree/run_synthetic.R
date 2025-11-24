rm(list=ls())
library(dplyr)
source("Functions-IT.R")


source("./util_synthetic.R")
# load the datasets
simulator_names = c(
  "simulate1",
  "simulate_imbalance_treatment",
  "simulate_long_rule"
)

get_rule <- function(id, tab) {
  pieces <- character()
  while (id != 1) {                          # climb to root
    parent <- id %/% 10
    row    <- tab[match(parent, tab$node), ]
    
    dir   <- if (id == parent * 10) "<=" else ">"
    cond  <- sprintf("%s %s %g", row$vname, dir, row$cut)  # now numeric
    pieces <- c(cond, pieces)
    id <- parent
  }
  paste(pieces, collapse = " & ")
}
col_type <- function(x) {
  u <- unique(x)
  
  if (length(u) == 2) {                # exactly two values → binary
    "nominal"
  } else {                             # everything else
    "numeric"
  }
}

# simulate low signal data
num_simulators = length(simulator_names)
ns = seq(1000,5000,500)
iters = 1:50
num_exp_res = length(ns) * length(iters)


gt_te_init = rep(0, num_exp_res)
max_te = rep(0, num_exp_res)
jaccard = rep(0, num_exp_res)
counter = 0

full_results = data.frame(
  simulator_name = rep("", num_exp_res),
  n = rep(0, num_exp_res),
  iter = rep(0, num_exp_res),
  jaccard = jaccard,
  max_te = max_te,
  gt_te = gt_te_init,
  diff_te = rep(0, num_exp_res),
  gt_te_of_gt_subgroup = rep(0, num_exp_res),
  estimated_te_of_gt_subgroup = rep(0, num_exp_res),
  gt_te_of_learned_subgroup = rep(0, num_exp_res),
  estimated_te_of_learned_subgroup = rep(0, num_exp_res)
)


for(simulator_name in simulator_names){
  for(n in ns){
    cat("running InteractionTree on simulator: ", simulator_name, " with n = ", n, "\n")
    for(iter_ in iters){
      counter = counter + 1
      set.seed(counter)
      file_path = paste0("../datasets/synthetic/", simulator_name, "/n_", n, "_iter_", iter_, ".csv")
      d = read.csv(file_path)
      colnames(d)[ncol(d)] <- "y"  # ensure the last column is Y
      colnames(d)[ncol(d)-1] = "trt"
      train_indices = sample(1:n, n * 0.5, replace = FALSE)
      d_train = d[train_indices, ]
      d_test = d[-train_indices, ]
      val_indices = sample(1:nrow(d_train), nrow(d_train) * 0.2, replace = FALSE)
      d_val = d_train[val_indices, ]
      d_train = d_train[-val_indices, ]
      
      
      # fit the model
      split.var = 1:(ncol(d_train) -2)
      cols.ctg = which(sapply(d_train[, split.var], col_type) == "nominal")
      

      tree0 <- grow.INT(data=d_train, test=d_val, min.ndsz=20, n0=5, 
                        split.var=split.var, ctg=cols.ctg, max.depth=5)
      prn <- prune.size.testsample(tree0)

      bsize <- prn$size[4] 
      # OBTAIN THE FINAL TREE STRUCTURE
      if(bsize > 0){
        btree <- obtain.btree(tree0, bsize=bsize); # btree  
      } else{
        btree = tree0
      }
      if(is.null(btree)){
        btree = tree0
      }
      
      btree <- btree %>%                           # your data-frame name
        mutate(
          cut   = as.numeric(as.character(cut)),   # numeric, not factor
          vname = as.character(vname),
          node = as.numeric(node) # plain character
        )
      
      leaves <- btree$node[is.na(btree$vname)]  
      
      if(nrow(btree)==1) leaves = numeric(0)
      leaf_rules <- setNames(lapply(leaves, get_rule, btree), leaves)
      # get the subgroup analysis results
      ## initialize 
      max_te = -Inf
      max_sg_bool = rep(T, nrow(d_test))
      gt_bool = get_gt_bool(dd=d_test, simulator_name=simulator_name)
      gt_te = mean(d_test$y[gt_bool & (d_test$trt == 1)]) - 
        mean(d_test$y[gt_bool & (d_test$trt == 0)])
      
      
      if(length(leaf_rules) == 0){
        bool_test = rep(T, nrow(d_test))
        
        subgroup_te = mean(d_test$y[bool_test & d_test$trt == 1]) - 
          mean(d_test$y[bool_test & d_test$trt == 0])
        if(subgroup_te > max_te){
          max_te = subgroup_te
          max_subgroup_jaccard = sum(bool_test & gt_bool) / 
            sum(bool_test | gt_bool)
          abs_error_te = abs(max_te - gt_te)
        }
        treatment_effects = subgroup_te
      } else{
        treatment_effects = rep(-Inf, length(leaf_rules))
        
        for(i in 1:length(leaf_rules)){
          rule = leaf_rules[[i]]
          items = strsplit(rule, " & ")[[1]]
          bool_test = rep(T, nrow(d_test))
          for(j in 1:length(items)){
            item = items[j]
            if(grepl(">", item)){
              feature_value = strsplit(item, " > ")[[1]]
              operator = ">"
            } else {
              feature_value = strsplit(item, " <= ")[[1]]
              operator = "<="
            }
            feature = feature_value[1]
            value = as.numeric(feature_value[2])
            if(operator == ">"){
              bool_test = bool_test & (d_test[[feature]] > value)
            } else {
              bool_test = bool_test & (d_test[[feature]] <= value)
            }
          }
          
          subgroup_te = mean(d_test$y[bool_test & d_test$trt == 1]) - 
            mean(d_test$y[bool_test & d_test$trt == 0])
          treatment_effects[i] = subgroup_te
          if(subgroup_te > max_te){
            max_te = subgroup_te
            max_subgroup_jaccard = sum(bool_test & gt_bool) / 
              sum(bool_test | gt_bool)
            abs_error_te = abs(max_te - gt_te)
          }
        }  
      }
      # collect the results
      which_max = which.max(treatment_effects)
      jaccard_similarity = sum(gt_bool & max_sg_bool) / sum(gt_bool | max_sg_bool)

      gt_te_theoretical_per_sample = gt_te_per_sample(d_test, simulator_name)
      gt_te_of_gt_subgroup_ = mean(gt_te_theoretical_per_sample[gt_bool])
      gt_te_of_learned_subgroup_ = mean(gt_te_theoretical_per_sample[max_sg_bool])
      estimated_te_of_learned_subgroup_ = mean(d_test$y[max_sg_bool & (d_test$trt == 1)]) - 
        mean(d_test$y[max_sg_bool & (d_test$trt == 0)])
      estimated_te_of_gt_subgroup_ = mean(d_test$y[gt_bool & (d_test$trt == 1)]) -
        mean(d_test$y[gt_bool & (d_test$trt == 0)])
      
      
      full_results[counter, ] = list(simulator_name, n, iter_, jaccard_similarity, 
                                     treatment_effects[which_max], gt_te,
                                     diff_te = abs(treatment_effects[which_max] - gt_te),
                                     gt_te_of_gt_subgroup = gt_te_of_gt_subgroup_,
                                     estimated_te_of_gt_subgroup = estimated_te_of_gt_subgroup_,
                                     gt_te_of_learned_subgroup = gt_te_of_learned_subgroup_,
                                     estimated_te_of_learned_subgroup = estimated_te_of_learned_subgroup_)
    }
    full_results$jaccard = round(full_results$jaccard, 3)
    full_results$max_te = round(full_results$max_te, 3)
    full_results$max_te = round(full_results$max_te, 3)
    full_results$gt_te = round(full_results$gt_te, 3)
    full_results$diff_te = round(full_results$diff_te, 3)
    
    write.csv(full_results, "res_InteractionTree_simulation.csv", row.names = FALSE)
  }
}



require(data.table)
full_res_dt = data.table(full_results)
summary_res = full_res_dt[, .(
  mean_jaccard = mean(jaccard),
  mean_max_te = mean(max_te),
  mean_gt_te = mean(gt_te),
  mean_diff_te = mean(diff_te)
), by = .(simulator_name, n)]

write.csv(summary_res, "summary_res_InteractionTree_simulations.csv", row.names = FALSE)



