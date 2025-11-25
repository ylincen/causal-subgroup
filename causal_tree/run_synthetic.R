rm(list = ls())
library(causalTree)


source("./util_synthetic.R")
# load the datasets
simulator_names = c(
  "simulate1",
  "simulate_imbalance_treatment",
  "simulate_long_rule"
)



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
    cat("running CausalTree on simulator: ", simulator_name, " with n = ", n, "\n")
    for(iter_ in iters){
      cat("iteration: ", iter_, "\n")
      counter = counter + 1
      set.seed(counter)
      
      file_path = paste0("../datasets/synthetic/", simulator_name, "/n_", n, "_iter_", iter_, ".csv")
      d = read.csv(file_path)
      
      train_indices = sample(1:n, n * 0.5, replace = FALSE)
      d_train = d[train_indices, ]
      d_test = d[-train_indices, ]
      
      # fit the causalTree model
      d_train_noT = d_train[, !colnames(d_train) %in% c("T")]
      d_test_noT = d_test[, !colnames(d_test) %in% c("T")]
      
      gt_max_bool = get_gt_bool(d_test, simulator_name)
      gt_te =
        mean(d_test$Y[gt_max_bool & (d_test$T == 1)]) - mean(d_test$Y[gt_max_bool & (d_test$T == 0)])
      
      
      tree_model <- causalTree(Y~., data = d_train_noT, treatment = d_train$T,
                         split.Rule = "CT", cv.option = "CT", split.Honest = T, cv.Honest = T, split.Bucket = F, 
                         xval = 5, cp = 0, minsize = 20)
      opcp <- tree_model$cptable[,1][which.min(tree_model$cptable[,4])]
      
      opfit <- prune(tree_model, opcp)
      tree_model = opfit
      
      tree_rules = rpart.rules(tree_model, cover=F)
      num_leaves = sum(tree_model$frame$var =="<leaf>")
      tree_treatment_effects = rep(0, nrow(tree_rules))
      tree_max_te = -Inf
      for(i in 1:num_leaves){
        rs = do.call(paste, tree_rules[i,])
        rs = strsplit(rs, "when")[[1]][2]
        if(is.na(rs)){
          bool_test = rep(T, nrow(d_test))
        } else{
          rs_vec = strsplit(rs, split="&")[[1]]
          bool_test = rep(T, nrow(d_test))
          for(j in 1:length(rs_vec)){
            r = rs_vec[j]
            
            operator <- regmatches(r, regexpr("(>=|<=|>|<|=)", r))
            if(length(operator) == 0){
              operator = "in"  
              feature = strsplit(r, "is")[[1]][1]
              feature = gsub(" ", "", feature)
              
              values = strsplit(r, "is")[[1]][2]
              values = as.numeric(strsplit(values, "to")[[1]])
              if(feature == "T"){
                next
              }
              bool_test = bool_test & (d_test[[feature]] >= values[1] & d_test[[feature]] <= values[2])
            } else{
              feature_value = strsplit(r, operator)[[1]]
              feature = feature_value[1]   
              feature = gsub(" ", "", feature)
              if(feature == "T"){
                next
              }
              value = as.numeric(feature_value[2])
              bool_test = bool_test & switch(
                operator,
                ">" = d_test[[feature]] > value,
                "<" = d_test[[feature]] < value,
                ">=" = d_test[[feature]] >= value,
                "<=" = d_test[[feature]] <= value,
                "=" = d_test[[feature]] == value
              )
            }
          }
        }
        if((sum(bool_test & (d_test$T == 1))==0) | 
           (sum(bool_test & (d_test$T == 0))==0)){
          next
        }
        te = mean(d_test$Y[bool_test & (d_test$T == 1)]) - mean(d_test$Y[bool_test & (d_test$T == 0)])
        tree_treatment_effects[i] = te
        if(tree_max_te < te){
          tree_max_te = te
          tree_max_sg_bool = bool_test
        }
      }
      which_max = which.max(tree_treatment_effects)
      
      jaccard_similarity = sum(gt_max_bool & tree_max_sg_bool) / sum(gt_max_bool | tree_max_sg_bool)
      
      
      gt_te_theoretical_per_sample = gt_te_per_sample(d_test, simulator_name)
      gt_te_of_gt_subgroup_ = mean(gt_te_theoretical_per_sample[gt_max_bool])
      gt_te_of_learned_subgroup_ = mean(gt_te_theoretical_per_sample[tree_max_sg_bool])
      estimated_te_of_learned_subgroup_ = mean(d_test$Y[tree_max_sg_bool & (d_test$T == 1)]) - 
        mean(d_test$Y[tree_max_sg_bool & (d_test$T == 0)])
      estimated_te_of_gt_subgroup_ = mean(d_test$Y[gt_max_bool & (d_test$T == 1)]) -
        mean(d_test$Y[gt_max_bool & (d_test$T == 0)])
      
      
      full_results[counter, ] = list(simulator_name, n, iter_, jaccard_similarity, 
                                     tree_treatment_effects[which_max], gt_te,
                                     diff_te = abs(tree_treatment_effects[which_max] - gt_te),
                                     gt_te_of_gt_subgroup=gt_te_of_gt_subgroup_,
                                     estimated_te_of_gt_subgroup=estimated_te_of_gt_subgroup_,
                                     gt_te_of_learned_subgroup=gt_te_of_learned_subgroup_,
                                     estimated_te_of_learned_subgroup=estimated_te_of_learned_subgroup_
                                     )
                                     
    }

    
    write.csv(full_results, "res_CausalTree_simulation.csv", row.names = FALSE)
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

write.csv(summary_res, "summary_res_CausalTree_simulations.csv", row.names = FALSE)
