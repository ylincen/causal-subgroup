rm(list = ls())
library(rpart)
library(rpart.plot)
require(purrr)
require(R.utils)
require(grf)
set.seed(1)

source("./util_synthetic.R")
source("./causalDT-main/causalDT/R/causalDT.R")
source("./causalDT-main/causalDT/R/causalDT-package.R")
source("./causalDT-main/causalDT/R/student.R")
source("./causalDT-main/causalDT/R/teacher.R")
source("./causalDT-main/causalDT/R/RcppExports.R")
source("./causalDT-main/causalDT/R/diagnostics-stability.R")
source("./causalDT-main/causalDT/R/utils-student.R")
source("./causalDT-main/causalDT/R/utils-crossfit.R")




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
gt_te = rep(0, num_exp_res)
tree_max_te = rep(0, num_exp_res)
tree_jaccard = rep(0, num_exp_res)
rsides_jaccard = rep(0, num_exp_res)
rsides_max_te = rep(0, num_exp_res)
counter = 0

full_results = data.frame(
  simulator_name = rep("", num_exp_res),
  n = rep(0, num_exp_res),
  iter = rep(0, num_exp_res),
  tree_jaccard = tree_jaccard,
  tree_max_te = tree_max_te,
  gt_te = gt_te,
  tree_diff_te = rep(0, num_exp_res),
  gt_te_of_gt_subgroup = rep(0, num_exp_res),
  estimated_te_of_gt_subgroup = rep(0, num_exp_res),
  gt_te_of_learned_subgroup = rep(0, num_exp_res),
  estimated_te_of_learned_subgroup = rep(0, num_exp_res)
)



for(simulator_name in simulator_names){
  for(n in ns){
    for(iter_ in iters){
      counter = counter + 1
      set.seed(counter)
      
      file_path = paste0("../datasets/synthetic/", simulator_name, "/n_", n, "_iter_", iter_, ".csv")
      d = read.csv(file_path)
      
      train_indices = sample(1:n, n * 0.5, replace = FALSE)
      d_train = d[train_indices, ]
      d_test = d[-train_indices, ]
      
      # get the gt info
      gt_max_bool = get_gt_bool(d_test, simulator_name)
      gt_te[counter] =
        mean(d_test$Y[gt_max_bool & (d_test$T == 1)]) - mean(d_test$Y[gt_max_bool & (d_test$T == 0)])
      
      
      
      # Fit a decision tree model
      x = d_train[, -which(names(d_train) %in% c("Y", "T"))]
      y = d_train$Y
      z = d_train$T
      causal_model = causalDT(X=x, Y=y, Z=z, rpart_prune="min")
      tree_model = causal_model$student_fit$fit

      # get the rules
      tree_rules = rpart.rules(tree_model, cover=F, roundint = F)
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
        if(tree_max_te < te){
          tree_max_te = te
          tree_max_sg_bool = bool_test
        }
      }
      which_max = which.max(tree_treatment_effects)
      # jaccard for that subgroup
      jaccard_similarity = sum(gt_max_bool & tree_max_sg_bool) / sum(gt_max_bool | tree_max_sg_bool)
      
      gt_te_theoretical_per_sample = gt_te_per_sample(d_test, simulator_name)
      gt_te_of_gt_subgroup_ = mean(gt_te_theoretical_per_sample[gt_max_bool])
      gt_te_of_learned_subgroup_ = mean(gt_te_theoretical_per_sample[tree_max_sg_bool])
      estimated_te_of_learned_subgroup_ = mean(d_test$Y[tree_max_sg_bool & (d_test$T == 1)]) - 
        mean(d_test$Y[tree_max_sg_bool & (d_test$T == 0)])
      estimated_te_of_gt_subgroup_ = mean(d_test$Y[gt_max_bool & (d_test$T == 1)]) -
        mean(d_test$Y[gt_max_bool & (d_test$T == 0)])
      
      full_results[counter, ] = list(simulator_name, n, iter_, 
                                     jaccard_similarity, tree_max_te, gt_te[counter],
                                     tree_diff_te = abs(tree_max_te - gt_te[counter]),
                                     gt_te_of_gt_subgroup = gt_te_of_gt_subgroup_,
                                     estimated_te_of_gt_subgroup = estimated_te_of_gt_subgroup_,
                                     gt_te_of_learned_subgroup = gt_te_of_learned_subgroup_,
                                     estimated_te_of_learned_subgroup = estimated_te_of_learned_subgroup_) 
                                     
      print(full_results[counter, ,drop=F])
    }
    # round up and save the results
    full_results$tree_jaccard = round(full_results$tree_jaccard, 3)
    full_results$tree_max_te = round(full_results$tree_max_te, 3)
    full_results$gt_te = round(full_results$gt_te, 3)
    full_results$tree_diff_te = round(full_results$tree_diff_te, 3)

    write.csv(full_results, "res_distill.csv", row.names = FALSE)
  }
}
require(data.table)
full_res_dt = data.table(full_results)
# summarize as mean and sd with different n
summary_res = full_res_dt[, .(
  mean_tree_jaccard = mean(tree_jaccard),
  mean_tree_max_te = mean(tree_max_te),
  mean_gt_te = mean(gt_te),
  mean_tree_diff_te = mean(tree_diff_te)
), by = .(n, simulator_name)]
summary_res

write.csv(summary_res, "summary_res_distill.csv", row.names = FALSE)
