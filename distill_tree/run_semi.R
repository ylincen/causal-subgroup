rm(list=ls())

library(rpart)
library(rpart.plot)
require(purrr)
require(R.utils)
require(grf)

source("./causalDT-main/causalDT/R/causalDT.R")
source("./causalDT-main/causalDT/R/causalDT-package.R")
source("./causalDT-main/causalDT/R/student.R")
source("./causalDT-main/causalDT/R/teacher.R")
source("./causalDT-main/causalDT/R/RcppExports.R")
source("./causalDT-main/causalDT/R/diagnostics-stability.R")
source("./causalDT-main/causalDT/R/utils-student.R")
source("./causalDT-main/causalDT/R/utils-crossfit.R")

source("./utils.R")

collect_exp_res = function(tree_model, gt_te_test, d_test){
  tree_rules = rpart.rules(tree_model, cover=F, roundint = F)
  num_leaves = sum(tree_model$frame$var =="<leaf>")
  tree_treatment_effects = rep(0, nrow(tree_rules))
  tree_treatment_effects_var = rep(0, nrow(tree_rules))
  diff_average = rep(0, nrow(tree_rules)) # the difference between the subgroup te and the average gt_tex
  diff_individual <- vector("list", nrow(tree_rules)) # the individual difference of the subgroup_te and the gt_te
  
  tree_max_te = -Inf
  skip_leaves_in_inference = rep(FALSE, num_leaves)
  for(i in 1:num_leaves){
    rs = do.call(paste, tree_rules[i,])
    rs = strsplit(rs, "when")[[1]][2]
    
    if(!grepl("treatment", rs)){
      skip_leaves_in_inference[i] = TRUE
    }
    
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
          feature = strsplit(r, " is ")[[1]][1]
          feature = gsub(" ", "", feature)
          
          values = strsplit(r, " is ")[[1]][2]
          values = as.numeric(strsplit(values, "to")[[1]])
          if(feature == "treatment"){
            next
          }
          if(length(values)==1){
            bool_test = bool_test & (d_test[[feature]] == values[1])
          } else{
            bool_test = bool_test & (d_test[[feature]] >= values[1] & 
                                       d_test[[feature]] <= values[2])  
          }
        } else{
          feature_value = strsplit(r, operator)[[1]]
          feature = feature_value[1]   
          feature = gsub(" ", "", feature)
          if(feature == "treatment"){
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
          if(any(is.na(bool_test))){
            browser()
          }
        }
      }
    }
    
    if((sum(bool_test & (d_test$treatment == 1))==0) | 
       (sum(bool_test & (d_test$treatment == 0))==0)){
      tree_treatment_effects[i] = NA
      next
    }
    
    te = mean(d_test$Y[bool_test & (d_test$treatment == 1)]) - mean(d_test$Y[bool_test & (d_test$treatment == 0)])
    tree_treatment_effects[i] = te
    tree_treatment_effects_var[i] = var(gt_te_test[bool_test])
    diff_average[i] = te - mean(gt_te_test[bool_test])
    diff_individual[[i]] = te - gt_te_test[bool_test] 

    if(tree_max_te < te){
      tree_max_te = te
      tree_max_sg_bool = bool_test  
    }
  }
    
  which_max = which.max(tree_treatment_effects)
  return(list(
    subgroup_treatment_effects = tree_treatment_effects,
    subgroup_gt_treatment_effects_var = tree_treatment_effects_var,
    subgroup_diff_average = diff_average,
    subgroup_diff_individual = diff_individual,
    max_subgroup_treatment_effect = tree_treatment_effects[which_max],
    which_max = which_max,
    max_subgroup_gt_treatment_effct = mean(gt_te_test[tree_max_sg_bool])
  ))
}




semi_datasets_folder = "../datasets/semi_synthetic/"
files = list.files(semi_datasets_folder, pattern = "*.csv", full.names = TRUE)

res_df = data.frame(
  dataset = character(),
  subgroup_treatment_effects = I(list()),
  subgroup_gt_treatment_effects_var = I(list()),
  subgroup_diff_average = I(list()),
  subgroup_diff_individual = I(list()),
  max_subgroup_treatment_effect = numeric(),
  which_max = integer(),
  max_subgroup_gt_treatment_effct = numeric()
)

res_df_summary = data.frame(
  dataset = character(),
  max_subgroup_treatment_effect = numeric(),
  subgroup_diff_average_weighted_mean = numeric(),
  subgroup_diff_average_unweighted_mean = numeric(),
  max_subgroup_diff_average = numeric(),
  max_subgroup_gt_treatment_effct = numeric()
)

for(i in 1:length(files)){
  data_file = files[i]
  cat("running Distll Tree on semi-synthetic dataset: ", data_file, "\n")
  
  read_and_process_data_res = read_and_process_semi_data(data_file)
  d = read_and_process_data_res$d
  gt_te = read_and_process_data_res$gt_te
  
  # train/test split
  set.seed(1)
  train_indices = sample(1:nrow(d), size = 0.5 * nrow(d), replace = FALSE)
  
  d_train = d[train_indices, ]
  d_test = d[-train_indices, ]
  
  gt_te_test = gt_te[-train_indices]
  
  x = d_train[, -which(names(d_train) %in% c("Y", "treatment"))]
  y = d_train$Y
  z = d_train$treatment
  causal_model = causalDT(X=x, Y=y, Z=z, rpart_prune="min")
  tree_model = causal_model$student_fit$fit
  
  exp_res = collect_exp_res(tree_model, gt_te_test, d_test)
  res_df = rbind(res_df, data.frame(
    dataset = basename(data_file),
    subgroup_treatment_effects = I(list(exp_res$subgroup_treatment_effects)),
    subgroup_gt_treatment_effects_var = I(list(exp_res$subgroup_gt_treatment_effects_var)),
    subgroup_diff_average = I(list(exp_res$subgroup_diff_average)),
    subgroup_diff_individual = I(list(exp_res$subgroup_diff_individual)),
    max_subgroup_treatment_effect = exp_res$max_subgroup_treatment_effect,
    which_max = exp_res$which_max, 
    max_subgroup_gt_treatment_effct = exp_res$max_subgroup_gt_treatment_effct
  ))
  
  subgroup_coverage = sapply(exp_res$subgroup_diff_individual, length)
  res_df_summary = rbind(res_df_summary, data.frame(
    dataset = basename(data_file),
    max_subgroup_treatment_effect = exp_res$max_subgroup_treatment_effect,
    subgroup_diff_average_weighted_mean = weighted.mean(exp_res$subgroup_diff_average, 
                                                        w = subgroup_coverage, 
                                                        na.rm = TRUE),
    subgroup_diff_average_unweighted_mean = mean(exp_res$subgroup_diff_average, na.rm = TRUE),
    max_subgroup_diff_average = 
      exp_res$subgroup_diff_average[exp_res$which_max],
    max_subgroup_gt_treatment_effct = exp_res$max_subgroup_gt_treatment_effct
  ))
  if((i %% 10 == 0) | (i == length(files))){
    write.csv(res_df,
              file = paste0("./semi_distill_results_full.csv"),
              row.names = FALSE,
              quote = TRUE)
    write.csv(res_df_summary,
              file = paste0("./semi_distill_results.csv"),
              row.names = FALSE,
              quote = TRUE)
  }
}




