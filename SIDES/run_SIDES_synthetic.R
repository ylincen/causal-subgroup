rm(list = ls())
library(rpart)
library(rpart.plot)

set.seed(1)

source("./util_synthetic.R")

require(rsides)
col_type <- function(x) {
  u <- unique(x)
  
  if (length(u) == 2) {                # exactly two values → binary
    "nominal"
  } else {                             # everything else
    "numeric"
  }
}


get_subgroup_bool = function(dd, rsides_results, gt_max_bool){
  # dd can be d_train or d_test
  # rsides_results is the output of SubgroupSearch of the RSIDES package
  # gt_max_bool is a boolean vector indicating the rows of dd that are in the ground_truth subgroup
  subgroups = rsides_results$patient_subgroups$Subgroups$subgroup
  treatment_effects = rep(0, length(subgroups))
  max_te = -Inf
  for(i in 1:length(subgroups)){
    sg = subgroups[[i]]
    if(sg == "Overall population"){
      bool = rep(TRUE, nrow(dd))
    } else {
      bool = rep(TRUE, nrow(dd))
      rules = strsplit(sg, ":")[[1]]
      for(j in 1:length(rules)){
        rule = rules[j]
        # get the operator from ">", "<", ">=", "<=", "="
        operator <- regmatches(rule, regexpr("(>=|<=|>|<|=)", rule))
        
        feature_value = strsplit(rule, operator)[[1]]
        feature = feature_value[1]
        value = as.numeric(feature_value[2])
        bool = bool & switch(
          operator,
          ">" = dd[[feature]] > value,
          "<" = dd[[feature]] < value,
          ">=" = dd[[feature]] >= value,
          "<=" = dd[[feature]] <= value,
          "=" = dd[[feature]] == value
        )
      }
    }
    # get the treatment effects
    bool_t1 = bool & (dd$T == 1)
    bool_t0 = bool & (dd$T == 0)
    if(sum(bool_t1) == 0 | sum(bool_t0) == 0){
      te = -Inf
    } else{
      te = mean(dd$Y[bool_t1]) - mean(dd$Y[bool_t0])
      
      
      if(max_te < te){
        max_te = te
        max_sg_bool = bool
      }  
    }
    
    treatment_effects[i] = te
  }
  
  which_max = which.max(treatment_effects)
  # jaccard for that subgroup
  jaccard_similarity = sum(gt_max_bool & max_sg_bool) / sum(gt_max_bool | max_sg_bool)
  return(list(
    jaccard_similarity = jaccard_similarity,
    treatment_effect = treatment_effects[which_max],
    max_sg_bool = max_sg_bool
  )
  )
}

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
rsides_jaccard = rep(0, num_exp_res)
rsides_max_te = rep(0, num_exp_res)
counter = 0

full_results = data.frame(
  simulator_name = rep("", num_exp_res),
  n = rep(0, num_exp_res),
  iter = rep(0, num_exp_res),
  rsides_jaccard = rsides_jaccard,
  rsides_max_te = rsides_max_te,
  gt_te = gt_te,
  rsides_diff_te = rep(0, num_exp_res),
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
      
      
      analysis_method = "Z-test for proportions"
      
      d_feature = d[, -c(ncol(d)-1, ncol(d))]  # exclude T and Y
      col_types = sapply(d_feature, col_type)
      
      class_cols = names(col_types[col_types == "nominal"])
      cont_cols = names(col_types[col_types == "numeric"])
      endpoint_parameters = list(outcome_variable = "Y",
                                 type = "binary",
                                 label = "Y",
                                 analysis_method = analysis_method,
                                 cont_covariates = paste(cont_cols, collapse = ", "),
                                 class_covariates = paste(class_cols, collapse = ", "),
                                 direction = 1)
      
      biomarker_names = colnames(d)[1:(ncol(d)-2)]
      biomarker_types = sapply(d[, biomarker_names], col_type)
      
      data_set_parameters = list(data_set = d_train,
                                 treatment_variable_name = "T",
                                 treatment_variable_control_value = "0",
                                 biomarker_names = biomarker_names,
                                 biomarker_types = biomarker_types)
      
      algorithm_parameters = list(
        n_perms_mult_adjust = 10,
        min_subgroup_size = 10,
        subgroup_search_algorithm = "SIDES procedure",
        ncores = 1,
        random_seed = 42,
        depth = ncol(d)-2, 
        width = 2)
      
      parameters = list(endpoint_parameters = endpoint_parameters,
                        data_set_parameters = data_set_parameters,
                        algorithm_parameters = algorithm_parameters)
      results = SubgroupSearch(parameters)
      # results
      
      # run the code
      gt_max_bool = get_gt_bool(d_test, simulator_name)
      res_list = get_subgroup_bool(d_test, results, gt_max_bool)
      gt_te[counter] =
        mean(d_test$Y[gt_max_bool & (d_test$T == 1)]) - mean(d_test$Y[gt_max_bool & (d_test$T == 0)])
  
      
            gt_te_theoretical_per_sample = gt_te_per_sample(d_test, simulator_name)
      gt_te_of_gt_subgroup_ = mean(gt_te_theoretical_per_sample[gt_max_bool])
      gt_te_of_learned_subgroup_ = mean(gt_te_theoretical_per_sample[res_list$max_sg_bool])
      estimated_te_of_gt_subgroup_ =
        mean(d_test$Y[gt_max_bool & (d_test$T == 1)]) - mean(d_test$Y[gt_max_bool & (d_test$T == 0)])
      estimated_te_of_learned_subgroup_ = 
        mean(d_test$Y[res_list$max_sg_bool & (d_test$T == 1)]) - mean(d_test$Y[res_list$max_sg_bool & (d_test$T == 0)])

      diff_gt_te = abs(gt_te_of_gt_subgroup_ - gt_te_of_learned_subgroup_)
      diff_estimated_te = 
        abs(estimated_te_of_gt_subgroup_ - estimated_te_of_learned_subgroup_)
      
      full_results[counter, ] = list(simulator_name, n, iter_, res_list$jaccard_similarity, 
                                     res_list$treatment_effect, gt_te[counter],
                                     rsides_diff_te = abs(res_list$treatment_effect - gt_te[counter]),
                                     gt_te_of_gt_subgroup = gt_te_of_gt_subgroup_,
                                     estimated_te_of_gt_subgroup = estimated_te_of_gt_subgroup_,
                                     gt_te_of_learned_subgroup = gt_te_of_learned_subgroup_,
                                     estimated_te_of_learned_subgroup = estimated_te_of_learned_subgroup_
      )
      print(full_results[counter, ,drop=F])
    }
    
    write.csv(full_results, "res_SIDES.csv", row.names = FALSE)
  }
}
require(data.table)
full_res_dt = data.table(full_results)
# summarize as mean and sd with different n
summary_res = full_res_dt[, .(
  mean_rsides_jaccard = mean(rsides_jaccard),
  mean_rsides_max_te = mean(rsides_max_te),
  mean_gt_te = mean(gt_te),
  mean_rsides_diff_te = mean(rsides_diff_te),
), by = .(n, simulator_name)]

summary_res

write.csv(summary_res, "summary_res_SIDES.csv", row.names = FALSE)
