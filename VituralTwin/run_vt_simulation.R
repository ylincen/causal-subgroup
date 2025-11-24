rm(list=ls())
require(aVirtualTwins)

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
    cat("running VT on simulator: ", simulator_name, " with n = ", n, "\n")
    for(iter_ in iters){
      counter = counter + 1
      set.seed(counter)
      
      file_path = paste0("../datasets/synthetic/", simulator_name, "/n_", n, "_iter_", iter_, ".csv")
      d = read.csv(file_path)
      
      train_indices = sample(1:n, n * 0.5, replace = FALSE)
      d_train = d[train_indices, ]
      d_test = d[-train_indices, ]
      
      tree_type = ifelse(length(unique(d$Y)) > 2, "regression", "class")
      vt.obj <- vt.data(dataset         = d_train,
                        outcome.field   = "Y",
                        treatment.field = "T",
                        interactions    = TRUE)
      # Run the VT method
      ## First step : create random forest model
      vt.for <- vt.forest(forest.type  = "one",
                          vt.data      = vt.obj,
                          interactions = TRUE,
                          ntree        = 500)
      ## Second step : find rules in data 
      vt.trees <- vt.tree(tree.type = tree_type,
                          vt.difft  = vt.for, 
                          threshold = quantile(vt.for$difft, seq(.5,.8,.1))
                          )
      ## Print results
      vt.sbgrps <- vt.subgroups(vt.trees)
      
      # get the subgroup analysis results
      ## initialize 
      max_te = -Inf
      max_sg_bool = rep(T, nrow(d_test))
      gt_bool = get_gt_bool(dd=d_test, simulator_name=simulator_name)
      gt_te = mean(d_test$Y[gt_bool & (d_test$T == 1)]) - 
        mean(d_test$Y[gt_bool & (d_test$T == 0)])
      
      treatment_effects = rep(-Inf, length(vt.sbgrps$Subgroup))
      for(i in 1:length(vt.sbgrps$Subgroup)){
        sg = vt.sbgrps$Subgroup[i]
        items = strsplit(sg, " & ")[[1]]
        bool_test = rep(T, nrow(d_test))
        for(j in 1:length(items)){
          item = items[j]
          if(grepl(">", item)){
            feature_value = strsplit(item, ">=")[[1]]
            operator = ">="
          } else {
            feature_value = strsplit(item, "< ")[[1]]
            operator = "<"
          }
          feature = feature_value[1]
          value = as.numeric(feature_value[2])
          if(operator == ">="){
            bool_test = bool_test & (d_test[[feature]] >= value)
          } else {
            bool_test = bool_test & (d_test[[feature]] < value)
          }
        }
        subgroup_te = mean(d_test$Y[bool_test & d_test$T == 1]) - 
          mean(d_test$Y[bool_test & d_test$T == 0])
        treatment_effects[i] = subgroup_te
        if(subgroup_te > max_te){
          max_te = subgroup_te
          max_subgroup_jaccard = sum(bool_test & gt_bool) / 
            sum(bool_test | gt_bool)
          abs_error_te = abs(max_te - gt_te)
          max_sg_bool = bool_test
        }
      }
      # collect the results
      which_max = which.max(treatment_effects)
      jaccard_similarity = sum(gt_bool & max_sg_bool) / sum(gt_bool | max_sg_bool)

      gt_te_theoretical_per_sample = gt_te_per_sample(d_test, simulator_name)
      gt_te_of_gt_subgroup_ = mean(gt_te_theoretical_per_sample[gt_bool])
      gt_te_of_learned_subgroup_ = mean(gt_te_theoretical_per_sample[max_sg_bool])
      estimated_te_of_learned_subgroup_ = mean(d_test$Y[max_sg_bool & (d_test$T == 1)]) - 
        mean(d_test$Y[max_sg_bool & (d_test$T == 0)])
      estimated_te_of_gt_subgroup_ = mean(d_test$Y[gt_bool & (d_test$T == 1)]) -
        mean(d_test$Y[gt_bool & (d_test$T == 0)])
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
    write.csv(full_results, "res_vt_simulation.csv", row.names = FALSE)
    
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

write.csv(summary_res, "summary_res_VT_simulations.csv", row.names = FALSE)
