rm(list = ls())
library(quint)

source("./util_synthetic.R")
# load the datasets
simulator_names = c(
  "simulate1",
  "simulate_imbalance_treatment",
  "simulate_long_rule"
)

get_rules <- function(df) {
  # split the “2,3” text into integer vectors
  kids <- strsplit(df$childnodes, ",")
  names(kids) <- df$parentnode
  
  rules <- list()
  walk_tree <- function(node, path = character()) {
    row <- match(node, df$parentnode)          # where does this node split?
    
    if (is.na(row)) {                          # no split ⇒ we reached a leaf
      rules[[length(rules) + 1]] <<- paste(path, collapse = " & ")
      return(invisible())
    }
    
    v  <- df$splittingvar[row]
    sp <- df$truesplitpoint[row]
    ch <- as.integer(kids[[row]])              # left, right
    
    walk_tree(ch[1], c(path, sprintf("%s < %s", v, sp)))
    walk_tree(ch[2], c(path, sprintf("%s >= %s", v, sp)))
  }
  
  walk_tree(1)                                 # start at the root
  unlist(rules, use.names = FALSE)
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
  diff_te = rep(0, num_exp_res)
)



for(simulator_name in simulator_names){
  # if(simulator_name != "simulate_long_rule") next()
  for(n in ns){
    cat("running QUINT on simulator: ", simulator_name, " with n = ", n, "\n")
    for(iter_ in iters){
      cat("iteration: ", iter_, "\n")
      counter = counter + 1
      set.seed(counter)
      
      file_path = paste0("../datasets/synthetic/", simulator_name, "/n_", n, "_iter_", iter_, ".csv")
      d = read.csv(file_path)
      
      train_indices = sample(1:n, n * 0.5, replace = FALSE)
      d_train = d[train_indices, ]
      d_test = d[-train_indices, ]
      
      # fit the quint model
      quint_model = quint(formula = Y ~ T | .,
                          data=d_train)
      quint_model_pruned <- prune(quint_model)

      if(is.null(quint_model_pruned$si)){
        subgroups = NULL
      } else{
        subgroups = get_rules(quint_model_pruned$si)  
      }
      # get the subgroup analysis results
      ## initialize 
      max_te = -Inf
      max_sg_bool = rep(T, nrow(d_test))
      gt_bool = get_gt_bool(dd=d_test, simulator_name=simulator_name)
      gt_te = mean(d_test$Y[gt_bool & (d_test$T == 1)]) - 
        mean(d_test$Y[gt_bool & (d_test$T == 0)])
      
      treatment_effects = rep(-Inf, length(subgroups))
      if(is.null(subgroups)){
        bool_test = rep(T, nrow(d_test))
        subgroup_te = mean(d_test$Y[bool_test & d_test$T == 1]) - 
          mean(d_test$Y[bool_test & d_test$T == 0])
        treatment_effects = subgroup_te
        if(subgroup_te > max_te){
          max_te = subgroup_te
          max_subgroup_jaccard = sum(bool_test & gt_bool) / 
            sum(bool_test | gt_bool)
          abs_error_te = abs(max_te - gt_te)
          max_sg_bool = bool_test
        }
      } else{
        for(i in 1:length(subgroups)){
          sg = subgroups[i]
          items = strsplit(sg, " & ")[[1]]
          bool_test = rep(T, nrow(d_test))
          for(j in 1:length(items)){
            item = items[j]
            if(grepl(">=", item)){
              feature_value = strsplit(item, " >= ")[[1]]
              operator = ">="
            } else {
              feature_value = strsplit(item, " < ")[[1]]
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
      }
      
      # collect the results
      which_max = which.max(treatment_effects)
      jaccard_similarity = sum(gt_bool & max_sg_bool) / sum(gt_bool | max_sg_bool)
      
      full_results[counter, ] = list(simulator_name, n, iter_, jaccard_similarity, 
                                     treatment_effects[which_max], gt_te,
                                     diff_te = abs(treatment_effects[which_max] - gt_te))
    }
    full_results$jaccard = round(full_results$jaccard, 3)
    full_results$max_te = round(full_results$max_te, 3)
    full_results$max_te = round(full_results$max_te, 3)
    full_results$gt_te = round(full_results$gt_te, 3)
    full_results$diff_te = round(full_results$diff_te, 3)
    
    write.csv(full_results, "res_quint_simulation.csv", row.names = FALSE)
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

write.csv(summary_res, "summary_res_quint_simulations.csv", row.names = FALSE)
