rm(list=ls())
library(quint)
source("./utils.R")

collect_exp_res_quint = function(subgroups, gt_te_test, d_test){
  treatment_effects = rep(0, length(subgroups))
  gt_subgroup_effects = rep(0, length(subgroups))
  subgroup_gt_treatment_effects_var = rep(0, length(subgroups))
  
  quint_treatment_effects = rep(0, length(subgroups))
  quint_treatment_effects_var = rep(0, length(subgroups))
  diff_average = rep(0, length(subgroups)) # the difference between the subgroup te and the average gt_tex
  diff_individual <- vector("list", length(subgroups)) # the individual difference of the subgroup_te and the gt_te
  
  
  
  
  max_te = -Inf
  max_sg_bool = rep(T, nrow(d_test))
  treatment_effects = rep(-Inf, length(subgroups))
  
  if(is.null(subgroups)){
    bool_test = rep(T, nrow(d_test))
    subgroup_te = mean(d_test$Y[bool_test & d_test$treatment == 1]) - 
      mean(d_test$Y[bool_test & d_test$treatment == 0])
    treatment_effects = subgroup_te
    subgroup_gt_treatment_effects_var = var(gt_te_test[bool_test])
    diff_average = subgroup_te - mean(gt_te_test[bool_test])
    diff_individual[[1]] = subgroup_te - gt_te_test[bool_test] 
    if(subgroup_te > max_te){
      max_te = subgroup_te
      abs_error_te = abs(max_te - gt_te_test)
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
      subgroup_te = mean(d_test$Y[bool_test & d_test$treatment == 1]) - 
        mean(d_test$Y[bool_test & d_test$treatment == 0])
      if((sum(d_test$y[bool_test & d_test$treatment == 1])==0) | 
         (sum(d_test$y[bool_test & d_test$treatment == 0])==0)){
        treatment_effects[i] = -Inf
        subgroup_gt_treatment_effects_var[i] = NA
        diff_average[i] = NA
        diff_individual[[i]] = NA
      } else{
        treatment_effects[i] = subgroup_te
        subgroup_gt_treatment_effects_var[i] = var(gt_te_test[bool_test])
        diff_average[i] = subgroup_te - mean(gt_te_test[bool_test])
        diff_individual[[i]] = subgroup_te - gt_te_test[bool_test] 
      }
      
      
      # treatment_effects[i] = subgroup_te
      # subgroup_gt_treatment_effects_var[i] = var(gt_te_test[bool_test])
      # diff_average[i] = te - mean(gt_te_test[bool_test])
      # diff_individual[[i]] = te - gt_te_test[bool_test] 
      
      
      if(subgroup_te > max_te){
        max_te = subgroup_te
        abs_error_te = abs(max_te - gt_te_test)
        max_sg_bool = bool_test
      }
    }
  }
  
  # collect the results
  which_max = which.max(treatment_effects)
  
  return(list(
    subgroup_treatment_effects = treatment_effects,
    subgroup_gt_treatment_effects_var = subgroup_gt_treatment_effects_var,
    subgroup_diff_average = diff_average,
    subgroup_diff_individual = diff_individual,
    max_subgroup_treatment_effect = treatment_effects[which_max],
    which_max = which_max,
    max_subgroup_gt_treatment_effct = mean(gt_te_test[max_sg_bool])
  ))
}

fit_quint_model = function(d_train){
  quint_model = quint(formula = Y ~ treatment | .,
                      data=d_train)
  quint_model_pruned <- prune(quint_model)
  # sink()
  
  if(is.null(quint_model_pruned$si)){
    subgroups = NULL
  } else{
    subgroups = get_rules(quint_model_pruned$si)  
  }
  return(subgroups)
}


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

  cat("running quint on semi-synthetic dataset: ", data_file, "\n")
  
  read_and_process_data_res = read_and_process_semi_data(data_file)
  d = read_and_process_data_res$d
  gt_te = read_and_process_data_res$gt_te
  
  # train/test split
  set.seed(1)
  train_indices = sample(1:nrow(d), size = 0.5 * nrow(d), replace = FALSE)
  
  d_train = d[train_indices, ]
  d_test = d[-train_indices, ]
  gt_te_test = gt_te[-train_indices]
  
  # fit a Quint model
  subgroups = fit_quint_model(d_train)

  exp_res = collect_exp_res_quint(subgroups, gt_te_test, d_test)

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
  if(i %% 10 == 0){
    write.csv(res_df,
              file = paste0("./semi_quint_results_full.csv"),
              row.names = FALSE,
              quote = TRUE)
    write.csv(res_df_summary,
              file = paste0("./semi_quint_results.csv"),
              row.names = FALSE,
              quote = TRUE)
  }
  write.csv(res_df,
            file = paste0("./semi_quint_results_full_new.csv"),
            row.names = FALSE,
            quote = TRUE)
  write.csv(res_df_summary,
            file = paste0("./semi_quint_results_new.csv"),
            row.names = FALSE,
            quote = TRUE)
}

# require(data.table)
# exp_df = fread("./semi_quint_results.csv", 
#                stringsAsFactors = FALSE, 
#                na.strings = c("", "NA", "NULL"))
# round up the numbers



