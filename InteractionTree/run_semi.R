rm(list=ls())
library(dplyr)
source("./Functions-IT.R")
source("./utils.R")

files = list.files("../datasets/semi_synthetic/", pattern="*.csv", full.names=T)

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

collect_exp_res_IT = function(leaf_rules, gt_te_test, d_test){
  subgroups = leaf_rules
  treatment_effects = rep(0, length(subgroups))
  gt_subgroup_effects = rep(0, length(subgroups))
  max_te = -Inf
  subgroup_gt_treatment_effects_var = rep(0, length(subgroups))
  
  
  diff_average = rep(0, length(subgroups)) # the difference between the subgroup te and the average gt_tex
  diff_individual <- vector("list", length(subgroups)) 
  max_sg_bool = rep(T, nrow(d_test)) # the boolean vector for the max subgroup, used to calculate the max subgroup treatment effect
  
  if(length(leaf_rules) == 0){
    bool_test = rep(T, nrow(d_test))
    
    subgroup_te = mean(d_test$y[bool_test & d_test$trt == 1]) - 
      mean(d_test$y[bool_test & d_test$trt == 0])
    treatment_effects = subgroup_te
    subgroup_gt_treatment_effects_var = var(gt_te_test[bool_test])
    diff_average = subgroup_te - mean(gt_te_test[bool_test])
    diff_individual[[1]] = subgroup_te - gt_te_test[bool_test] 
    if(subgroup_te > max_te){
      max_te = subgroup_te
      abs_error_te = abs(max_te - gt_te)
      max_sg_bool = bool_test
    }
    
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
      if((sum(d_test$y[bool_test & d_test$trt == 1])==0) | 
         (sum(d_test$y[bool_test & d_test$trt == 0])==0)){
        treatment_effects[i] = -Inf
        subgroup_gt_treatment_effects_var[i] = NA
        diff_average[i] = NA
        diff_individual[[i]] = NA
      } else{
        treatment_effects[i] = subgroup_te
        subgroup_gt_treatment_effects_var[i] = var(gt_te_test[bool_test])
        diff_average[i] = subgroup_te - mean(gt_te_test[bool_test])
        diff_individual[[i]] = subgroup_te - gt_te_test[bool_test] 
        if(subgroup_te > max_te){
          max_te = subgroup_te
          abs_error_te = abs(max_te - gt_te)
          max_sg_bool = bool_test
        } 
      }
      
    }
  }
  which_max=which.max(treatment_effects)

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

  cat("running IT on semi-synthetic dataset: ", data_file, "\n")
  read_and_process_data_res = read_and_process_semi_data(data_file)
  d = read_and_process_data_res$d
  gt_te = read_and_process_data_res$gt_te
  
  n = nrow(d)
  colnames(d)[ncol(d)] = "y"
  colnames(d)[ncol(d)-1] = "trt"  # ensure the last column is Y
  
  train_indices = sample(1:n, n * 0.5, replace = FALSE)
  d_train = d[train_indices, ]
  d_test = d[-train_indices, ]
  gt_te_test = gt_te[-train_indices]
  
  val_indices = sample(1:nrow(d_train), nrow(d_train) * 0.2, replace = FALSE)
  d_val = d_train[val_indices, ]
  d_train = d_train[-val_indices, ]
  
  dat = d_train
  test.dat = d_val
  

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
  
  exp_res = collect_exp_res_IT(leaf_rules, 
                              gt_te_test = gt_te_test, 
                              d_test = d_test
                              )
  # collect the results
  which_max = which.max(exp_res$treatment_effects)
  
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
  if((i %% 10 == 0)|(i == length(files))){
    write.csv(res_df,
              file = paste0("./NEW_semi_IT_results_full.csv"),
              row.names = FALSE,
              quote = TRUE)
    write.csv(res_df_summary,
              file = paste0("./NEW_semi_IT_results.csv"),
              row.names = FALSE,
              quote = TRUE)
  }
}

