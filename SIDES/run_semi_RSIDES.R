rm(list=ls())
options(warn = 2)


source("./utils.R")

require(rsides)
col_type <- function(x) {
  u <- unique(x)
  
  if (length(u) == 2) {                # exactly two values → binary
    "nominal"
  } else {                             # everything else
    "numeric"
  }
}
collect_exp_res_RSIDES = function(rsides_results, gt_te_test, dd){
  subgroups = rsides_results$patient_subgroups$Subgroups$subgroup
  treatment_effects = rep(0, length(subgroups))
  gt_subgroup_effects = rep(0, length(subgroups))
  max_te = -Inf
  
  rsides_treatment_effects = rep(0, length(subgroups))
  rsides_treatment_effects_var = rep(0, length(subgroups))
  diff_average = rep(0, length(subgroups)) # the difference between the subgroup te and the average gt_tex
  diff_individual <- vector("list", length(subgroups)) # the individual difference of the subgroup_te and the gt_te
  
  
  
  
  
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
    bool_t1 = bool & (dd$treatment == 1)
    bool_t0 = bool & (dd$treatment == 0)
    if(sum(bool_t1) == 0 | sum(bool_t0) == 0){
      te = NA
      gt_subgroup_effects[i] = NA
    } else{
      te = mean(dd$Y[bool_t1]) - mean(dd$Y[bool_t0])
      
      gt_subgroup_effects[i] = mean(gt_te[bool])
      
      if(max_te < te){
        max_te = te
        max_sg_bool = bool
      }  
    }
    rsides_treatment_effects[i] = te
    rsides_treatment_effects_var[i] = var(gt_te_test[bool])
    
    # get the ground truth subgroup te: 
    gt_te_subgroup = gt_te_test[bool]
    diff_average[i] = te - mean(gt_te_subgroup, na.rm = TRUE)
    diff_individual[[i]] = gt_te_subgroup - te
    
    
  }
  
  which_max = which.max(treatment_effects)
  
  return(list(
    subgroup_treatment_effects = rsides_treatment_effects,
    subgroup_gt_treatment_effects_var = rsides_treatment_effects_var,
    subgroup_diff_average = diff_average,
    subgroup_diff_individual = diff_individual,
    max_subgroup_treatment_effect = rsides_treatment_effects[which_max],
    which_max = which_max,
    max_subgroup_gt_treatment_effct = mean(gt_te_test[max_sg_bool])
  )
  )
  
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
  cat("running RSIDES on semi-synthetic dataset: ", data_file, "\n")
  
  read_and_process_data_res = read_and_process_semi_data(data_file)
  d = read_and_process_data_res$d
  gt_te = read_and_process_data_res$gt_te
  
  # train/test split
  set.seed(1)
  train_indices = sample(1:nrow(d), size = 0.5 * nrow(d), replace = FALSE)
  
  d_train = d[train_indices, ]
  d_test = d[-train_indices, ]
  
  gt_te_test = gt_te[-train_indices]
  
  # Fit a RSIDES model
  if(length(unique(d$y_factual)) == 2){
    output_type = "binary"
    analysis_method = "Z-test for proportions"
    
  } else{
    output_type = "continuous"
    analysis_method = "T-test"
  }
  colnames(d)[ncol(d)] = "Y"
  colnames(d)[ncol(d)-1] = "treamtment"
  
  
  d_feature = d[, -c(ncol(d)-1, ncol(d))]  # exclude T and Y
  col_types = sapply(d_feature, col_type)
  
  class_cols = names(col_types[col_types == "nominal"])
  cont_cols = names(col_types[col_types == "numeric"])
  endpoint_parameters = list(outcome_variable = "Y",
                             type = output_type,
                             label = "Y",
                             analysis_method = analysis_method,
                             cont_covariates = paste(cont_cols, collapse = ", "),
                             class_covariates = paste(class_cols, collapse = ", "),
                             direction = 1)
  
  biomarker_names = colnames(d)[1:(ncol(d)-2)]
  biomarker_types = sapply(d[, biomarker_names], col_type)
  
  data_set_parameters = list(data_set = d_train,
                             treatment_variable_name = "treatment",
                             treatment_variable_control_value = "0",
                             biomarker_names = biomarker_names,
                             biomarker_types = biomarker_types)
  
  algorithm_parameters = list(
    n_perms_mult_adjust = 10,
    min_subgroup_size = 10,
    subgroup_search_algorithm = "SIDES procedure",
    ncores = 1,
    random_seed = 42,
    depth = min(5, ncol(d) - 2),
    width = 2)
  
  parameters = list(endpoint_parameters = endpoint_parameters,
                    data_set_parameters = data_set_parameters,
                    algorithm_parameters = algorithm_parameters)
  results = SubgroupSearch(parameters)
  
  exp_res = collect_exp_res_RSIDES(results, gt_te_test, d_test)
  
  
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
  if((i %% 10 == 0) | (i == length(files))) {
    write.csv(res_df,
              file = paste0("./FINAL_semi_rsides_results_full.csv"),
              row.names = FALSE,
              quote = TRUE)
    write.csv(res_df_summary,
              file = paste0("./FINAL_semi_rsides_results.csv"),
              row.names = FALSE,
              quote = TRUE)
  }
  # write.csv(res_df, 
  #           file = paste0("./semi_rsides_results_full_IHDP.csv"), 
  #           row.names = FALSE, 
  #           quote = TRUE)
  # write.csv(res_df_summary,
  #           file = paste0("./semi_rsides_results_IHDP.csv"), 
  #           row.names = FALSE, 
  #           quote = TRUE)
}
