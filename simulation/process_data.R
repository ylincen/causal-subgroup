# Most, if not all, datasets come from https://github.com/rguo12/awesome-causality-data
rm(list = ls())
require(dplyr)
require(data.table)
require(ggplot2)




process_acic2016 = function(){
  # if (require("remotes", quietly = TRUE) == FALSE) {
  #   install.packages("remotes")
  #   require("remotes")
  # }
  # remotes::install_github("vdorie/aciccomp/2016")
  require(aciccomp2016)
  for(i in 1:nrow(parameters_2016)){
    cat("Processing ACIC2016 dataset ", i, " of ", nrow(parameters_2016), "\n")
    acic2016 = dgp_2016(input_2016, parameters_2016[i,], random.seed = 1)
    dd2 = input_2016
    # one-hot encode dd$x2 and dd$x_24
    library(fastDummies)
    
    dd2$x_2 = as.factor(dd2$x_2)
    dd2$x_24 = as.factor(dd2$x_24)
    dd <- dummy_cols(dd2, select_columns = c("x_2", "x_21" ,"x_24"), remove_first_dummy = F)
    # remove the original columns as their one-hot encoded version have been included.
    dd$x_2 = NULL
    dd$x_24 = NULL
    dd$x_21 = NULL
    
    dd$treatment = acic2016$z
    dd$y0 = acic2016$y.0
    dd$y1 = acic2016$y.1
    
    dd$y_factual = ifelse(dd$treatment == 1, dd$y1, dd$y0)
    dd$y_cfactual = ifelse(dd$treatment == 1, dd$y0, dd$y1)
    
    dd$y0 = NULL
    dd$y1 = NULL
    
    write.csv(dd, file = paste0("clean_data/ACIC2016_", i, ".csv"), 
              row.names = FALSE, quote = FALSE)
  }
}


# run the simulation
process_acic2016()


