import test    
import os  
import pandas as pd    
from datetime import datetime      
import warnings    
import sys   
warnings.filterwarnings("ignore")    
  

# run with python evaluation.py config.txt
def main(config_file):
    """
    config_file: str,including the parameters separated by comma in the following order:  :
    - Model for causal inference (CURLS, CRE, CT, CF, RIPPER, BRCG, PYS, DT)  
    - Name of the data file (name.csv) 
    - Covariates for the model (separated by space) (auto for all covariates) 
    - Treatment variable  
    - Outcome variable  
    - need kfold or not (True or False)
    - Additional parameters (as key=value pairs separated by space): max_rule_length,num_rules_to_find,rules_least_cov_ratio,pru_coef,random_state,etc.
    """

    # Set the display options  
    pd.set_option('display.max_colwidth', None)   
    
    # Set the path to the data  
    data_dir = "../data/"  
  
    # Create a directory to save the results  
    rule_list = []    
    metrics_list = []  
  
    # Read the configuration file  
    with open(config_file, 'r') as file:  
        lines = file.readlines()  
  
    # The first line of the config file is the file name  
    filename = lines[0].strip()  
    lines = lines[1:] 
    
    #run the evaluation for each line in the file
    for params in lines:  
        params = params.strip().split(",")    
        model = params[0].strip()    
        data_file = params[1].strip()    
        data_path = os.path.join(data_dir, data_file)    
        covariates = params[2].strip().split()    
        treatment = params[3].strip()    
        outcome = params[4].strip()    
        kfold = params[5].strip()
  
        if kfold.lower() == "true":  
            kfold = True  
        else:  
            kfold = False  
            
        # Handle additional parameters    
        additional_params = {}    
        if len(params) > 6:    
            for param in params[6].strip().split():    
                key, value = param.split('=')    
                try:    
                    additional_params[key] = int(value)    
                except ValueError:    
                    additional_params[key] = float(value)     
  
        #get the column names
        if covariates == ["auto"]:
            covariates = list(pd.read_csv(data_path).columns)
            covariates = [covariate for covariate in covariates if covariate not in ['e', 'wt', 'v', 'TE']]
            covariates.remove(treatment)
            covariates.remove(outcome)
  
        result = test.causal_inference(model=model, data_path=data_path, covariates=covariates, treatment=treatment, outcome=outcome, kfold = kfold, **additional_params)    
        rules = result['rules']
        metrics = result['metrics']
        
        print(rules.to_string())  
        print(metrics.to_string())  
          
        rule_list.append(rules)  
        metrics_list.append(metrics)  
  
    # Save the results to a file
    rules = pd.concat(rule_list)  
    metrics = pd.concat(metrics_list)  
    
    rules.to_csv(f"./result/rule_metric_{filename}.csv", index=True)   

    empty_df = pd.DataFrame([''])  
    empty_df.to_csv(f"./result/rule_metric_{filename}.csv", mode='a', index=False, header=False)  
    
    metrics.to_csv(f"./result/rule_metric_{filename}.csv", mode='a', index=True)  

  
if __name__ == "__main__":    
    # Get the config file name from the command line argument  
    config_file = sys.argv[1]  
    main(config_file)  