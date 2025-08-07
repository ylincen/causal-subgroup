import pandas as pd  
import numpy as np
from dowhy import CausalModel
import re  
import ast
import warnings  

warnings.filterwarnings('ignore', category=DeprecationWarning) 

def get_Metric( train_df, test_df, model, rule_num = 2, thre = 0.05, outcome_col = None, treatment_col = None, covariates = None, 
               te_col = None, te_pre_col = None, rules = None, cate = True, ate = None): 
    """  
    This function calculates the metrics for the given model.  
      
    Args:  
        train_df (pd.DataFrame): The training dataset used for the model.  
        test_df (pd.DataFrame): The testing dataset used for the model evaluation.  
        model (str): The model used to get rules.  
        rule_num (int, optional): The number of rules to be used. Defaults to 2.  
        thre (float, optional): The threshold of samples covered for the rule. Defaults to 0.05.  
        outcome_col (str, optional): The name of the column in the dataframe that represents the outcome variable.  
        treatment_col (str, optional): The name of the column in the dataframe that represents the treatment variable.  
        covariates (list, optional): The list of covariates to be used in the model.  
        te_col (str, optional): The name of the column in the dataframe that represents the true treatment effect.  
        te_pre_col (str, optional): The name of the column in the dataframe that represents the predicted treatment effect.  
        rules (list, optional): The list of rules get from the model.  
        cate (bool, optional): If True, calculate the Conditional Average Treatment Effect (CATE). Defaults to True.  
        ate (float, optional): Get from the 'CRE'. Defaults to None.  
          
    Returns:  
        dict: two dataframe containing the calculated metrics. 
    """
    
    metrics = pd.DataFrame(index = [model])
    
    if rules is not None and len(rules) != 0:

        #preprocess the rules and get samples that satisfy the rules
        rule_list = preprocess_rules(rules)
        train_met = rule_metrics(train_df , rule_list)
        train_dfs = train_met[3]
    
        #calculate the average treatment effect of the rule group
        if cate is False:
            cate_list = []
            for rdf in train_dfs:
                cate = cal_ate(rdf, treatment_col, outcome_col, covariates)
                cate_list.append(cate)
        elif isinstance(cate, list):
            cate_list = cate
        
        #get the top k rules with the highest effect  
        if model != "CF":
            train_dfs, rule_list, cate_list = topk_rule(rule_list, cate_list, train_df, train_dfs, rule_num, thre)
        else:
            train_dfs, rule_list, cate_list = topk_rule(rule_list, cate_list, train_df, train_dfs, rule_num, thre/2)

        rule_df = pd.DataFrame(rule_list, columns=['rules'])
        rule_df['cate_pre'] = cate_list

        #calculate the coverage of the rule group
        cov_list_train = []
    
        for idx, tdf in enumerate(train_dfs):

            cov = len(tdf)/len(train_df)
            cov_list_train.append(round(cov,3))        

        #evaluation with the test set
        test_met = rule_metrics(test_df , rule_list)
        metrics["rule_count"] = [test_met[0]]
        metrics["avg_length"] = [test_met[1]]      
        metrics["test_overlap"] = [test_met[2]]
        test_dfs = test_met[3]
        
        # Assign the predicted treatment effect to the test samples 
        selected_df = pd.concat(test_dfs)  
        selected_df[te_pre_col] = np.nan  
        temp_dfs = [selected_df.copy() for _ in test_dfs]  

        for i, sub_df in enumerate(test_dfs):  
            sub_df[te_pre_col] = cate_list[i]
            temp_dfs[i][te_pre_col] = temp_dfs[i].index.isin(sub_df.index) * cate_list[i]  

        if model == "CRE":
            selected_df[te_pre_col] = pd.concat(temp_dfs).groupby(level=0)[te_pre_col].sum()    
            counts = selected_df.index.value_counts()  
            selected_df = selected_df.drop_duplicates() 
            selected_df[te_pre_col] = selected_df[te_pre_col] - (counts - 1) * ate   

        else:
            selected_df[te_pre_col] = pd.concat(temp_dfs).groupby(level=0)[te_pre_col].mean() 
            selected_df = selected_df.drop_duplicates()  

        # Assign the predicted treatment effect to the train samples 
        covered_train = pd.concat(train_dfs)  
        covered_train[te_pre_col] = np.nan  
        temp_dfs = [covered_train.copy() for _ in train_dfs]

        for i, sub_df in enumerate(train_dfs):  
             
            sub_df[te_pre_col] = cate_list[i]
            temp_dfs[i][te_pre_col] = temp_dfs[i].index.isin(sub_df.index) * cate_list[i] 

        if model == "CRE":
            covered_train[te_pre_col] = pd.concat(temp_dfs).groupby(level=0)[te_pre_col].sum()    
            counts = covered_train.index.value_counts()  
            covered_train = covered_train.drop_duplicates() 
            covered_train[te_pre_col] = covered_train[te_pre_col] - (counts - 1) * ate  

        else:
            covered_train[te_pre_col] = pd.concat(temp_dfs).groupby(level=0)[te_pre_col].mean() 
            covered_train = covered_train.drop_duplicates()        

        metrics["test_coverage"] = (len(selected_df) / len(test_df))

        cov_list_test = []   
        rule_df["test_avg_ite(no weight)"] = np.nan  

        for idx, rdf in enumerate(test_dfs):  

            #calculate the coverage of the rule group
            cov = len(rdf)/len(test_df)
            cov_list_test.append(round(cov,3))
            
            #calculate the average treatment effect of the rule group
            rule_df.iloc[idx, rule_df.columns.get_loc("test_avg_ite(no weight)")] = (rdf[te_col]).sum() / len(rdf) 

        #get the variance of the rules
        var = cal_var(train_dfs, outcome_col, treatment_col, te_col)

        var = cal_var(test_dfs, outcome_col, treatment_col, te_col)
        rule_df['test_tr_var_y'] = var["tr_var_list"] 
 
        
    else:
        print(f"Warning: rule not provided.")
        
    
    # calculate pehe and mape
    if te_pre_col and te_col is not None:
        
        pehe_list = []
        mape_list = []

        for cdf in train_dfs:
            ite = cdf[te_col]
            te_pre = cdf[te_pre_col]
            pehe = np.sqrt(np.nanmean((ite - te_pre) ** 2))  
            mape = np.nanmean(np.abs((ite - te_pre) / ite))
            pehe_list.append(pehe)
            mape_list.append(mape)

        pehe_list = []
        mape_list = []

        for tdf in test_dfs:
            ite = tdf[te_col]
            te_pre = tdf[te_pre_col]
            pehe = np.sqrt(np.nanmean((ite - te_pre) ** 2))  
            mape = np.nanmean(np.abs((ite - te_pre) / ite))
            pehe_list.append(pehe)
            mape_list.append(mape)

        rule_df["test_PEHE"] = pehe_list
        rule_df["test_MAPE"] = mape_list

    rule_df.index = [model + " rule" + str(i) for i in range(len(rule_df))]

    return metrics, rule_df


#calculate the variance of a list of rule groups
def cal_var( test_dfs, outcome_col, treatment_col, te_col):
    tr_var_list = []  
    con_var_list = []  
    te_var_list = []
    for df in test_dfs:  
        result = cal_var_df(df, outcome_col, treatment_col, te_col)  
        tr_var_list.append(result["tr_var"])  
        con_var_list.append(result["con_var"])      
        te_var_list.append(result["te_var"])  

    result = {
        "tr_var_list": tr_var_list,
        "con_var_list": con_var_list,
        "te_var_list": te_var_list
    }

    return result

#calculate the variance of one rule group 
def cal_var_df(df, outcome_col, treatment_col, te_col):  
    
    temp1 = np.sum(df[outcome_col] * df['wt'] * df[treatment_col])  
    temp0 = np.sum(df[outcome_col] * df['wt'] * (1 - df[treatment_col]))  
    twt = np.sum(df['wt'])  
    ttreat = np.sum(df['wt'] * df[treatment_col])  
    tr_sqr_sum = np.sum(df[outcome_col]**2 * df['wt'] * df[treatment_col])  
    con_sqr_sum = np.sum(df[outcome_col]**2 * df['wt'] * (1 - df[treatment_col]))  
  
    tr_var = tr_sqr_sum / ttreat - (temp1 / ttreat)**2  
    if(twt - ttreat == 0):
        print("warning: weight of control group:" + str(temp0))
    con_var = con_sqr_sum / (twt - ttreat) - (temp0 / (twt - ttreat))**2  
    
    tempTE = np.sum(df[te_col] * df['wt'])  
    twtTE = np.sum(df['wt'])  
    
    te_sqr_sum = np.sum(df[te_col]**2 * df['wt'])  
    te_var = te_sqr_sum / twtTE - (tempTE / twtTE)**2  
    
    result = {   
        "tr_var": tr_var,  
        "con_var": con_var,
        "te_var": te_var
    }  
  
    return result  


#get the top k rules with the highest effect
def topk_rule( rule_list, cate, train_df, train_dfs, k, thre):
    
    #order the rules by effect
    zipped = list(zip(cate, rule_list, train_dfs))       
    zipped_sorted = sorted(zipped, key=lambda x: x[0], reverse=True)  
    cate, rule_list, train_dfs = zip(*zipped_sorted) 
    cate = list(cate)  
    rule_list = list(rule_list)  
    train_dfs = list(train_dfs)  

    for i in reversed(range(len(train_dfs))):  
        if len(train_dfs[i]) < len(train_df) * thre: 
            del train_dfs[i]  
            del rule_list[i]
            del cate[i]

    if len(rule_list) < k:  
        print(f"Warning: There are only {len(rule_list)} in the rules, less than {k}.")
        return train_dfs, rule_list, cate
    
    train_dfs_k = train_dfs[:k]  
    rule_list_k = rule_list[:k]
    cate_k = cate[:k]

    return train_dfs_k, rule_list_k, cate_k
    
  
def preprocess_rules(raw_rules):  
    if isinstance(raw_rules, list):
        rules_list = [rule for rule in raw_rules if 't' not in rule]  
        return rules_list
        
    # Remove brackets, 'if\n' and parentheses  
    rules = raw_rules.replace("[", "").replace("]", "").replace("if\n", "").replace("(", "").replace(")", "")  
      
    # Replace connectors with 'and' and 'or'  
    rules = rules.replace("^", " and ").replace(" v\n", " or ")
      
    # Split rules into a list  
    rules_list = rules.split(" or ")  
      
    # Remove the outcome from the last rule  
    rules_list = [rule.split("\nthen\n")[0].strip() for rule in rules_list if rule.strip()]  

    # Remove the rules with 't'
    rules_list = [rule for rule in rules_list if 't' not in rule]  

    return rules_list  
    
def rule_metrics(dataset, rules = None):  
    if rules is None:  
        print("No rules provided.")  
        return  
    
    total_length = 0  
    overlap = 0
    rule_count = len(rules)
    
    #the samlpe list
    rule_dfs = []
    for rule in rules:            
        # Assuming the format is "condition"  
        conditions = re.split(r' and | or ', rule)  # Split by 'and' or 'or'  
        
        # Calculate rule length  
        rule_length = len(conditions)  
        total_length += rule_length  
  
        #get samples that satisfy the rule
        rule_df = dataset.copy()  
        for condition in conditions:  
            # Handle different types of conditions  
            if "==" in condition:  
                column, value = condition.split("==")  
                value = value.strip()  
                if value.isnumeric(): 
                    rule_df.loc[:, column.strip()] = rule_df.loc[:, column.strip()].astype(float)   
                    rule_df = rule_df[rule_df[column.strip()] == float(value)]  
                else: 

                    try:  
                        value = ast.literal_eval(value)  

                    except ValueError:   
                        pass  
                    rule_df = rule_df[rule_df[column.strip()] == value]  

            elif "!=" in condition:  
                column, value = condition.split("!=")  
                value = value.strip()  

                if value.isnumeric():  
                    rule_df.loc[:, column.strip()] = rule_df.loc[:, column.strip()].astype(float)  
                    rule_df = rule_df[rule_df[column.strip()] != float(value)]  

                else:  
                    try:   
                        value = ast.literal_eval(value)  

                    except ValueError:   
                        pass 

                    rule_df = rule_df[rule_df[column.strip()] != value]  
            elif "<=" in condition:  
                column, value = condition.split("<=")  
                rule_df.loc[:, column.strip()] = rule_df.loc[:, column.strip()].astype(float)  
                rule_df = rule_df[rule_df[column.strip()] <= float(value)]  
            elif ">=" in condition:  
                column, value = condition.split(">=")  
                rule_df.loc[:, column.strip()] = rule_df.loc[:, column.strip()].astype(float)  
                rule_df = rule_df[rule_df[column.strip()] >= float(value)]  
            elif "<" in condition:  
                column, value = condition.split("<")  
                rule_df.loc[:, column.strip()] = rule_df.loc[:, column.strip()].astype(float)  
                rule_df = rule_df[rule_df[column.strip()] < float(value)]  
            elif ">" in condition:  
                column, value = condition.split(">")  
                rule_df.loc[:, column.strip()] = rule_df.loc[:, column.strip()].astype(float)  
                rule_df = rule_df[rule_df[column.strip()] > float(value)]  
            elif " in " in condition:  
                column, value = condition.split(" in ")  
                value_list = ast.literal_eval(value.strip())  
                rule_df = rule_df[rule_df[column.strip()].isin(value_list)] 
            else:  
                raise ValueError(f"Unknown condition {condition} in rule {rule}")  
           
        rule_dfs.append(rule_df)

    for i in range(rule_count):  
        for j in range(i+1, rule_count):  
            overlap += len(pd.merge(rule_dfs[i], rule_dfs[j], how='inner')) 
    if rule_count not in [0, 1]:
        overlap = overlap/len(dataset)/((rule_count-1)*rule_count/2)
    else:
        overlap = overlap/len(dataset)
    
    if rule_count == 0:
        print(str(rules)) 
    avg_length = total_length / rule_count    
  
    return rule_count, avg_length,overlap,rule_dfs

def cal_ate(df, treatment, outcome, covariates):

    """
    if df[treatment].nunique() == 1:
        return np.nan
    model= CausalModel(data = df, treatment = [treatment], outcome = outcome, common_causes = covariates)    
    identified_estimand = model.identify_effect(proceed_when_unidentifiable=True)
    causal_estimate_ipw = model.estimate_effect(identified_estimand,
                                            method_name="backdoor.propensity_score_weighting",
                                            target_units = "ate",
                                            method_params={"weighting_scheme":"ips_normalized_weight"})
    return causal_estimate_ipw.value
    """
    
    twt = np.sum(df['wt'])  
    ttreat = np.sum(df['wt'] * df[treatment])  
    temp1 = np.sum(df[outcome] * df['wt'] * df[treatment])
    temp0 = np.sum(df[outcome] * df['wt'] * (1 - df[treatment]))

    if ttreat == 0 or (twt - ttreat) == 0:
        return np.nan
    effect = (temp1 / ttreat) - (temp0 / (twt - ttreat))

    return effect
