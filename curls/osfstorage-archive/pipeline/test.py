import pandas as pd 
import numpy as np
from sklearn.model_selection import train_test_split
import metric
import baseline 
from sklearn.model_selection import KFold  

#the main function to get the rules and metrics
def causal_inference(model, data_path, covariates, treatment, outcome,kfold = False, **kwargs):  
    
    # process the data
    df = pd.read_csv(data_path)   
    random_state = kwargs.pop('random_state', 42)  
    train_df, test_df = train_test_split(df, test_size = 0.2, random_state = random_state)  
    
    k = kwargs.pop('num_rules_to_find', 2)
    thre = kwargs.pop('rules_least_cov_ratio', 0.05) 
  
    kf = KFold(n_splits=5, shuffle=True, random_state=random_state) 

    #if cross validation 
    if kfold == True:
        metrics_results = []  
        rules_results = []  

        for train_index, test_index in kf.split(df):  
            train_df, test_df = df.iloc[train_index], df.iloc[test_index]   

            if model == "CT":
                import rpy2.robjects as ro  
                from rpy2.robjects import pandas2ri  
                from rpy2.robjects.conversion import localconverter         
                with localconverter(ro.default_converter + pandas2ri.converter):  
                    train_rdf = pandas2ri.py2rpy(train_df)
                    test_rdf = pandas2ri.py2rpy(test_df)
                    
                ro.r('''Sys.setlocale("LC_ALL", "en_US.UTF-8")''')    
                # load the r file  
                ro.r['source']('r_baseline.R')  

                # load and exec r function  
                r_function = ro.globalenv['causal_tree']   
                result = r_function(train_rdf, test_rdf, covariates, treatment, outcome, k, thre,**kwargs)
                
                # process the error
                if result[0][0] == "rules not found enough":
                    print("rules not found enough")
                    continue

                # change dataframe to python 
                with localconverter(ro.default_converter + pandas2ri.converter):  
                    rules = ro.conversion.rpy2py(result[0])

                # change rules to python grammar
                import re    
                rules = rules.reset_index(drop=True)
                for i in range(len(rules)):   
                    rules.loc[i, "rule"] = re.sub(r'%in% c\((.*?)\)', r'in [\1]', rules.loc[i, "rule"])   
                rule_list = list(rules["rule"])
                rule_list = [rule.replace("&", "and") for rule in rule_list]
                        
                cate_list = rules["treatment_effect"].tolist() 

                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                    treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rule_list, cate = cate_list)
                metrics["train_time"] = result[1]
                
            if model == "RIPPER":
                rules, train_time  = baseline.ripper(train_df , covariates, treatment, outcome) 
                if len(rules) == 0:
                    print("no rules found")
                    continue
                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                        treatment_col = treatment,covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = str(rules), cate = False)
                metrics["train_time"] = train_time

            if model == "BRCG":
                rules, train_time = baseline.brcg(train_df , covariates, treatment, outcome)
                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                    treatment_col = treatment,covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = str(rules), cate = False)
                metrics["train_time"] = train_time
                
            if model == "PYS":
                
                rules, train_time = baseline.pys(train_df , covariates, treatment, outcome)
                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                        treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rules, cate = False)
                metrics["train_time"] = train_time

            if model == "CRE":
                import rpy2.robjects as ro  
                from rpy2.robjects import pandas2ri  
                from rpy2.robjects.conversion import localconverter         
                with localconverter(ro.default_converter + pandas2ri.converter):  
                    train_rdf = pandas2ri.py2rpy(train_df)
                    test_rdf = pandas2ri.py2rpy(test_df)
                    
                ro.r('''Sys.setlocale("LC_ALL", "en_US.UTF-8")''')    
                # load the r file  
                ro.r['source']('r_baseline.R')  
                
                # load and exec r function  
                r_function = ro.globalenv['Causal_rule_ensemble']   
                result = r_function(train_rdf, test_rdf, covariates, treatment, outcome)
                
                # change dataframe to python
                with localconverter(ro.default_converter + pandas2ri.converter):  
                    rule_df = ro.conversion.rpy2py(result[0])  

                ate = result[1]
                rules = rule_df["Rule"].tolist()    
                cate_list = rule_df["Estimate"].tolist() 
                
                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rules, cate = cate_list, ate = ate)
                metrics["train_time"] = result[2]
                
            if model == "DT":
                rules, train_time = baseline.dt(train_df , covariates, treatment, outcome, k, thre)
                
                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                            treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rules, cate = False)
                metrics["train_time"] = train_time

            if model == "CF":
                import rpy2.robjects as ro  
                from rpy2.robjects import pandas2ri  
                from rpy2.robjects.conversion import localconverter         
                with localconverter(ro.default_converter + pandas2ri.converter):  
                    train_rdf = pandas2ri.py2rpy(train_df)
                    test_rdf = pandas2ri.py2rpy(test_df)
                    
                ro.r('''Sys.setlocale("LC_ALL", "en_US.UTF-8")''')    
                # load the r file  
                ro.r['source']('r_baseline.R')  

                # load and exec r function  
                print(ro.globalenv)
                r_function = ro.globalenv['causal_forest']   
                result = r_function(train_rdf, test_rdf, covariates, treatment, outcome, k, thre,**kwargs)

                # process the error
                if result[0][0] == "rules not found enough":
                    print("rules not found enough")
                    continue
                
                # change dataframe to python 
                with localconverter(ro.default_converter + pandas2ri.converter):  
                    rules = ro.conversion.rpy2py(result[0])
                
                print(str(rules))
                # change rules to python grammar
                import re    
                rules = rules.reset_index(drop=True)
                for i in range(len(rules)):   
                    rules.loc[i, "rule"] = re.sub(r'%in% c\((.*?)\)', r'in [\1]', rules.loc[i, "rule"])   
                rule_list = list(rules["rule"])
                rule_list = [rule.replace("&", "and") for rule in rule_list]
                        
                cate_list = rules["treatment_effect"].tolist() 
                
                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                    treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rule_list, cate = False)
                metrics["train_time"] = result[1]

            if model == "CURLS":
                rule_list, ite_pre, train_time = baseline.CURLS( train_df, covariates, treatment, outcome, k, thre, **kwargs)
                
                var_wt = kwargs.get('variance_weight', 0)
                num_thresh = kwargs.get('num_thresh', 8)
                max_rule_length = kwargs.get('max_rule_length', 4)
                minimum_coverage = thre

                new_model = model + f" var_wt={var_wt}" + f" n_thre={num_thresh}" + f" max_rule_len={max_rule_length}" + f" min_cover={minimum_coverage}"

                metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = new_model, rule_num = k, thre = thre, outcome_col = outcome, 
                                    treatment_col = treatment,covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rule_list, cate = ite_pre)
                metrics["train_time"] = train_time  

            metrics_results.append(metrics)  
            rules_results.append(rules.drop(columns=['rules']))       
        
        if len(rules_results) == 1:
            print("only one fold result exist")
        
        metrics_df = pd.concat(metrics_results)  
        rules_df = pd.concat(rules_results)  
    
        # calculate the average and standard deviation from the folds
        avg_metrics = metrics_df.groupby(metrics_df.index).mean()  
        avg_rules = rules_df.groupby(rules_df.index).mean()  
        std_metrics = metrics_df.groupby(metrics_df.index).std()  
        std_rules = rules_df.groupby(rules_df.index).std()  
        avg_metrics = avg_metrics.round(3)
        avg_rules = avg_rules.round(3)
        std_metrics = std_metrics.round(3)
        std_rules = std_rules.round(3)   
        
        # generate the result
        result_rules = avg_rules.astype(str) + " (" + std_rules.astype(str) + ")"         
        result_metrics = avg_metrics.astype(str) + " (" + std_metrics.astype(str) + ")"

        result = {
            "rules": result_rules,
            "metrics": result_metrics,
            "avg_metrics": avg_metrics,
            "avg_rules": avg_rules
        }
        
        return result
    else:

        if model == "CT":
            import rpy2.robjects as ro  
            from rpy2.robjects import pandas2ri  
            from rpy2.robjects.conversion import localconverter         
            with localconverter(ro.default_converter + pandas2ri.converter):  
                train_rdf = pandas2ri.py2rpy(train_df)
                test_rdf = pandas2ri.py2rpy(test_df)
                
            ro.r('''Sys.setlocale("LC_ALL", "en_US.UTF-8")''')    
            # load the r file  
            ro.r['source']('r_baseline.R')  

            # load and exec r function  
            r_function = ro.globalenv['causal_tree']   
            result = r_function(train_rdf, test_rdf, covariates, treatment, outcome, k, thre,**kwargs)
            
            # change dataframe to python 
            with localconverter(ro.default_converter + pandas2ri.converter):  
                rules = ro.conversion.rpy2py(result[0])

            if result[0][0] == "rules not found enough":
                print("rules not found enough")
                metrics = pd.DataFrame(index = [model])
                rules = pd.DataFrame(index = [model])
                return rules, metrics 

            # change rules to python grammar
            import re    
            rules = rules.reset_index(drop=True)
            for i in range(len(rules)):   
                rules.loc[i, "rule"] = re.sub(r'%in% c\((.*?)\)', r'in [\1]', rules.loc[i, "rule"])   
            rule_list = list(rules["rule"])
            rule_list = [rule.replace("&", "and") for rule in rule_list]
                    
            cate_list = rules["treatment_effect"].tolist() 
            
            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rule_list, cate = cate_list)
            metrics["train_time"] = result[1]
            
        if model == "RIPPER":
            rules, train_time  = baseline.ripper(train_df , covariates, treatment, outcome) 
            if len(rules) == 0:
                print("no rules found")
                metrics = pd.DataFrame(index = [model])
                rules = pd.DataFrame(index = [model])
                return rules, metrics 
                        
            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                    treatment_col = treatment,covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = str(rules), cate = False)
            metrics["train_time"] = train_time

        if model == "BRCG":
            rules, train_time = baseline.brcg(train_df , covariates, treatment, outcome)
            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                treatment_col = treatment,covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = str(rules), cate = False)
            metrics["train_time"] = train_time
            
        if model == "PYS":
            
            rules, train_time = baseline.pys(train_df , covariates, treatment, outcome)
            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                    treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rules, cate = False)
            metrics["train_time"] = train_time

        if model == "CRE":
            import rpy2.robjects as ro  
            from rpy2.robjects import pandas2ri  
            from rpy2.robjects.conversion import localconverter         
            with localconverter(ro.default_converter + pandas2ri.converter):  
                train_rdf = pandas2ri.py2rpy(train_df)
                test_rdf = pandas2ri.py2rpy(test_df)
                
            ro.r('''Sys.setlocale("LC_ALL", "en_US.UTF-8")''')    
            # load the r file  
            ro.r['source']('r_baseline.R')  
            
            # load and exec r function  
            r_function = ro.globalenv['Causal_rule_ensemble']   
            result = r_function(train_rdf, test_rdf, covariates, treatment, outcome)
            
            # change dataframe to python
            with localconverter(ro.default_converter + pandas2ri.converter):  
                rule_df = ro.conversion.rpy2py(result[0])  

            ate = result[1]
            rules = rule_df["Rule"].tolist()    
            cate_list = rule_df["Estimate"].tolist() 
            
            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                            treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rules, cate = cate_list, ate = ate)
            metrics["train_time"] = result[2]
            
        if model == "DT":
            rules, train_time = baseline.dt(train_df , covariates, treatment, outcome, k, thre)
            
            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                        treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rules, cate = False)
            metrics["train_time"] = train_time

        if model == "CF":
            import rpy2.robjects as ro  
            from rpy2.robjects import pandas2ri  
            from rpy2.robjects.conversion import localconverter         
            with localconverter(ro.default_converter + pandas2ri.converter):  
                train_rdf = pandas2ri.py2rpy(train_df)
                test_rdf = pandas2ri.py2rpy(test_df)
                
            ro.r('''Sys.setlocale("LC_ALL", "en_US.UTF-8")''')    
            # load the r file  
            ro.r['source']('r_baseline.R')  

            # load and exec r function  
            print(ro.globalenv)
            r_function = ro.globalenv['causal_forest']   
            result = r_function(train_rdf, test_rdf, covariates, treatment, outcome, k, thre,**kwargs)
            
            # change dataframe to python 
            with localconverter(ro.default_converter + pandas2ri.converter):  
                rules = ro.conversion.rpy2py(result[0])

            if result[0][0] == "rules not found enough":
                print("rules not found enough")
                metrics = pd.DataFrame(index = [model])
                rules = pd.DataFrame(index = [model])
                return rules, metrics 
                        
            # change rules to python grammar
            import re    
            rules = rules.reset_index(drop=True)
            for i in range(len(rules)):   
                rules.loc[i, "rule"] = re.sub(r'%in% c\((.*?)\)', r'in [\1]', rules.loc[i, "rule"])   
            rule_list = list(rules["rule"])
            rule_list = [rule.replace("&", "and") for rule in rule_list]
                    
            cate_list = rules["treatment_effect"].tolist() 
            
            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = model, rule_num = k, thre = thre, outcome_col = outcome, 
                                treatment_col = treatment, covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rule_list, cate = False)
            metrics["train_time"] = result[1]

        if model == "CURLS":
            rule_list, ite_pre, train_time = baseline.CURLS( train_df, covariates, treatment, outcome, k, thre, **kwargs)

            var_wt = kwargs.get('variance_weight', 0)
            num_thresh = kwargs.get('num_thresh', 8)
            max_rule_length = kwargs.get('max_rule_length', 4)
            minimum_coverage = thre

            new_model = model + f" var_wt={var_wt}" + f" n_thre={num_thresh}" + f" max_rule_len={max_rule_length}" + f" min_cover={minimum_coverage}"

            metrics, rules = metric.get_Metric(train_df = train_df, test_df = test_df, model = new_model, rule_num = k, thre = thre, outcome_col = outcome, 
                                treatment_col = treatment,covariates = covariates, te_col = 'TE', te_pre_col = "te_pre", rules = rule_list, cate = ite_pre)
            metrics["train_time"] = train_time
        
        result = {
            "rules": rules,
            "metrics": metrics
        }

        return result
