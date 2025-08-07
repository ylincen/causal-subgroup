from sklearn.model_selection import train_test_split
from aix360.algorithms.rbm import FeatureBinarizer
import pandas as pd
import numpy as np
import time
import warnings  
warnings.filterwarnings("ignore", category=UserWarning)  

def ripper( train_df, covariates, treatment, outcome):

    from aix360.algorithms.rule_induction.ripper import RipperExplainer

    # prepare data
    train_df = binary_outcome(train_df, outcome)
    columns = covariates + [outcome]
    train_df = train_df[columns]
    TARGET_COLUMN = outcome
    POS_VALUE = 1
    y_train = train_df[TARGET_COLUMN]
    x_train = train_df.drop(columns=[TARGET_COLUMN])

    #train and get the result
    estimator = RipperExplainer()
    start_time = time.time()
    estimator.fit(x_train, y_train, target_label=POS_VALUE)
    end_time = time.time()
    train_time = end_time - start_time
    
    rules = estimator.explain()
    
    return rules, train_time

def brcg( train_df, covariates, treatment, outcome):
    
    from aix360.algorithms.rule_induction.rbm.boolean_rule_cg import BooleanRuleCG
    
    # prepare data
    train_df = binary_outcome(train_df, outcome)
    columns = covariates + [outcome]
    train_df = train_df[columns]
    TARGET_COLUMN = outcome
    POS_VALUE = 1
    y_train = train_df[TARGET_COLUMN]
    x_train = train_df.drop(columns=[TARGET_COLUMN])
    
    #train and get the result
    fb = FeatureBinarizer(negations=True)
    X_train_fb = fb.fit_transform(x_train)  
    explainer = BooleanRuleCG(silent=True)
    start_time = time.time()
    explainer.fit(X_train_fb, y_train)
    end_time = time.time()
    train_time = end_time - start_time
    
    rules = explainer.explain()
    
    return rules, train_time

def pys( train_df, covariates, treatment, outcome):
   
    import pysubgroup as ps
    train_df = binary_outcome(train_df, outcome)
    target = ps.BinaryTarget (outcome, True)
    non_covariate_columns = [col for col in train_df.columns if col not in covariates]  
    searchspace = ps.create_selectors(train_df, ignore=non_covariate_columns)
    task = ps.SubgroupDiscoveryTask (
        train_df,
        target,
        searchspace,
        result_set_size=10,
        depth=5,
        qf=ps.WRAccQF())
    
    start_time = time.time()
    result = ps.DFS().execute(task)
    end_time = time.time()
    train_time = end_time - start_time
    
    result = result.to_dataframe()
    
    # process the rules
    rules = result["subgroup"].astype(str).tolist() 
    processed_rules = pys_process(rules)
    
    return processed_rules, train_time

def dt( train_df, covariates, treatment, outcome, k, thre, max_depth = 4):
    from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier   
    from sklearn.preprocessing import LabelEncoder 

    # prepare data
    train_df = train_df.copy()
    encoders = {}     
    for column in train_df.columns:    
        if train_df[column].dtype == 'object':  # 如果是类别型特征    
            le = LabelEncoder()  
            train_df[column] = le.fit_transform(train_df[column])  # 进行标签编码    
            encoders[column] = le  # 保存LabelEncoder对象  
     
    X_train = train_df[covariates]
    y_train = train_df[outcome]  

    # train the model 
    model = DecisionTreeRegressor(max_depth=max_depth, random_state=42)  
    start_time = time.time()
    model.fit(X_train, y_train)  
    end_time = time.time()
    train_time = end_time - start_time

    # get the rules
    rules = get_reg_rules(model.tree_, covariates, thre) 
    rule_list = [rule[:-5] for prob, rule in rules]
    rule_list = [decode_rule(rule, encoders) for rule in rule_list]

    return rule_list, train_time

def CURLS( train_df, covariates, treatment, outcome, k, thre, max_rule_length = 4, num_thresh=8, variance_weight = 0.0):
    import sys
    sys.path.append("/Users/yanglincen/projects/causal_sub_R/competitors/curls/osfstorage-archive")
    from CURLS import DataPreprocessor, CausalRuleLearner
    
    train_df = train_df.reset_index(drop=True)  

    col_categ = []
    for var in covariates:   
        if pd.api.types.is_object_dtype(train_df[var]):   
            col_categ.append(var)

    data_preprocessor = DataPreprocessor(covariates, col_categ, num_thresh = num_thresh)
    df_binarized = data_preprocessor.fit_transform(train_df)
    causal_learner = CausalRuleLearner(df_binarized, max_rule_length = max_rule_length,variance_weight = variance_weight, minimum_coverage = thre)
    
    start_time = time.time()
    causal_learner.find_rules(num_rules_to_find = k)
    end_time = time.time()
    train_time = end_time - start_time
    
    causal_learner.print_rules()
    rule_list = causal_learner.conditions_str
    ite_list = causal_learner.treatment_effect_values

    return rule_list, ite_list, train_time
    

#process the rules from pysubgroup
def pys_process(rules):  
    processed_rules = []  
    for rule in rules:  
        if "AND" in rule:  
            sub_rules = rule.split("AND")  
            processed_sub_rules = [pys_process([sub_rule.strip()])[0] for sub_rule in sub_rules]  
            processed_rules.append(" AND ".join(processed_sub_rules))  
        elif rule.count(":") == 2:  
            column, range_str = rule.split(":", 1)
            lower, upper = range_str.strip("[] ").split(":")
            processed_rules.append(f"{column.strip()} >= {lower.strip()} AND {column.strip()} < {upper.strip()}")  
        else:  
            processed_rules.append(rule)  
    return [rule.replace("AND", "and") for rule in processed_rules]

# get the rules from decision tree
def get_reg_rules(tree, feature_names, thre, node_id=0, decision_path=None):  

    decision_paths = []  
      
    if tree.children_left[node_id] == tree.children_right[node_id]:  
  
        if tree.n_node_samples[node_id] > tree.n_node_samples[0] * thre: 

            avg_value = tree.value[node_id][0]
            decision_paths.append((avg_value, decision_path))  
    else:  
        # 如果不是叶子节点，递归处理左右子节点  
        feature_name = feature_names[tree.feature[node_id]]  
        threshold = round(tree.threshold[node_id], 2)  # 保留m位有效数字  
          
        if decision_path is None:  
            decision_path = ''  
          
        left_decision_path = decision_path + f"{feature_name} <= {threshold} and "  
        decision_paths.extend(get_reg_rules(tree, feature_names, thre, tree.children_left[node_id], left_decision_path))  
          
        right_decision_path = decision_path + f"{feature_name} > {threshold} and "  
        decision_paths.extend(get_reg_rules(tree, feature_names, thre, tree.children_right[node_id], right_decision_path))  
      
    return decision_paths  
    
# Function to convert the outcome data to binary 0 1  
def binary_outcome(df, outcome):    
    df_copy = df.copy()  
      
    # Check if outcome exists in the dataframe  
    if outcome not in df_copy.columns:    
        print(f"Error: No '{outcome}' column in the dataframe")    
        return df_copy    
    
    column_dtype = df_copy[outcome].dtype    
    
    # Check if the outcome column is numeric  
    if np.issubdtype(column_dtype, np.number):    
        unique_values = df_copy[outcome].unique()    
        # Check if the outcome is continuous or non-binary  
        if len(unique_values) > 2 or (len(unique_values) == 2 and set(unique_values) != {0, 1}):    
            # Convert continuous variable to binary  
            median = df_copy[outcome].median()    
            df_copy[outcome] = df_copy[outcome].apply(lambda x: 1 if x > median else 0)    
            print(f"Continuous variable converted to binary. Values greater than the median {median} are marked as 1, others as 0")    
    else:    
        # Handle non 0 1 binary variables  
        value_counts = df_copy[outcome].value_counts()    
        majority_class = value_counts.idxmax()    
        minority_class = value_counts.idxmin()    
    
        # Check if the outcome is non-binary  
        if len(value_counts) > 2 or (len(value_counts) == 2 and set(df_copy[outcome].unique()) != {0, 1}):    
            df[outcome] = df[outcome].apply(lambda x: 1 if x == majority_class else 0)    
            print(f"Binary variable adjusted, majority class {majority_class} marked as 1, minority class {minority_class} marked as 0")    
        
    return df_copy    
  
# Function to convert the encoded values in the rules back to their categorical form  
def decode_rule(rule, encoders):    
    # Create a dictionary to map all possible category values of each feature to their corresponding encoded values  
    category_maps = {feature: list(encoder.classes_) for feature, encoder in encoders.items()}    
    
    # Split the rule into a list of conditions  
    conditions = rule.split(' and ')    
    
    # Process each condition  
    for i in range(len(conditions)):    
        # Split the condition into feature, comparator, and value  
        feature, comp, value = conditions[i].split(' ')    
    
        # Only convert when the feature exists in the encoders dictionary  
        if feature in encoders:    
            # Convert the value to integer  
            value = int(float(value))    
    
            # Convert the value back to the original category label  
            value = encoders[feature].inverse_transform([value])[0]    
    
            # Get all possible category values of this feature  
            categories = category_maps[feature]    
    
            # Determine all possible category values included in this condition based on the comparator  
            if comp == '<=':    
                included_categories = categories[:categories.index(value) + 1]    
                conditions[i] = f'{feature} in {included_categories}'    
            elif comp == '>':    
                included_categories = categories[categories.index(value) + 1:]    
                conditions[i] = f'{feature} in {included_categories}'    
    
    # Reassemble the rule  
    rule = ' and '.join(conditions)    
    
    return rule    