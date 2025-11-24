import re
import pandas as pd
import numpy as np

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV

import baseline 

import sys
sys.path.append('./CURLS/osfstorage-archive/')
# sys.path.append('./CURLS/osfstorage-archive/pipeline/')



# pre-compile condition pattern once
_COND = re.compile(r'^\s*(\w+)\s*(>=|<=|!=|=|>|<)\s*(.+?)\s*$')


def fit_rf_propensity(
    X_train: pd.DataFrame | np.ndarray,
    X_test: pd.DataFrame | np.ndarray,
    treatment: np.ndarray,
    *,
    n_estimators: int = 100,
    max_depth: int | None = None,
    min_samples_leaf: int = 20,
    random_state: int | None = None,
):
    """
    Random-forest propensity-score model with isotonic calibration.
    Returns the calibrated model and propensity scores.
    """
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_samples_leaf=min_samples_leaf,
        random_state=random_state,
    )
    calibrated = make_pipeline(
        StandardScaler(with_mean=False),  # keeps sparse matrices sparse
        CalibratedClassifierCV(rf, cv=5, method="isotonic"),
    )
    calibrated.fit(X_train, treatment)
    ps_train = calibrated.predict_proba(X_train)[:, 1]
    ps_test = calibrated.predict_proba(X_test)[:, 1]
    return calibrated, ps_train, ps_test

def _cond_mask(df: pd.DataFrame, text: str) -> pd.Series:
    """
    Convert a *single* condition (no “and”) into a boolean mask.
    Handles:
        • comparisons, e.g. “X7 > 1.21”
        • bare columns,  e.g. “X13”  (treated as truthy if non-zero / True)
    """
    text = text.strip()
    m = _COND.match(text)
    if m:                                  # comparison
        col, op, val = m.groups()

        # try numeric, else leave as string
        try:
            val = float(val)
        except ValueError:
            val = val.strip('\'"')          # unquote bare strings

        s = df[col]
        return {">":  s >  val,
                "<":  s <  val,
                ">=": s >= val,
                "<=": s <= val,
                "=": s == val,
                "!=": s != val}[op]
    else:                                   # bare column
        return df[text].astype(bool)


def subgroup(df: pd.DataFrame, rule: str) -> pd.DataFrame:
    """
    Return the rows of *df* that satisfy a single rule string.

    Example
    -------
    >>> subgroup(df, 'X7 > 1.21 and X9 > -1.22 and X10 <= 1.21 and X13')
    """
    mask = pd.Series(True, index=df.index)
    for part in rule.split('and'):
        mask &= _cond_mask(df, part)
    return mask.to_numpy()


def read_and_process_data(path, simulator_name, n, iter):
    
    data_name = path + simulator_name + "/n_" + str(n) + "_iter_" + str(iter) + ".csv"

    print("running curls on data:", data_name)
    df = pd.read_csv(path)

    # change the column name 'treatment' to 't'
    df = df.rename(columns={"T": "t"})
    df = df.rename(columns={"Y": "y"})

    np.random.seed(1)
    train_indices = np.random.choice(df.index, size=int(0.5 * len(df)), replace=False)
    train_df = df.loc[train_indices].reset_index(drop=True)
    test_df = df.drop(train_indices).reset_index(drop=True)

    # use train_df to fit the propensity model
    X_train = train_df.drop(columns=['t', 'y'])
    treatment_train = train_df['t'].values
    X_test = test_df.drop(columns=['t', 'y'])
    rf_model, propensity_scores_train, propensity_scores_test = fit_rf_propensity(
        X_train, X_test, treatment_train, n_estimators=100, max_depth=None, min_samples_leaf=20, random_state=42
    )
    # calculate the "weights" as T/propensity_scores + (1-T)/(1-propensity_scores)
    train_df['wt'] = train_df['t'] / propensity_scores_train + (1 - train_df['t']) / (1 - propensity_scores_train)
    test_df['wt'] = test_df['t'] / propensity_scores_test + (1 - test_df['t']) / (1 - propensity_scores_test)

    return df, train_df, test_df, train_indices


def get_ground_truth_subgroup_bool(df, simulator_name):
    simulator_names = (
        "simulate1",
        "simulate_imbalance_treatment",
        "simulate_long_rule"
    )
    assert simulator_name in simulator_names, f"Simulator name {simulator_name} not recognized."

    if simulator_name in ("simulate1", "simulate_imbalance_treatment"):
        return df["X1"] > 1
    elif simulator_name == "simulate_long_rule":
        return (df["X1"] > -1) & (df["X2"] > -1) & (df["X3"] > -1)
    else:
        raise ValueError("Unknown simulator name")


def gt_te_per_sample(d, simulator_name):
    """
    Calculate ground truth treatment effect per sample.
    
    Parameters:
    -----------
    d : pandas.DataFrame
        DataFrame containing covariates X1, X2, X3
    simulator_name : str
        Name of the simulator ("simulate1", "simulate_imbalance_treatment", or "simulate_long_rule")
    
    Returns:
    --------
    numpy.ndarray
        Array of treatment effects (tau) for each sample
    """
    if simulator_name in ["simulate1", "simulate_imbalance_treatment"]:
        # Same outcome mechanism; only treatment assignment differs.
        mu1 = np.where(d['X1'] > 1, 0.8, 0.2)   # E[Y | T=1, X]
        mu0 = np.where(d['X1'] < -1, 0.75, 0.2) # E[Y | T=0, X]
        tau = mu1 - mu0
        return tau
    
    elif simulator_name == "simulate_long_rule":
        # Recreate the rules exactly as in the simulator
        rule1 = (d['X1'] > -1) & (d['X2'] > -1) & (d['X3'] > -1)
        rule2 = (d['X1'] > -1) & (d['X2'] > -1) & (~rule1)
        rule3 = (d['X1'] > -1) & (~rule1) & (~rule2)
        
        # E[Y | T=1, X] by rule; baseline is 0.2 when no rule fires
        mu1 = np.where(rule1, 0.8,
                      np.where(rule2, 0.6,
                              np.where(rule3, 0.4, 0.2)))
        
        # Under T=0 everything is baseline 0.2 in this simulator
        mu0 = np.full(len(d), 0.2)
        
        tau = mu1 - mu0
        return tau
    
    else:
        raise ValueError(f"Unknown simulator_name: {simulator_name}")

def fit_curls(train_df, test_df, simulator_name):
    colnames = train_df.columns.tolist()
    covariates = [col for col in colnames if col not in ['t', 'y', "wt",  "TE"]]
    treatment = 't'
    outcome = 'y'

    num_rules_to_find=2
    rules_least_cov_ratio=0.03 

    variance_weights = np.linspace(0.1, 1.5, num=10)
    best_ite_pre = -np.inf
    for variance_weight in variance_weights:
        rule_list, ite_pre, train_time = baseline.CURLS(train_df, covariates, treatment, outcome, k=num_rules_to_find, thre=rules_least_cov_ratio, 
                                                        variance_weight=variance_weight, num_thresh=14, max_rule_length=train_df.shape[1])
        if np.max(ite_pre) > best_ite_pre:
            which_best = np.argmax(ite_pre)

            best_ite_pre = ite_pre[which_best]
            best_rule_list = rule_list
            best_rule = best_rule_list[which_best]
            best_train_time = train_time
            best_variance_weight = variance_weight
    
    rule = best_rule
    mask = subgroup(test_df, rule)
    train_time = best_train_time
    

    # we need the following for our experiments: 
    #   - the treatment effect of the subgroup estimated from the test set
    #   - the ground truth treatment effect of the subgroup, average and variance. 
    #   - the difference between the gt_average and the estimated treatment effect. 
    #   - the difference between the estimated treatment of the subgroup and the individual treatment effect of the subgroup.
    mask_t1 = mask & (test_df['t'] == 1)
    mask_t0 = mask & (test_df['t'] == 0)

    if np.sum(mask_t1) == 0:
        mask_t1 = (test_df['t'] == 1)  # if no treatment in the subgroup, use all treated
    if np.sum(mask_t0) == 0:
        mask_t0 = (test_df['t'] == 0)

    # get the ground truth subgroup 
    gt_bool = get_ground_truth_subgroup_bool(test_df, simulator_name)  # change this to the appropriate simulator name

    treatment_effect_estimated = np.mean(test_df.loc[mask_t1, 'y']) - np.mean(test_df.loc[mask_t0, 'y'])

    jaccard_sim = np.sum(mask & gt_bool) / np.sum(mask | gt_bool) if np.sum(mask | gt_bool) > 0 else 0

    theoretical_gt_te_per_sample = gt_te_per_sample(test_df, simulator_name)
    gt_te_of_gt_subgroup_ = np.mean(theoretical_gt_te_per_sample[gt_bool])
    gt_te_of_learned_subgroup_ = np.mean(theoretical_gt_te_per_sample[mask])
    estimated_te_of_learned_subgroup_ = treatment_effect_estimated
    estimated_te_of_gt_subgroup_ = np.mean(test_df.loc[gt_bool & (test_df['t'] == 1), 'y']) - np.mean(test_df.loc[gt_bool & (test_df['t'] == 0), 'y'])

    return {
        'treatment_effect_estimated': treatment_effect_estimated,
        'train_time': train_time,
        'jaccard_similarity': jaccard_sim,
        'gt_te_of_gt_subgroup': gt_te_of_gt_subgroup_,
        'gt_te_of_learned_subgroup': gt_te_of_learned_subgroup_,
        'estimated_te_of_learned_subgroup': estimated_te_of_learned_subgroup_,
        'estimated_te_of_gt_subgroup': estimated_te_of_gt_subgroup_
    }


if __name__ == '__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser(description="Run CURLS on semi-synthetic data.")
    parser.add_argument('--data_path', type=str, default=None, help='Path to the data file.')
    args = parser.parse_args()

    if args.data_path is None:
        args.data_path = "../../../datasets/synthetic/simulate1/n_1000_iter_1.csv"
    
    # Extract simulator name, n, and iter from the data path
    fname  = os.path.basename(args.data_path)              # n_1000_iter_1.csv
    match  = re.match(r'n_(\d+)_iter_(\d+)\.csv', fname)
    sim    = os.path.basename(os.path.dirname(args.data_path))

    n = int(match.group(1)) if match else 1000  # default to 1000 if not found
    iter = int(match.group(2)) if match else 1  # default to 1 if not found
    simulator_name = sim 

    # Read and process the data
    df, train_df, test_df, train_indices = read_and_process_data(args.data_path, simulator_name, n, iter)

    # Fit CURLS and get results
    results = fit_curls(train_df, test_df, simulator_name=simulator_name)

    # save the results to a csv file
    results_df = pd.DataFrame([results])
    data_name = args.data_path.split("/")[-1].split(".")[0] if args.data_path else "IHDP"
    save_name = data_name + "_curls_results.csv"
    save_folder = "../pipeline/myresults/semi-synthetic/"
    date_today = pd.Timestamp.now().strftime('%Y-%m-%d')
    save_folder =save_folder + date_today + "/"
    # create the folder if it does not exist
    import os
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    results_df.to_csv(os.path.join(save_folder, save_name), index=False)





