import re
import pandas as pd
import numpy as np

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV

import baseline 

import sys
sys.path.append('./curls/osfstorage-archive/')

# from CURLS import DataPreprocessor, CausalRuleLearner


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


def read_and_process_data(path):
    if path is None:
        path = "../../../datasets/semi_synthetic/ACIC2016_17.csv"
    
    data_name = path.split("/")[-1].split(".")[0]

    print("running curls on data:", data_name)
    df = pd.read_csv(path)
    y1 = np.array(df["y_factual"])
    y1[df['treatment'] == 0] = df["y_cfactual"][df['treatment'] == 0]
    y0 = np.array(df["y_factual"])
    y0[df['treatment'] == 1] = df["y_cfactual"][df['treatment'] == 1]

    df['TE'] = y1 - y0
    df["y"] = df["y_factual"]
    df = df.drop(columns = ["y_factual", "y_cfactual"])
    # print(df.head())

    # change the column name 'treatment' to 't'
    df = df.rename(columns={"treatment": "t"})

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


def fit_curls(train_df, test_df):
    colnames = train_df.columns.tolist()
    covariates = [col for col in colnames if col not in ['t', 'y', "wt",  "TE"]]
    treatment = 't'
    outcome = 'y'

    num_rules_to_find=2
    rules_least_cov_ratio=0.03 

    variance_weights = np.linspace(0.1, 1.5, num=5)
    variance_weights = [0.2]
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

    mask_t1 = mask & (test_df['t'] == 1)
    mask_t0 = mask & (test_df['t'] == 0)

    if np.sum(mask_t1) == 0:
        mask_t1 = (test_df['t'] == 1)  # if no treatment in the subgroup, use all treated
    if np.sum(mask_t0) == 0:
        mask_t0 = (test_df['t'] == 0)

    treatment_effect_estimated = np.mean(test_df.loc[mask_t1, 'y']) - np.mean(test_df.loc[mask_t0, 'y'])
    treatment_effect_gt = np.mean(test_df.loc[mask, 'TE'])
    treatment_effect_gt_var = np.var(test_df.loc[mask, 'TE'])
    diff_gt_estimated = treatment_effect_gt - treatment_effect_estimated
    diff_estimated_ite_mae = np.mean(abs(treatment_effect_estimated - test_df.loc[mask, 'y']))  # estimated ITE for the subgroup
    diff_esimatedd_ite_var = np.var(treatment_effect_estimated - test_df.loc[mask, 'y'])
    return {
        'treatment_effect_estimated': treatment_effect_estimated,
        'treatment_effect_gt': treatment_effect_gt,
        'treatment_effect_gt_var': treatment_effect_gt_var,
        'diff_gt_estimated': diff_gt_estimated,
        'diff_estimated_ite_mae': diff_estimated_ite_mae,
        'diff_estimated_ite_var': diff_esimatedd_ite_var,
        'train_time': train_time
    }

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Run CURLS on semi-synthetic data.")
    parser.add_argument('--data_path', type=str, default=None, help='Path to the data file.')
    args = parser.parse_args()

    # Read and process the data
    df, train_df, test_df, train_indices = read_and_process_data(args.data_path)

    # Fit CURLS and get results
    results = fit_curls(train_df, test_df)

    # save the results to a csv file
    results_df = pd.DataFrame([results])
    data_name = args.data_path.split("/")[-1].split(".")[0] if args.data_path else "IHDP"
    save_name = data_name + "_curls_results.csv"
    save_folder = "./pipeline/myresults/semi-synthetic/"
    # create the folder if it does not exist
    import os
    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    results_df.to_csv(os.path.join(save_folder, save_name), index=False)





