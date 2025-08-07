import pandas as pd 
import numpy as np
from sklearn.model_selection import train_test_split
import metric
import baseline 
from sklearn.model_selection import KFold  


from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import RandomForestClassifier
from sklearn.calibration import CalibratedClassifierCV

from run_semi import *

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

df = pd.read_csv("../pipeline/mydata/semi-synthetic/IHDP.csv")
y1 = np.array(df["y_factual"])
y1[df['treatment'] == 0] = df["y_cfactual"][df['treatment'] == 0]
y0 = np.array(df["y_factual"])
y0[df['treatment'] == 1] = df["y_cfactual"][df['treatment'] == 1]

df['TE'] = y1 - y0
df_y_cfactual = df["y_cfactual"]
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


colnames = train_df.columns.tolist()
covariates = [col for col in colnames if col not in ['t', 'y', "wt",  "TE"]]
treatment = 't'
outcome = 'y'

num_rules_to_find=1
rules_least_cov_ratio=0.03 
variance_weight=0.43

rule_list, ite_pre, train_time = baseline.CURLS(train_df, covariates, treatment, outcome, k=num_rules_to_find, thre=rules_least_cov_ratio, 
                                                variance_weight=0.43, num_thresh=14, max_rule_length=df.shape[1])

# save the rule_list, ite_pre, and train_time to a file
import pickle
with open('rule_list.pkl', 'wb') as f:
    pickle.dump(rule_list, f)
with open('ite_pre.pkl', 'wb') as f:
    pickle.dump(ite_pre, f)
with open('train_time.pkl', 'wb') as f:
    pickle.dump(train_time, f)

import pickle
with open('rule_list.pkl', 'rb') as f:
    rule_list = pickle.load(f)

# check the mask for the first rule
mask = subgroup(test_df, rule_list[0])
