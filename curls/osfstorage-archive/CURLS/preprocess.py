import pandas as pd
from .features import FeatureBinarizer

class DataPreprocessor:
    def __init__(
        self,
        covariate_columns,
        col_categ=None,
        num_thresh=4,
        negations=True,
        thresh_str=True,
    ):
        """
        Initialize the DataPreprocessor with given columns and binarization parameters.

        :param covariate_columns: List of names of the covariate columns.
        :param col_categ: List of categorical column names.
        :param num_thresh: The number of thresholds for binarizing numerical features.
        :param negations: Whether to include negations for binarized features.
        :param thresh_str: Whether to represent thresholds as strings.
        """
        self.covariate_columns = covariate_columns
        self.col_categ = col_categ
        self.num_thresh = num_thresh
        self.negations = negations
        self.thresh_str = thresh_str
        self.binarizer = FeatureBinarizer(
            colCateg=self.col_categ,
            numThresh=self.num_thresh,
            negations=self.negations,
            threshStr=self.thresh_str,
        )

    def fit_transform(self, df):
        """
        Fit the binarizer to the data and transform the specified covariate columns,
        then return the original DataFrame with transformed covariates.

        :param df: DataFrame containing the data to be preprocessed.
        :return: DataFrame with binarized covariate columns merged with the original data.
        """

        # Temporary: Check if 'v' column exists; if not, create it as wt * y
        if "v" not in df.columns:
            df["v"] = df["wt"] * df["y"]

        # Select only the specified covariate columns for binarization
        covariate_df = df[self.covariate_columns]

        # Perform binarization on the covariates
        df_binarized = self.binarizer.fit_transform(covariate_df)
        df_binarized.columns = [
            " ".join(col).strip() for col in df_binarized.columns.values
        ]

        # Drop the original covariate columns from the input df
        df_dropped = df.drop(columns=self.covariate_columns)

        # Concatenate the binarized covariates back into the original DataFrame
        df_final = pd.concat([df_binarized, df_dropped], axis=1)

        return df_final
