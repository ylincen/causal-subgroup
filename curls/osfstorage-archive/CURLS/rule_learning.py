import pandas as pd
import numpy as np

class CausalRuleLearner:
    def __init__(
        self,
        data,
        treatment_col_name="t",
        outcome_col_name="y",
        treatment_effect_col_name="TE",
        propensity_col_name="e",
        weight_col_name="wt",
        weighted_outcome_col_name="v",
        max_rule_length=None,
        variance_weight=0.0,
        minimum_coverage=0.1,
    ):
        self.data = data
        self.num_samples = data.shape[0]
        self.treatment_col = treatment_col_name
        self.outcome_col = outcome_col_name
        self.treatment_effect_col = treatment_effect_col_name
        self.propensity_col = propensity_col_name
        self.weight_col = weight_col_name
        self.weighted_outcome_col = weighted_outcome_col_name
        self.covariate_cols = [
            col
            for col in data.columns
            if col
            not in [
                self.treatment_col,
                self.outcome_col,
                self.treatment_effect_col,
                self.propensity_col,
                self.weight_col,
                self.weighted_outcome_col,
            ]
        ]
        self.max_rule_length = max_rule_length
        self.rules = []
        self.conditions_str = []
        self.objective_values = []
        self.treatment_effect_values = []
        self.variance_values = []
        self.covered_data_indices = []
        self.mu_t = None
        self.variance_weight = variance_weight
        self.used_features = set()
        self.minimum_coverage = minimum_coverage

    def marginal_gain(self, current_rule, candidate_feature_idx, Q_type, add=True):
        """
        Calculate the marginal gain of adding or removing a single feature to/from the rule set.

        :param current_rule: The current rule represented by a binary array of selected features.
        :param candidate_feature_idx: The index of the feature to be considered for addition or removal.
        :param Q_type: The type of Q value for which to compute the marginal gain
        :param add: Whether to calculate the gain for adding (True) or removal (False).
        :return: The marginal gain value.
        """
        new_rule = current_rule.copy()
        new_rule[candidate_feature_idx] = 1 if add else 0

        Q_current = self.calculate_Q(current_rule, Q_type)
        Q_new = self.calculate_Q(new_rule, Q_type)

        return Q_new - Q_current if add else Q_current - Q_new

    def true_objective_local_search_candidate(
        self, current_solution, max_iteratinos=5, tolerance=1e-6
    ):
        """
        Generate a candidate solution from the current solution using a local search strategy based on the true objective function.

        :param current_solution: The current rule represented by a binary array of selected features.
        :param max_iterations: Maximum number of iterations.
        :param tolerance: Tolerance for convergence.
        :return: A candidate solution and its type of operation ('add', 'remove', or 'replace').
        """
        best_true_objective = self.objective_function(current_solution)[0]
        best_candidate = current_solution.copy()
        improved = True
        iteration_count = 0

        while improved and iteration_count < max_iteratinos:
            improved = False
            best_rule_length = best_candidate.sum()

            # Try adding features
            for i in range(len(best_candidate)):
                if best_candidate[i] == 0 and (
                    self.max_rule_length is None
                    or best_rule_length < self.max_rule_length
                ):
                    candidate_add = current_solution.copy()
                    candidate_add[i] = 1
                    true_obj_add = self.objective_function(candidate_add)[0]
                    if true_obj_add > best_true_objective + tolerance:
                        best_true_objective = true_obj_add
                        best_candidate = candidate_add
                        improved = True

            # Try removing features
            for i in range(len(best_candidate)):
                if best_candidate[i] == 1:
                    candidate_remove = best_candidate.copy()
                    candidate_remove[i] = 0
                    true_obj_remove = self.objective_function(candidate_remove)[0]
                    if true_obj_remove > best_true_objective + tolerance:
                        best_true_objective = true_obj_remove
                        best_candidate = candidate_remove
                        improved = True

            # Try replacing features
            for i in range(len(best_candidate)):
                for j in range(
                    i + 1, len(best_candidate)
                ):  # Avoid duplicate replacements
                    if best_candidate[i] == 1 and best_candidate[j] == 0:
                        candidate_replace = best_candidate.copy()
                        candidate_replace[i] = 0
                        candidate_replace[j] = 1
                        true_obj_replace = self.objective_function(candidate_replace)[0]
                        if true_obj_replace > best_true_objective + tolerance:
                            best_true_objective = true_obj_replace
                            best_candidate = candidate_replace
                            improved = True

            iteration_count += 1

        return best_candidate, best_true_objective

    def local_search_candidate(self, current_solution):
        """
        Generate a candidate solution from the current solution using a local search strategy.

        :param current_solution: The current rule represented by a binary array of selected features.
        :return: A candidate solution and its type of operation ('add', 'remove', or 'replace').
        """
        best_surrogate_object = self.surrogate_objective(
            current_solution, current_solution
        )
        best_candidate = current_solution.copy()
        operation = None
        current_rule_length = current_solution.sum()

        # Try adding features
        for i in range(len(current_solution)):
            if current_solution[i] == 0 and (
                self.max_rule_length is None
                or current_rule_length < self.max_rule_length
            ):
                candidate_add = current_solution.copy()
                candidate_add[i] = 1
                surrogate_obj_add = self.surrogate_objective(
                    current_solution, candidate_add
                )
                if surrogate_obj_add > best_surrogate_object:
                    best_surrogate_object = surrogate_obj_add
                    best_candidate = candidate_add
                    operation = "add"

        # Try removing features
        for i in range(len(current_solution)):
            if current_solution[i] == 1:
                candidate_remove = current_solution.copy()
                candidate_remove[i] = 0
                surrogate_obj_remove = self.surrogate_objective(
                    current_solution, candidate_remove
                )
                if surrogate_obj_remove > best_surrogate_object:
                    best_surrogate_object = surrogate_obj_remove
                    best_candidate = candidate_remove
                    operation = "remove"

        # Try replacing features
        for i in range(len(current_solution)):
            for j in range(
                i + 1, len(current_solution)
            ):  # Avoid duplicate replacements
                if current_solution[i] == 1 and current_solution[j] == 0:
                    candidate_replace = current_solution.copy()
                    candidate_replace[i] = 0
                    candidate_replace[j] = 1
                    surrogate_obj_replace = self.surrogate_objective(
                        current_solution, candidate_replace
                    )
                    if surrogate_obj_replace > best_surrogate_object:
                        best_surrogate_object = surrogate_obj_replace
                        best_candidate = candidate_replace
                        operation = "replace"

        return best_candidate, operation

    def mm_iteration(self, current_solution, max_iterations=10):
        """
        Perform an iteration of the MM algorithm using multiple steps of local search until no further improvement is found or the maximum number of iterations is reached.

        :param current_solution: The current rule represented by a binary array of selected features.
        :param max_iterations: Maximum number of iterations.
        :return: The best solution found within the maximum number of iterations.
        """
        improved = True
        iteration_count = 0
        while improved and iteration_count < max_iterations:
            candidate_solution, operation = self.local_search_candidate(
                current_solution
            )
            if operation is None:
                improved = False
            else:
                current_solution = candidate_solution
                self.mu_t = self.calculate_Q(current_solution, "Q1") / self.calculate_Q(
                    current_solution, "Q3"
                )
                iteration_count += 1
        # print("local_search_count: ", iteration_count)
        return current_solution

    def optimize(self, max_iterations=10, tolerance=1e-6):
        """
        Optimize the causal rule learning using the MM algorithm and submodular optimization.

        :param max_iterations: Maximum number of iterations.
        :param tolerance: Tolerance for convergence.
        :return: Best found solution and its objective function value, treatment effect, variance and covered samples.
        """
        # Initialize rule using the greedy method
        current_rule = self.greedy_initialization()
        self.mu_t = self.calculate_Q(current_rule, "Q1") / self.calculate_Q(
            current_rule, "Q3"
        )
        (
            current_objective_value,
            current_tau,
            current_variance,
        ) = self.objective_function(current_rule)

        for iteration in range(max_iterations):
            # Perform one MM iteration to get a new rule
            print(iteration)
            new_rule = self.mm_iteration(current_rule)
            new_objective_value, new_tau, new_variance = self.objective_function(
                new_rule
            )

            # Check for improvement towards convergence
            if new_objective_value <= current_objective_value + tolerance:
                break  # Exit if no significant improvement

            # Update current rule and objective value if there is an improvement
            current_rule = new_rule
            current_objective_value = new_objective_value
            current_tau = new_tau
            current_variance = new_variance

        final_rule, final_objective_value = self.true_objective_local_search_candidate(
            current_rule
        )
        if final_objective_value > current_objective_value:
            current_rule = final_rule
            current_objective_value = final_objective_value
            current_tau, current_variance = self.objective_function(current_rule)[1:3]

        # final_rule, operation = self.true_objective_local_search_candidate(current_rule)
        # final_objective_value, final_tau, final_variance = self.objective_function(
        #     final_rule
        # )
        # if final_objective_value > current_objective_value:
        #     current_rule = final_rule
        #     current_objective_value = final_objective_value
        #     current_tau = final_tau
        #     current_variance = final_variance

        covered_indices = self.data.index[self.calculate_rule_cover(current_rule)]
        # Return the best rule, its objective function value and the covered samples
        return (
            current_rule,
            current_objective_value,
            current_tau,
            current_variance,
            covered_indices,
        )

    def update_used_features(self):
        """
        Update the set of used features based on the current rules.
        """
        self.used_features = set(
            feature_idx
            for rule in self.rules
            for feature_idx, used in enumerate(rule)
            if used
        )

    def greedy_initialization(self):
        """
        Initialize the rule by greedily adding features that maximize the objective function,
        while ensuring that features not already used in existing rules are considered.
        """
        current_rule = np.zeros(len(self.covariate_cols), dtype=int)
        best_objective_value, _, _ = self.objective_function(current_rule)

        while True:
            best_feature = None
            for feature_idx in range(len(self.covariate_cols)):
                # Check if the feature is not used yet and not in the current rule
                if (
                    feature_idx not in self.used_features
                    and current_rule[feature_idx] == 0
                ):
                    # Temporarily add the feature to the rule
                    current_rule[feature_idx] = 1
                    # Calculate the new objective value with this feature added
                    objective_value, _, _ = self.objective_function(current_rule)
                    # Check if this is the best objective value so far
                    if objective_value > best_objective_value:
                        best_objective_value = objective_value
                        best_feature = feature_idx
                    # Remove the feature again to test the next one
                    current_rule[feature_idx] = 0

            # If no feature improves the objective value, stop the loop
            if best_feature is None:
                break

            # Add the best feature to the rule
            current_rule[best_feature] = 1
            self.used_features.add(best_feature)  # Update the set of used features

            # Optional: If there's a max rule length, check if we've reached it
            if (
                self.max_rule_length is not None
                and current_rule.sum() >= self.max_rule_length
            ):
                break

        return current_rule

    def calculate_rule_cover(self, rule):
        """
        Calculate which samples are covered by a given rule.

        :param rule: A binary array representing the presence (1) or absence (0) of conditions in the rule.
        :return: A boolean Series indicating which samples are covered by the rule.
        """
        # If the rule is empty, return a Series that covers all samples
        if not rule.any():
            return pd.Series([True] * self.num_samples)

        # Only consider the features included in the rule (where rule is 1)
        included_features = [
            self.covariate_cols[i] for i in range(len(rule)) if rule[i] == 1
        ]

        # Compute which samples satisfy the rule's conditions
        conditions = self.data[included_features].all(axis=1)

        return conditions

    def calculate_Q1(self, rule):
        """
        计算 Q1
        Calculate the Q1 value for a given rule.

        :param rule: The binary array representing the current rule.
        :return: The Q1 value.
        """
        cover = self.calculate_rule_cover(rule)
        Q1 = self.data.loc[
            cover & self.data[self.treatment_col].astype(bool),
            self.weighted_outcome_col,
        ].sum()
        return Q1

    def calculate_Q2(self, rule):
        """
        计算 Q2
        Calculate the Q2 value for a given rule.

        :param rule: The binary array representing the current rule.
        :return: The Q2 value.
        """
        cover = self.calculate_rule_cover(rule)
        Q2 = self.data.loc[
            cover & ~self.data[self.treatment_col].astype(bool),
            self.weighted_outcome_col,
        ].sum()
        return Q2

    def calculate_Q3(self, rule):
        """
        计算 Q3 - 对于接受干预(T=1)的样本子集的权重的和.

        :param rule: The binary array representing the current rule.
        :return: The Q3 value.
        """
        cover = self.calculate_rule_cover(rule)
        Q3 = self.data.loc[
            cover & self.data[self.treatment_col].astype(bool),
            self.weight_col,
        ].sum()
        return Q3

    def calculate_Q4(self, rule):
        """
        计算 Q4 - 对于未接受干预(T=0)的样本子集的权重的和.

        :param rule: The binary array representing the current rule.
        :return: The Q4 value.
        """
        cover = self.calculate_rule_cover(rule)
        Q4 = self.data.loc[
            cover & ~self.data[self.treatment_col].astype(bool),
            self.weight_col,
        ].sum()
        return Q4

    def calculate_Q5(self, rule):
        """
        Calculate Q5, the weighted sum of squared differences from the previous mean for treated samples covered by the rule.
        :param rule: The binary array representing the current rule.
        :return: The Q5 value.
        """
        if self.mu_t is None:  # 如果 mu_t 未初始化，则抛出错误
            raise ValueError("mu_t has not been initialized.")
        cover = self.calculate_rule_cover(rule)
        # Calculate the weighted sum of squared differences from the stored mu_t
        Q5 = (
            self.data.loc[
                cover & self.data[self.treatment_col].astype(bool), self.weight_col
            ]
            * (
                self.data.loc[
                    cover & self.data[self.treatment_col].astype(bool), self.outcome_col
                ]
                - self.mu_t
            )
            ** 2
        )
        Q5 = Q5.sum()
        return Q5

    def calculate_Q6(self, rule):
        return self.calculate_Q3(rule)

    def calculate_Q(self, rule, Q_type):
        """
        Calculate the Q value for a given rule.

        :param rule: The binary array representing the current rule.
        :param Q_type: A string indicating whether to compute.
        :return: The Q value.
        """
        # Choose the correct function to calculate Q
        calculate_Q_func = (
            self.calculate_Q1
            if Q_type == "Q1"
            else (
                self.calculate_Q2
                if Q_type == "Q2"
                else self.calculate_Q3
                if Q_type == "Q3"
                else self.calculate_Q4
                if Q_type == "Q4"
                else self.calculate_Q5
                if Q_type == "Q5"
                else self.calculate_Q6
            )
        )

        return calculate_Q_func(rule)

    def objective_function(self, rule):
        """
        Objective function based on the causal effect estimation.

        :param rule: The binary array representing the current rule.
        :return: The objective function value.
        """
        # Calculate Q for the given rule
        Q1 = self.calculate_Q(rule, "Q1")
        Q2 = self.calculate_Q(rule, "Q2")
        Q3 = self.calculate_Q(rule, "Q3")
        Q4 = self.calculate_Q(rule, "Q4")
        # Q5 = self.calculate_Q(rule, "Q5")
        # Q6 = self.calculate_Q(rule, "Q6")

        # Avoid division by zero by setting a lower bound on Q3 and Q4
        Q3 = max(Q3, 1e-10)
        Q4 = max(Q4, 1e-10)

        # 用的是 risk difference
        tau_R = (Q1 / Q3) - (Q2 / Q4)
        # variance_R = Q5 / Q6

        cover = self.calculate_rule_cover(rule)
        covered_df = self.data.loc[cover]

        # if not meet the minimum coverage, return -inf
        if covered_df.shape[0] / self.data.shape[0] < self.minimum_coverage:
            return -np.inf, tau_R, 0

        if "TE" in covered_df.columns:
            # If treatment effect column exists, use it to calculate the true treatment effect
            true_treatment_effect = covered_df[self.treatment_effect_col].mean()
        else:
            true_treatment_effect = np.nan # Placeholder for true treatment effect if not available
        # print("estimated tau_R: ", tau_R, "true_tau: ", true_treatment_effect)

        treated_cover_df = self.data.loc[
            cover & self.data[self.treatment_col].astype(bool)
        ]

        weighted_treated_mean_y = Q1 / Q3
        weighted_treated_variance_y = (
            treated_cover_df[self.weight_col]
            * (treated_cover_df[self.outcome_col] - weighted_treated_mean_y) ** 2
        ).sum() / Q3

        # print("estimated variance: ", Q5 / Q6, "true_variance: ", weighted_treated_variance_y)

        # weighted_objective = tau_R
        weighted_objective = tau_R - self.variance_weight * weighted_treated_variance_y

        return weighted_objective, tau_R, weighted_treated_variance_y

    def calculate_submodular_bounds(self, X, Y, Q_type):
        """
        Calculate the submodular bounds m1 and m2 at solution X given candidate Y.

        :param X: Current solution.
        :param Y: Candidate solution.
        :param Q_type: A string indicating whether to compute the bound.
        """
        # Choose the correct function to calculate Q
        calculate_Q_func = (
            self.calculate_Q1
            if Q_type == "Q1"
            else self.calculate_Q4
            if Q_type == "Q4"
            else self.calculate_Q6
        )  # For Q6, use calculate_Q6

        # Initialize the values of m1 and m2 using the chosen function passing X as parameter
        m1 = m2 = calculate_Q_func(X)

        # Compute m1 and m2 using the definitions provided
        for j in range(len(X)):
            if X[j] == 1 and Y[j] == 0:
                # Feature present in X but not in Y (X \ {j})
                marginal_gain_X_without_j = self.marginal_gain(X, j, Q_type, add=False)
                m1 -= marginal_gain_X_without_j  # Subtract from m1 using first bound

                # For m2, need to compute marginal gain with respect to V without {j}
                full_set = np.ones_like(X)
                marginal_gain_full_set_with_j = self.marginal_gain(
                    full_set, j, Q_type, add=False
                )
                m2 -= (
                    marginal_gain_full_set_with_j  # Subtract from m2 using second bound
                )

            elif X[j] == 0 and Y[j] == 1:
                # Feature not in X but present in Y (Y \ {j})
                marginal_gain_empty_set = self.marginal_gain(
                    np.zeros_like(X), j, Q_type, add=True
                )
                marginal_gain_X = self.marginal_gain(X, j, Q_type, add=True)
                m1 += marginal_gain_empty_set  # Add to m1 using first bound
                m2 += marginal_gain_X  # Add to m2 using second bound

        return m1, m2

    def surrogate_objective(self, current_solution, candidate_solution):
        """
        Calculate a surrogate lower bound of the objective function for use in the MM algorithm.
        """
        # If the candidate solution covers no samples, return -inf
        if self.calculate_rule_cover(candidate_solution).sum() == 0:
            return -np.inf

        # Calculate the submodular bounds for Q1 using the current solution
        m1_current_Q1, m2_current_Q1 = self.calculate_submodular_bounds(
            current_solution, candidate_solution, Q_type="Q1"
        )
        # Calculate the submodular bounds for Q4 using the current solution
        m1_current_Q4, m2_current_Q4 = self.calculate_submodular_bounds(
            current_solution, candidate_solution, Q_type="Q4"
        )
        # Calculate the submodular bounds for Q6 using the current solution
        m1_current_Q6, m2_current_Q6 = self.calculate_submodular_bounds(
            current_solution, candidate_solution, Q_type="Q6"
        )

        # Choose the larger of the two bounds
        m_chosen_Q1 = max(m1_current_Q1, m2_current_Q1)
        m_chosen_Q4 = max(m1_current_Q4, m2_current_Q4)
        m_chosen_Q6 = max(m1_current_Q6, m2_current_Q6)

        # Calculate Q2, Q3, Q5 for the current and candidate solutions
        Q2_current = self.calculate_Q2(current_solution)
        Q2_candidate = self.calculate_Q2(candidate_solution)
        Q3_current = self.calculate_Q3(current_solution)
        Q3_candidate = self.calculate_Q3(candidate_solution)
        Q5_current = self.calculate_Q5(current_solution)
        Q5_candidate = self.calculate_Q5(candidate_solution)

        # Avoid taking log of zero by setting a lower bound on Q2
        Q2_current = max(Q2_current, 1e-10)
        # Q2_candidate = max(Q2_candidate, 1e-10)
        Q3_current = max(Q3_current, 1e-10)
        # Q3_candidate = max(Q3_candidate, 1e-10)
        Q5_current = max(Q5_current, 1e-10)

        taylor_approximation_Q2 = (
            np.log(Q2_current) + (Q2_candidate - Q2_current) / Q2_current
        )
        taylor_approximation_Q3 = (
            np.log(Q3_current) + (Q3_candidate - Q3_current) / Q3_current
        )
        taylor_approximation_Q5 = (
            np.log(Q5_current) + (Q5_candidate - Q5_current) / Q5_current
        )

        # Avoid taking log of zero by using a small positive value or setting a default value
        if m_chosen_Q1 <= 0:
            # log_m_chosen_Q1 = -1
            log_m_chosen_Q1 = -np.inf
        else:
            log_m_chosen_Q1 = np.log(m_chosen_Q1)

        if m_chosen_Q4 <= 0:
            # log_m_chosen_Q4 = -1
            log_m_chosen_Q4 = -np.inf
        else:
            log_m_chosen_Q4 = np.log(m_chosen_Q4)

        if m_chosen_Q6 <= 0:
            # log_m_chosen_Q4 = -1
            log_m_chosen_Q6 = -np.inf
        else:
            log_m_chosen_Q6 = np.log(m_chosen_Q6)

        # if 1 - self.variance_weight == 0:
        #     effect_component = 0
        # else:
        #     effect_component = (1 - self.variance_weight) * (
        #         log_m_chosen_Q1
        #         + log_m_chosen_Q4
        #         - taylor_approximation_Q2
        #         - taylor_approximation_Q3
        #     )

        effect_component = (
            log_m_chosen_Q1
            + log_m_chosen_Q4
            - taylor_approximation_Q2
            - taylor_approximation_Q3
        )

        if self.variance_weight == 0:
            variance_component = 0
        else:
            variance_component = self.variance_weight * (
                log_m_chosen_Q6 - taylor_approximation_Q5
            )

        result = effect_component + variance_component

        return result

    def set_treated_weighted_outcome_to_small_number(self, covered_indices):
        """
        Set the weighted outcome(v) of covered T=1 samples to a small positive number.

        :param covered_indices: Indexes of the samples to be updated.
        """
        small_positive_number = 1e-6
        covered_mask = (self.data.index.isin(covered_indices)) & (
            self.data[self.treatment_col] == 1
        )
        self.data.loc[covered_mask, self.weighted_outcome_col] = small_positive_number

    def find_rules(self, num_rules_to_find=1):
        """
        Find a given number of rules, setting the weighted outcome of covered T=1 samples to a small positive number.

        :param num_rules_to_find: The number of rules to find.
        """
        for i in range(num_rules_to_find):
            print("Finding rule #", i + 1)
            rule, objective_value, tau, variance, covered_indices = self.optimize()

            # If no more significant rules are found, stop
            if not np.any(rule):
                break

            self.rules.append(rule)
            self.update_used_features()
            self.objective_values.append(objective_value)
            self.treatment_effect_values.append(tau)
            self.variance_values.append(variance)
            self.covered_data_indices.append(covered_indices)

            # Set the weighted outcome of covered T=1 samples to a small positive number
            self.set_treated_weighted_outcome_to_small_number(covered_indices)

    def print_rules(self):
        """
        Print out all the learned rules with their related values.
        """
        for idx, rule in enumerate(self.rules):
            # Convert rule from binary array to conditions
            conditions = [
                self.covariate_cols[i] for i, included in enumerate(rule) if included
            ]
            conditions_str = " and ".join(conditions)
            self.conditions_str.append(conditions_str)
            print(f"Rule {idx + 1}: IF {conditions_str}")
            print(f"Objective Value: {self.objective_values[idx]}")
            print(f"Treatment Effect: {self.treatment_effect_values[idx]}")
            print(f"Treated Variance: {self.variance_values[idx]}")
            print(f"Covered Data Count: {len(self.covered_data_indices[idx])}")
            # print(f"Covered Data Indices: {self.covered_data_indices[idx]}")
            print("-" * 40)
