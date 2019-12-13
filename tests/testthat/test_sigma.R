library(bacon)
# Create fake data
df_test <- bacon::castle
formula <- l_homicide ~ post + l_pop + l_income
id_var <- "state"
time_var <- "year"

vars <- unpack_variable_names(formula)
outcome_var <- vars$outcome_var
treated_var <- vars$treated_var
control_vars <- vars$control_vars

df_test <- rename_vars(df_test, id_var, time_var, outcome_var, treated_var)

# Create grid of treatment groups
treatment_group_calc <- create_treatment_groups(df_test,
                                                return_merged_df = TRUE)
# extract data to test
test_data <- treatment_group_calc$data
# dummy control var
control_formula <- update(formula,
                          paste0("treated ~ . -", treated_var))
# Calculate the Ds
test_data <- calculate_ds(test_data, control_formula)

# Calculate Sigma
Sigma <- calculate_Sigma(test_data)
one_minus_Sigma <- calculate_one_minus_Sigma(test_data)

# Defined in Eqn 25, page 26 of Working Paper
# https://cdn.vanderbilt.edu/vu-my/wp-content/uploads/sites/2318/2019/07/29170757/ddtiming_7_29_2019.pdf
test_that("Sigma and One Minus Sigma Sum to 1", {
  expect_equal(Sigma + one_minus_Sigma, 1)
})


beta_hat_w <- calculate_beta_hat_w(test_data)
beta_hat_d <- calculate_beta_hat_d(test_data)
beta_hat_dd <- Sigma * beta_hat_w + one_minus_Sigma * beta_hat_d 



beta_hat_dd_two_way <- lm(l_homicide ~ 0 + post + l_pop + l_income + factor(state) + factor(year),
                          data = bacon::castle)

two_way_coef <- beta_hat_dd_two_way$coefficients["post"]
names(two_way_coef) <- NULL



test_that("Recover Two Way Coef Using EQN 25, pg 26 WP", {
  expect_equal(beta_hat_dd, two_way_coef)
})
