library(bacon)

# Create fake data
df_test_orig <- bacon::castle
formula <- l_homicide ~ post + l_pop + l_income
id_var <- "state"
time_var <- "year"

vars <- unpack_variable_names(formula)
outcome_var <- vars$outcome_var
treated_var <- vars$treated_var
control_vars <- vars$control_vars

df_test <- rename_vars(df_test_orig, id_var, time_var, outcome_var, treated_var)

# Create grid of treatment groups
treatment_group_calc <- create_treatment_groups(df_test, control_vars, 
                                                return_merged_df = TRUE)
# extract data to test
test_data <- treatment_group_calc$data
# dummy control var
control_formula <- update(formula,
                          paste0("treated ~ .  + factor(time) + factor(id) -",
                                 treated_var))
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

# Beta Hats
beta_hat_dd_23 <- cov(test_data$outcome, test_data$d_it_til)/var(test_data$d_it_til)

beta_hat_w <- calculate_beta_hat_w(test_data)
beta_hat_b <- calculate_beta_hat_b(test_data)
beta_hat_dd_25 <- Sigma * beta_hat_w + one_minus_Sigma * beta_hat_b

two_way_model <- lm(l_homicide ~ 0 + post + l_pop + l_income + factor(state) + 
                      factor(year),
                    data = df_test_orig)

beta_hat_dd_two_way <- two_way_model$coefficients["post"]
names(beta_hat_dd_two_way) <- NULL

ret_bacon <- bacon(formula,
                   data = bacon::castle,
                   id_var = "state",
                   time_var = "year")

two_by_twos <- ret_bacon$two_by_twos
beta_hat_b_two_by_twos <- sum(two_by_twos$estimate*two_by_twos$weight)

test_that("Recover Two Way Coef Using EQN 23, pg 25", {
  expect_equal(beta_hat_dd_23, beta_hat_dd_two_way)
})

test_that("Recover Two Way Coef Using EQN 25, pg 26 WP", {
  expect_equal(beta_hat_dd_25, beta_hat_dd_two_way)
})

# This will never work - the weighted mean of beta_hat_2x2 doesn't equal
# beta_hat_b. We need to include other terms
test_that("Equation 26 Correct", {
  expect_equal(beta_hat_b, beta_hat_b_two_by_twos)
})


beta_hat_dd_final <- Sigma * beta_hat_w + one_minus_Sigma * sum(two_by_twos$estimate * two_by_twos$weight)

test_that("Equation 27 Correct", {
  expect_equal(beta_hat_dd_final, beta_hat_dd_two_way)
})
