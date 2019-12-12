library(bacon)
# Create fake data
df_test <- data.frame(time = rep(seq(1, 100), 3),
                      id = c(rep(1, 100), rep(2, 100), rep(3, 100)),
                      treated = c(rep(0, 100), rep(0, 33), rep(1, 67),
                                  rep(0, 84), rep(1, 16)),
                      y = rnorm(300), 
                      x1 = rcauchy(300))
# Create grid of treatment groups
treatment_group_calc <- create_treatment_groups(df_test,
                                                return_merged_df = TRUE)
# extract data to test
test_data <- treatment_group_calc$data
# dummy control var
control_vars <- "x1"
# Formulas
control_formula <- paste(control_vars, collapse = " + ")
control_formula <- as.formula(paste("treated ~", control_formula))
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
