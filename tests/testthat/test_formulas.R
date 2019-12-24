# Testing with "+ ." 
data <- castle[, c("state", "year", "l_homicide", "post", "l_pop", "l_income")]
formula <- l_homicide ~ post + . - state - year
id_var <- "state"
time_var <- "year"

formula <- formula(terms(formula, data = data))

vars <- unpack_variable_names(formula)
outcome_var <- vars$outcome_var
treated_var <- vars$treated_var
control_vars <- vars$control_vars
data <- rename_vars(data, id_var, time_var, outcome_var, treated_var)

control_formula <- update(formula,
                          paste0("treated ~ . + factor(time) + factor(id) -", treated_var)
)

test_that("Controll formula with .", {
  expect_equal(control_formula, treated ~ l_pop + l_income + factor(time) + factor(id))
})

# Testing with polynomials
data <- castle[, c("state", "year", "l_homicide", "post", "l_pop", "l_income")]
formula <- l_homicide ~ post + poly(l_pop, 2)
id_var <- "state"
time_var <- "year"

formula <- formula(terms(formula, data = data))

vars <- unpack_variable_names(formula)
outcome_var <- vars$outcome_var
treated_var <- vars$treated_var
control_vars <- vars$control_vars
data <- rename_vars(data, id_var, time_var, outcome_var, treated_var)

control_formula <- update(formula,
                          paste0("treated ~ . + factor(time) + factor(id) -", treated_var)
)

test_that("Controll formula with polynomial", {
  expect_equal(control_formula, treated ~ poly(l_pop, 2) + factor(time) + factor(id))
})

# Testing with interactions
data <- castle[, c("state", "year", "l_homicide", "post", "l_pop", "l_income")]
formula <- l_homicide ~ post + l_pop*l_income
id_var <- "state"
time_var <- "year"

formula <- formula(terms(formula, data = data))

vars <- unpack_variable_names(formula)
outcome_var <- vars$outcome_var
treated_var <- vars$treated_var
control_vars <- vars$control_vars
data <- rename_vars(data, id_var, time_var, outcome_var, treated_var)

control_formula <- update(formula,
                          paste0("treated ~ . + factor(time) + factor(id) -", treated_var)
)

test_that("Controll formula with interactions", {
  expect_equal(control_formula, 
               treated ~ l_pop + l_income + factor(time) + factor(id) + 
                 l_pop:l_income)
})
