# Formulas
variables_unpacked <- unpack_variable_names(l_homicide ~ post + l_pop + l_income)

variables_unpacked
test_that("Variables Unpacked Correctly", {
  expect_equal(variables_unpacked$outcome_var, "l_homicide")
  expect_equal(variables_unpacked$treated_var, "post")
  expect_equal(variables_unpacked$control_vars, c("l_pop", "l_income"))
})
