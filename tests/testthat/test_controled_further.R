df_test_orig <- bacondecomp::castle 
df_test_orig <- df_test_orig[, c("l_homicide",
                                 "l_pop",
                                 "l_income",
                                 "post",
                                 "state",
                                 "year")]

compare_fe_and_bacon <- function(formula){
  env_formula <- formula(terms(formula, data = df_test_orig))
  two_way_fe <- lm(data = df_test_orig,
                   formula = update(env_formula, . ~ . + factor(state) + factor(year)))
  two_way_fe_coef <- two_way_fe$coefficients["post"]
  names(two_way_fe_coef) <- NULL
  bacon_result <- bacon(formula = formula,
                        data = df_test_orig,
                        id_var = "state",
                        time_var = "year")
  bacon_weighted_mean <- weighted.mean(bacon_result$two_by_twos$estimate,
                                       bacon_result$two_by_twos$weight)
  bacon_fe_coef <- bacon_weighted_mean * (1 - bacon_result$Omega) + bacon_result$beta_hat_w * bacon_result$Omega
  
  return(list(bacon = bacon_fe_coef,
              two_way = two_way_fe_coef))
}


simple_formula <- l_homicide ~ post + l_pop + l_income  
simple_formula_results <- compare_fe_and_bacon(simple_formula)

interaction_formula <- l_homicide ~ post + l_pop*l_income
interaction_formula_results <- compare_fe_and_bacon(interaction_formula)


dot_formula <- l_homicide ~ post + . - state - year
dot_formula_results <- compare_fe_and_bacon(dot_formula) 

test_that("Two Way FE recovered", {
  expect_equal(simple_formula_results$bacon, simple_formula_results$two_way)
  expect_equal(interaction_formula_results$bacon, interaction_formula_results$two_way)
  expect_equal(dot_formula_results$bacon, dot_formula_results$two_way)
  expect_equal(dot_formula_results$bacon, simple_formula_results$two_way)
})
