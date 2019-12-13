library(bacon)

# This test calculates the bacon decomposition and then takes the weighted sum
# of the decomposition which should return the traditional two-way FE estimates

# Bacon decomp
df_bacon <- bacon(l_homicide ~ post,
                  data = bacon::castle,
                  id_var = "state",
                  time_var = "year")
# Weighted Sum
two_way_bacon_coef <- sum(df_bacon$estimate * df_bacon$weight)


# Two way FE
two_way_fe <- lm(l_homicide ~ post + factor(state) + factor(year),
                 data = bacon::castle)
two_way_fe_coef <- two_way_fe$coefficients["post"]
names(two_way_fe_coef) <- NULL


test_that("Two Way FE recovered", {
  expect_equal(two_way_bacon_coef, two_way_fe_coef)
})



##### Multivariate #####

ret_bacon <- bacon(l_homicide ~ post + income + police,
                   data = bacon::castle,
                   id_var = "state",
                   time_var = "year")

test_that("Multivariate bacon returns working object", {
  expect_equal(class(ret_bacon), "list")
  expect_true(length(ret_bacon) == 3)
})


# Bacon estimate
beta_hat_w <- ret_bacon$beta_hat_w
beta_hat_b <- weighted.mean(ret_bacon$two_by_twos$estimate, 
                            ret_bacon$two_by_twos$weight)

Sigma <- ret_bacon$Sigma

two_way_bacon_coef_multi <- Sigma*beta_hat_w + (1 - Sigma)*beta_hat_b

# Two way FE estimate
two_way_fe_multi <- lm(l_homicide ~ post + income + police + factor(state) + factor(year),
                 data = bacon::castle)
two_way_fe_coef_multi <- two_way_fe_multi$coefficients["post"]
names(two_way_fe_coef_multi) <- NULL

test_that("Two Way FE multi recovered", {
  expect_equal(two_way_bacon_coef_multi, two_way_fe_coef_multi)
})

