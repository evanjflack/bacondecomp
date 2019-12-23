library(bacon)

# These tests calculates the bacon decomposition and then takes the weighted sum
# of the decomposition which should return the traditional two-way FE estimates

# Uncontrolled -----------------------------------------------------------------
# Bacon decomp
df_bacon <- bacon(l_homicide ~ post,
                  data = bacon::castle,
                  id_var = "state",
                  time_var = "year")

two_way_bacon_coef <- sum(df_bacon$estimate * df_bacon$weight)

# Two way FE
two_way_fe <- lm(l_homicide ~ post + factor(state) + factor(year),
                 data = bacon::castle)
two_way_fe_coef <- two_way_fe$coefficients["post"]
names(two_way_fe_coef) <- NULL

test_that("Two-way FE recovered", {
  expect_equal(two_way_bacon_coef, two_way_fe_coef)
})

# Controlled -------------------------------------------------------------------
ret_bacon <- bacon(l_homicide ~ post + l_pop + l_income,
                   data = bacon::castle,
                   id_var = "state",
                   time_var = "year")

hi <- ret_bacon$two_by_twos

library(ggplot2)
ggplot(hi) + 
  aes(x = weight, y = estimate) + 
  geom_point()

test_that("Controlled bacon returns working object", {
  expect_equal(class(ret_bacon), "list")
  expect_true(length(ret_bacon) == 3)
})

# Bacon estimate
beta_hat_w <- ret_bacon$beta_hat_w
beta_hat_b <- weighted.mean(ret_bacon$two_by_twos$estimate, 
                            ret_bacon$two_by_twos$weight)
Sigma <- ret_bacon$Sigma
two_way_bacon_coef_cont <- Sigma*beta_hat_w + (1 - Sigma)*beta_hat_b

# Two way FE estimate
two_way_fe_cont <- lm(l_homicide ~ post + l_pop + l_income + factor(state) + 
                        factor(year), data = bacon::castle)
two_way_fe_coef_cont <- two_way_fe_cont$coefficients["post"]
names(two_way_fe_coef_cont) <- NULL

test_that("Two Way FE multi recovered", {
  expect_equal(two_way_bacon_coef_cont, two_way_fe_coef_cont)
})
