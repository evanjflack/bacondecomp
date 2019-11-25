library(bacon)

# This test calculates the bacon decomposition and then takes the weighted sum
# of the decomposition which should return the traditional two-way FE estimates


# Bacon decomp
decomp <- bacon(l_homicide ~ post,
                data = bacon::castle,
                id_var = "state",
                time_var = "year")
df_bacon <- decomp$decomp
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


