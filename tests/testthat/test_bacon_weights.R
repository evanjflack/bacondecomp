# This test mimics the setup of figure 2 and tests that calculate_weights()
# computes the correct weights

# Weights from the paper
exp_weights <- c(.365, .222, .278, .135)

# Data frame that mimics set up from figure 2, but with random outcomes
df_test <- data.frame(time = rep(seq(1, 100), 3),
                      id = c(rep(1, 100), rep(2, 100), rep(3, 100)),
                      treated = c(rep(0, 100), rep(0, 33), rep(1, 67),
                                  rep(0, 84), rep(1, 16)),
                      y = rnorm(300))

# Find weights
df_bacon <- bacon(y ~ treated,
                  data = df_test,
                  id_var = "id",
                  time_var = "time")

weights <- df_bacon[order(-df_bacon$untreated, df_bacon$treat), "weight"]
weights <- round(weights, 3)


test_that("Example weights recovered", {
  expect_equal(weights, exp_weights)
})
