library(bacon)

df_test <- data.frame(time = rep(seq(1, 100), 3),
                      id = c(rep(1, 100), rep(2, 100), rep(3, 100)),
                      treated = c(rep(0, 100), rep(0, 33), rep(1, 67),
                                  rep(0, 84), rep(1, 16)),
                      y = rnorm(300), 
                      x1 = rcauchy(300))

control_vars <- "x1"

control_formula <- paste(control_vars, collapse = " + ")
control_formula <- as.formula(paste("treated ~", control_formula))

data <- calculate_ds(df_test, control_formula)
Sigma <- calculate_Sigma(data)
one_minus_Sigma <- calculate_one_minus_Sigma(data)

