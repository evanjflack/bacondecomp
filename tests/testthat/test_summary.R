# This tests that the summary of average estimate by type is the same as the 
# stata output using the castle data set

two_by_twos <- bacon(l_homicide ~ post,
                  data = bacondecomp::castle,
                  id_var = "state",
                  time_var = "year", 
                  quietly = T)

sum_df <- print_summary(two_by_twos, return_df = TRUE)
sum_df <- sum_df[order(sum_df$type), ]
sum_df$weight <- round(sum_df$weight, 3)
sum_df$avg_est <- round(sum_df$avg_est, 3)

exp_sum_df <- data.frame(type = c("Earlier vs Later Treated", 
                                  "Later vs Earlier Treated", 
                                  "Treated vs Untreated"), 
                         weight = c(0.060, 0.032, 0.908), 
                         avg_est = c(-0.006, 0.070, 0.088))
exp_sum_df <- exp_sum_df[order(exp_sum_df$type), ]

test_that("Print summary is correct", {
  expect_equal(mean(exp_sum_df == sum_df), 1)
})
