# These tests check that the correct ewrror message is returned

# Unbalanced Panel
castle <- bacondecomp::castle
castle <- castle[-sample(1:nrow(castle), 1), ]

test_that("Unbalanced panel error", {
  expect_error(bacon(l_homicide ~ post,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), "Unbalanced Panel")
})

# NA Observations (Uncontrolled)
castle <- bacondecomp::castle
castle[sample(1:nrow(castle), 1), "post"] <- NA

test_that("NA error uncontrolled", {
  expect_error(bacon(l_homicide ~ post,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), "NA observations")
})

# NA Observations (Controlled)
castle <- bacondecomp::castle
castle[sample(1:nrow(castle), 1), "l_pop"] <- NA

test_that("NA error uncontrolled", {
  expect_error(bacon(l_homicide ~ post + l_income + l_pop,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), "NA observations")
})

# Weakly Increasing Treatment
castle <- bacondecomp::castle
castle[castle$state == "Alabama" & castle$year == 2009, "post"] <- 0
test_that("Weakly Increasing Treatment", {
  expect_error(bacon(l_homicide ~ post,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), 
               "Treatment not weakly increasing with time")
})
