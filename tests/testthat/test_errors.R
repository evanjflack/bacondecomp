library(bacon)

# These tests check that the correct ewrror message is returned

# Unbalanced Panel
castle <- bacon::castle
castle <- castle[-sample(1:nrow(castle), 1), ]

test_that("Unbalanced panel error", {
  expect_error(bacon(l_homicide ~ post,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), "Unbalanced Panel")
})

# NA Observations
castle <- bacon::castle
castle[-sample(1:nrow(castle), 1), "post"] <- NA

test_that("NA error", {
  expect_error(bacon(l_homicide ~ post,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), "NA observations")
})
