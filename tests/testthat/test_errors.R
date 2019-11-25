library(bacon)

# This test calculates the bacon decomposition and then takes the weighted sum
# of the decomposition which should return the traditional two-way FE estimates

castle <- bacon::castle
castle <- castle[-sample(1:nrow(castle), 1), ]

test_that("Unbalanced panel error", {
  expect_error(bacon(l_homicide ~ post,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), "Unbalanced Panel")
})

castle <- bacon::castle
castle[-sample(1:nrow(castle), 1), "post"] <- NA

test_that("NA error", {
  expect_error(bacon(l_homicide ~ post,
                     data = castle,
                     id_var = "state",
                     time_var = "year"), "NA observations")
})



