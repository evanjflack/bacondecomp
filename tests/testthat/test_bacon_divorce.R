divorce_df <- bacondecomp::divorce

stevenson_wolfson_replicates <- c(
  1.6,
  -1.5,
  -1.5,
  -3.0,
  -8.0,
  -10.0,
  -11.9,
  -12.8,
  -13.3,
  -16.4,
  -18.7
)
model_fit <- lm(data = divorce_df[divorce_df$sex == 2, ],
                suiciderate_elast_jag ~ factor(chyrspos) + factor(st) + factor(year))

model_coefs <- (model_fit$coefficients)[2:12]
names(model_coefs) <- NULL
model_coefs <- round(model_coefs, 3)
model_coefs <- model_coefs * 100

test_that("Table 1 col 1 Stevenson and Wolfers  Replicates",
          {
            expect_equal(
              model_coefs,
              stevenson_wolfson_replicates
            )
          })