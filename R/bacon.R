#' Decompose two-way fixed effects models
#'
#' bacon() is a function that perfroms the Goodman-Bacon decomposition for
#'  differences-in-differences with variation in treatment timing.
#'
#' @param df a data.frame
#' @param id_var character, name of id variable for units
#' @param time_var character, name of time variable
#' @param treat_time_var character, name of variable indicating when treatment
#'  switched on for unit. Set to NA if unit is never treated.
#' @param treated_var character, name of treatment variable (binary)
#' @param outcome_var character, name of outcome variable
#'
#' @return data.frame of all 2x2 estimates and weights
#'
#' @examples
#' math_reform <- bacon::math_reform
#' math_reform[is.na(math_reform$reformyr_math), "reformyr_math"] <- 99999
#'
#' bacon(math_reform, id_var = "state", time_var = "class",
#'       treat_time_var = "reformyr_math",  treated_var = "reform_math",
#'      outcome_var = "incearn_ln")
#'
#' @export
bacon <- function(df, id_var = "id", time_var = "time",
                  treat_time_var = "treat_time", treated_var = "treated",
                  outcome_var = "outcome") {
  df <- df %>%
    rename("id" = id_var, "time" = time_var, "treat_time" = treat_time_var,
           "treated" = treated_var, "outcome" = outcome_var)

  two_by_twos <- expand.grid(unique(df$treat_time), unique(df$treat_time)) %>%
    rename("treated" = "Var1", "untreated" = "Var2") %>%
    subset(!(treated == untreated | treated == 99999)) %>%
    mutate(estimate = 0, weight = 0)

  for (i in 1:nrow(two_by_twos)) {
    treated_group <- two_by_twos[i, "treated"]
    untreated_group <- two_by_twos[i, "untreated"]

    df1 <- subset(df, treat_time %in% c(treated_group, untreated_group))

    if (untreated_group == 99999) {
      n_u <- sum(df1$treat_time == untreated_group)
      n_t <- sum(df1$treat_time == treated_group)
      p_t <- mean(df1[df1$treat_time == treated_group, "treated"])
      weight1 <- n_u*n_t*p_t*(1-p_t)
    } else if (treated_group < untreated_group) {
      df1 <- subset(df1, time < untreated_group)
      n_e <- sum(df1$treat_time == treated_group)
      n_l <- sum(df1$treat_time == untreated_group)
      p_e <- mean(df1[df1$treat_time == treated_group, "treated"])
      p_l <- mean(df1[df1$treat_time == untreated_group, "treated"])
      weight1 <- n_e*n_l*(p_e - p_l)*(1 - (p_e - p_l))*((1 - p_e)/(1 - p_e + p_l))
    } else if (treated_group > untreated_group) {
      df1 <- subset(df1, time >= untreated_group)
      n_e <- sum(df1$treat_time == untreated_group)
      n_l <- sum(df1$treat_time == treated_group)
      p_e <- mean(df1[df1$treat_time == untreated_group, "treated"])
      p_l <- mean(df1[df1$treat_time == treated_group, "treated"])
      weight1 <- n_e*n_l*(p_e - p_l)*(1 - (p_e - p_l))*(p_l/(1 - p_e + p_l))
    }

    estimate1 <- lm(outcome ~ treated + factor(time) + factor(id),
                   data = df1)$coefficients[2]

    two_by_twos[i, ] <- two_by_twos[i, ] %>%
      mutate(estimate = estimate1, weight = weight1)
  }

  two_by_twos %>%
    mutate(weight = weight/sum(weight))
}
