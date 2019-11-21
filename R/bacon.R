#' Goodman-Bacon decomposition
#'
#' bacon() is a function that perfroms the Goodman-Bacon decomposition for
#'  differences-in-differences with variation in treatment timing.
#'
#' @param df a data.frame
#' @param id_var character, name of id variable for units
#' @param time_var character, name of time variable
#' @param treated_var character, name of treatment variable (binary)
#' @param outcome_var character, name of outcome variable
#' @param quiet logical, if TRUE does not print summary of 2x2 estimates by
#'  type
#'
#' @return data.frame of all 2x2 estimates and weights
#'
#' @examples
#' # Math Reform ---------------------------------------------------------------
#' df_bacon <- bacon(df = bacon::math_reform,
#'                   id_var = "state",
#'                   time_var = "class",
#'                   treated_var = "reform_math",
#'                   outcome_var = "incearn_ln")
#' ggplot(df_bacon) +
#'   aes(x = weight, y = estimate, shape = factor(type)) +
#'   labs(x = "Weight", y = "Estimate", shape = "Type") +
#'   geom_point()
#'
#' # Castle Doctrine -----------------------------------------------------------
#' df_bacon <- bacon(df = bacon::castle,
#'                   id_var = "state",
#'                   time_var = "year",
#'                   treated_var = "post",
#'                   outcome_var = "l_homicide")
#' ggplot(df_bacon) +
#'   aes(x = weight, y = estimate, shape = factor(type)) +
#'   labs(x = "Weight", y = "Estimate", shape = "Type") +
#'   geom_point()
#'
#' @export
bacon <- function(df, id_var = "id",
                  time_var = "time",
                  treated_var = "treated",
                  outcome_var = "outcome",
                  quiet = FALSE) {
  # Rename variables
  df <- df %>%
    rename("id" = id_var, "time" = time_var, "treated" = treated_var,
           "outcome" = outcome_var)

  df_treat <- df %>%
    group_by(id) %>%
    filter(treated == 1) %>%
    filter(row_number() ==1) %>%
    select(id, time) %>%
    rename("treat_time" = "time")

  df <- df %>%
    merge(df_treat, by = "id", all.x = T) %>%
    arrange(id, time)

  # Set NAS for treat_time to 99999
  df[is.na(df$treat_time), "treat_time"] <- 99999

  # create data.frame of all posible 2x2 estimates
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
      weight1 <- n_u * n_t * p_t * (1 - p_t)
    } else if (treated_group < untreated_group) {
      df1 <- subset(df1, time < untreated_group)
      n_e <- sum(df1$treat_time == treated_group)
      n_l <- sum(df1$treat_time == untreated_group)
      p_e <- mean(df1[df1$treat_time == treated_group, "treated"])
      p_l <- mean(df1[df1$treat_time == untreated_group, "treated"])
      weight1 <- n_e * n_l * (p_e - p_l) * (1 - (p_e - p_l))
      weight1 <- weight1 * (1 - p_e) / (1 - p_e + p_l)
    } else if (treated_group > untreated_group) {
      df1 <- subset(df1, time >= untreated_group)
      n_e <- sum(df1$treat_time == untreated_group)
      n_l <- sum(df1$treat_time == treated_group)
      p_e <- mean(df1[df1$treat_time == untreated_group, "treated"])
      p_l <- mean(df1[df1$treat_time == treated_group, "treated"])
      weight1 <- n_e * n_l * (p_e - p_l) * (1 - (p_e - p_l))
      weight1 <- weight1 * (p_l / (1 - p_e + p_l))
    }
    estimate1 <- lm(outcome ~ treated + factor(time) + factor(id),
                   data = df1)$coefficients[2]
    two_by_twos[i, ] <- two_by_twos[i, ] %>%
      mutate(estimate = estimate1, weight = weight1)
  }
  two_by_twos %<>%
    mutate(weight = weight / sum(weight)) %>%
    mutate(type = ifelse(untreated == 99999, "Treated vs Unteated",
                         ifelse(treated < untreated, "Early vs Late",
                                "Late vs Early")))
  if (quiet == F) {
    two_by_twos %>%
      group_by(type) %>%
      summarise(weight = sum(weight), avg_estimate = mean(estimate)) %>%
      print()
  }
  return(two_by_twos)
}
