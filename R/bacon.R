#' Goodman-Bacon decomposition
#'
#' bacon() is a function that perfroms the Goodman-Bacon decomposition for
#'  differences-in-differences with variation in treatment timing.
#'
#' @param formula a symbolic representation of the
#'  model to be fitted.
#' @param data a data.frame containing the variables in the model.
#' @param id_var character, the name of id variable for units.
#' @param time_var character, the name of time variable.
#' @param quiet logical, if TRUE does not print overal two-way fixed effects
#'  estimate or summary of 2x2 estimates by type.
#'
#' @return data.frame of all 2x2 estimates and weights
#'
#' @examples
#' # Math Reform ---------------------------------------------------------------
#' df_bacon <- bacon(incearn_ln ~ reform_math,
#'                   data = bacon::math_reform,
#'                   id_var = "state",
#'                   time_var = "class")
#' ggplot(df_bacon) +
#'   aes(x = weight, y = estimate, shape = factor(type)) +
#'   labs(x = "Weight", y = "Estimate", shape = "Type") +
#'   geom_point()
#'
#' # Castle Doctrine -----------------------------------------------------------
#' df_bacon <- bacon(l_homicide ~ post,
#'                   data = bacon::castle,
#'                   id_var = "state",
#'                   time_var = "year")
#' ggplot(df_bacon) +
#'   aes(x = weight, y = estimate, shape = factor(type)) +
#'   labs(x = "Weight", y = "Estimate", shape = "Type") +
#'   geom_point()
#'
#' @export
bacon <- function(formula,
                  data,
                  id_var,
                  time_var,
                  quiet = FALSE) {
  # Rename variables
  outcome_var <- as.character(formula)[2]
  treated_var <- as.character(formula)[3]
  data <- data %>%
    rename("id" = id_var, "time" = time_var, "treated" = treated_var,
           "outcome" = outcome_var) %>%
    select(id, time, treated, outcome)

  # Check for NA observations
  nas <- sum(is.na(data))
  if (nas > 0) stop("NA observations")

  # Check for balanced panel
  bal <- data %>%
    group_by(id) %>%
    tally()
  balanced <- ifelse(mean(bal$n == bal$n[1]) == 1, 1, 0)
  if(!balanced) stop("Unbalanced Panel")

  df_treat <- data %>%
    group_by(id) %>%
    filter(treated == 1) %>%
    filter(row_number() ==1) %>%
    select(id, time) %>%
    rename("treat_time" = "time")
  data <- data %>%
    merge(df_treat, by = "id", all.x = T) %>%
    arrange(id, time) %>%
    mutate(treat_time = ifelse(is.na(treat_time), 99999, treat_time))

  # First period in the panel
  first_period <- min(data$time)

  # create data.frame of all posible 2x2 estimates
  two_by_twos <- expand.grid(unique(data$treat_time), unique(data$treat_time)) %>%
    rename("treated" = "Var1", "untreated" = "Var2") %>%
    subset(!(treated == untreated)) %>%
    subset(!(treated == 99999)) %>%
    subset(!(treated == first_period)) %>%
    mutate(estimate = 0, weight = 0)

  for (i in 1:nrow(two_by_twos)) {
    treated_group <- two_by_twos[i, "treated"]
    untreated_group <- two_by_twos[i, "untreated"]
    data1 <- subset(data, treat_time %in% c(treated_group, untreated_group))

    # Calculated weight
    # n_u - observations in untreated group
    # n_t - observations in treated group
    # p_t - proportion of the time treated group was treated
    # n_e - observations in early treated group
    # n_l - observations in late treated group
    # p_e - proportion of the time early treated group was treated
    # p_l - proportion of the time late treated group was treated
    if (untreated_group == 99999) {
      # Treated vs untreated
      n_u <- sum(data1$treat_time == untreated_group)
      n_t <- sum(data1$treat_time == treated_group)
      p_t <- mean(data1[data1$treat_time == treated_group, "treated"])
      weight1 <- n_u * n_t * p_t * (1 - p_t)
    } else if (treated_group < untreated_group) {
      # early vs late (before late is treated)
      data1 <- subset(data1, time < untreated_group)
      n_e <- sum(data1$treat_time == treated_group)
      n_l <- sum(data1$treat_time == untreated_group)
      p_e <- mean(data1[data1$treat_time == treated_group, "treated"])
      p_l <- mean(data1[data1$treat_time == untreated_group, "treated"])
      weight1 <- n_e * n_l * (p_e - p_l) * (1 - (p_e - p_l))
      weight1 <- weight1 * (1 - p_e) / (1 - p_e + p_l)
    } else if (treated_group > untreated_group) {
      # late vs early (after early is treated)
      data1 <- subset(data1, time >= untreated_group)
      n_e <- sum(data1$treat_time == untreated_group)
      n_l <- sum(data1$treat_time == treated_group)
      p_e <- mean(data1[data1$treat_time == untreated_group, "treated"])
      p_l <- mean(data1[data1$treat_time == treated_group, "treated"])
      weight1 <- n_e * n_l * (p_e - p_l) * (1 - (p_e - p_l))
      weight1 <- weight1 * (p_l / (1 - p_e + p_l))
    }

    # Estimate 2x2 diff-in-diff
    estimate1 <- lm(outcome ~ treated + factor(time) + factor(id),
                   data = data1)$coefficients[2]

    two_by_twos[i, ] <- two_by_twos[i, ] %>%
      mutate(estimate = estimate1, weight = weight1)
  }

  # Rescale weights to sum to 1
  two_by_twos <- two_by_twos %>%
    mutate(weight = weight / sum(weight)) %>%
    # Classify estimate type
    mutate(type = ifelse(untreated == 99999, "Treated vs Untreated",
                         ifelse(untreated == first_period,
                                "Always Treated vs Later Treated",
                                ifelse(treated < untreated, "Early vs Late",
                                       "Late vs Early"))))

  # Print two-way FE estimate and summary of 2x2 estimates by type
  if (quiet == F) {
    overall_est <- weighted.mean(two_by_twos$estimate, two_by_twos$weight)
    print(paste0("Two-way FE estimate = ", overall_est))

    two_by_twos %>%
      group_by(type) %>%
      summarise(weight = sum(weight), avg_estimate = mean(estimate)) %>%
      print()
  }
  return(two_by_twos)
}
