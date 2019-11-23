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
#'
#' @author Evan Flack
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
#'
#' \donttest{
#'
#' ggplot(df_bacon) +
#'   aes(x = weight, y = estimate, shape = factor(type)) +
#'   labs(x = "Weight", y = "Estimate", shape = "Type") +
#'   geom_point()
#'
#'   }
#'
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

  data <- data[, c(id_var, time_var, outcome_var, treated_var)]
  colnames(data) <- c("id", "time", "outcome", "treated")
  # Check for NA observations
  nas <- sum(is.na(data))
  if (nas > 0) stop("NA observations")

  # Check for balanced panel
  bal <- aggregate(time ~ id,  data = data, FUN = length)
  balanced <- ifelse(mean(bal$time == bal$time[1]) == 1, 1, 0)
  if (!balanced) stop("Unbalanced Panel")


  df_treat <- data[data$treated == 1, ]
  df_treat <- df_treat[, c("id", "time")]
  df_treat <- aggregate(time ~ id,  data = df_treat, FUN = min)
  colnames(df_treat) <- c("id", "treat_time")
  data <- merge(data, df_treat, by = "id", all.x = T)
  data[is.na(data$treat_time), "treat_time"] <- 99999

  # First period in the panel
  first_period <- min(data$time)

  # create data.frame of all posible 2x2 estimates
  two_by_twos <- expand.grid(unique(data$treat_time),
                             unique(data$treat_time))
  colnames(two_by_twos) <- c("treated", "untreated")
  two_by_twos <- two_by_twos[!(two_by_twos$treated == two_by_twos$untreated), ]
  two_by_twos <- two_by_twos[!(two_by_twos$treated == 99999), ]
  two_by_twos <- two_by_twos[!(two_by_twos$treated == first_period), ]
  two_by_twos[, c("estimate", "weight")] <- 0

  for (i in 1:nrow(two_by_twos)) {
    treated_group <- two_by_twos[i, "treated"]
    untreated_group <- two_by_twos[i, "untreated"]
    data1 <- data[data$treat_time %in% c(treated_group, untreated_group), ]

    weight1 <- calculate_weights(data = data1,
                                 treated_group = treated_group,
                                 untreated_group = untreated_group)

    # Estimate 2x2 diff-in-diff
    estimate1 <- lm(outcome ~ treated + factor(time) + factor(id),
                   data = data1)$coefficients[2]

    two_by_twos[i, "estimate"] <- estimate1
    two_by_twos[i, "weight"] <- weight1
  }

  # Rescale weights to 1
  total_weight <- sum(two_by_twos$weight)
  two_by_twos$weight <- two_by_twos$weight/total_weight
  # Classify estimate type
  two_by_twos$type <- ifelse(two_by_twos$untreated == 99999,
                             "Treated vs Untreated",
                             ifelse(two_by_twos$untreated == first_period,
                                    "Always Treated vs Later Treated",
                                    ifelse(two_by_twos$treated <
                                             two_by_twos$untreated,
                                           "Early vs Late", "Late vs Early")))

  # Print two-way FE estimate and summary of 2x2 estimates by type
  if (quiet == F) {
    # Print two way FE estimate
    overall_est <- weighted.mean(two_by_twos$estimate, two_by_twos$weight)
    print(paste0("Two-way FE estimate = ", overall_est))

    # print summary of 2x2 estimates by type
    avg_est <- aggregate(estimate ~ type, data = two_by_twos, FUN = mean)
    colnames(avg_est) <- c("type", "avg_estimate")
    sum_weight <- aggregate(weight ~ type, data = two_by_twos, FUN = sum)
    avg_est_weight <- merge(avg_est, sum_weight, by = "type")
    print(avg_est_weight)
  }
  return(two_by_twos)
}



#' Calculate Weights for 2x2 Grid
#'
#'  Calculated weights using:
#'  n_u - observations in untreated group
#'  n_t - observations in treated group
#'  p_t - proportion of the time treated group was treated
#'  n_e - observations in early treated group
#'  n_l - observations in late treated group
#'  p_e - proportion of the time early treated group was treated
#'  p_l - proportion of the time late treated group was treated
#'
#'
#' @param data
#' @param treated_group
#' @param untreated_group
#'
#' @return Scalar weight for 2x2 grid
calculate_weights <- function(data,
                              treated_group,
                              untreated_group){
  if (untreated_group == 99999) {
    # Treated vs untreated
    n_u <- sum(data$treat_time == untreated_group)
    n_t <- sum(data$treat_time == treated_group)
    p_t <- mean(data[data$treat_time == treated_group, "treated"])
    weight1 <- n_u * n_t * p_t * (1 - p_t)
  } else if (treated_group < untreated_group) {
    # early vs late (before late is treated)
    data <- subset(data, time < untreated_group)
    n_e <- sum(data$treat_time == treated_group)
    n_l <- sum(data$treat_time == untreated_group)
    p_e <- mean(data[data$treat_time == treated_group, "treated"])
    p_l <- mean(data[data$treat_time == untreated_group, "treated"])
    weight1 <- n_e * n_l * (p_e - p_l) * (1 - (p_e - p_l))
    weight1 <- weight1 * (1 - p_e) / (1 - p_e + p_l)
  } else if (treated_group > untreated_group) {
    # late vs early (after early is treated)
    data <- subset(data, time >= untreated_group)
    n_e <- sum(data$treat_time == untreated_group)
    n_l <- sum(data$treat_time == treated_group)
    p_e <- mean(data[data$treat_time == untreated_group, "treated"])
    p_l <- mean(data[data$treat_time == treated_group, "treated"])
    weight1 <- n_e * n_l * (p_e - p_l) * (1 - (p_e - p_l))
    weight1 <- weight1 * (p_l / (1 - p_e + p_l))
  }
  return(weight1)
}
