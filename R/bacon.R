#' Goodman-Bacon Decomposition
#'
#' bacon() is a function that perfroms the Goodman-Bacon decomposition for
#'  differences-in-differences with variation in treatment timing.
#'
#' @param formula a symbolic representation of the
#'  model to be fitted.
#' @param data a data.frame containing the variables in the model.
#' @param id_var character, the name of id variable for units.
#' @param time_var character, the name of time variable.
#'
#' @author Evan Flack
#' @return data.frame of all 2x2 estimates and weights
#'
#' @examples
#' # Castle Doctrine -----------------------------------------------------------
#' df_bacon <- bacon(l_homicide ~ post,
#'                   data = bacon::castle,
#'                   id_var = "state",
#'                   time_var = "year")
#'
#' \donttest{
#'
#' library(ggplot2)
#' ggplot(df_bacon) +
#'   aes(x = weight, y = estimate, shape = factor(type)) +
#'   labs(x = "Weight", y = "Estimate", shape = "Type") +
#'   geom_point()
#'
#'   }
#'
#'
#' @export

data <- bacon::castle
formula <- l_homicide ~ post + l_pop + l_income
id_var <- "state"
time_var <- "year"

bacon <- function(formula,
                  data,
                  id_var,
                  time_var) {
  vars <- unpack_variable_names(formula)
  outcome_var <- vars$outcome_var
  treated_var <- vars$treated_var
  control_vars <- vars$control_vars

  data <- rename_vars(data, id_var, time_var, outcome_var, treated_var)


  # Check for NA observations
  nas <- sum(is.na(data[, c("id", "time", "outcome", "treated", control_vars)]))
  if (nas > 0) {
    stop("NA observations")
  }

  # Check for balanced panel
  bal <- stats::aggregate(time ~ id,  data = data, FUN = length)
  balanced <- ifelse(all(bal$time == bal$time[1]), 1, 0)
  if (!balanced) {
    stop("Unbalanced Panel")
  }

  # Create grid of treatment groups
  treatment_group_calc <- create_treatment_groups(data,
                                                  return_merged_df = TRUE)
  two_by_twos <- treatment_group_calc$two_by_twos
  data <- treatment_group_calc$data

  # First period in the panel
  first_period <- min(data$time)

  # Uncontrolled ---------------------------------------------------------------
  if (length(control_vars == 0)) {
    for (i in 1:nrow(two_by_twos)) {
      treated_group <- two_by_twos[i, "treated"]
      untreated_group <- two_by_twos[i, "untreated"]

      data1 <- subset_data(data, treated_group, untreated_group)

      weight1 <- calculate_weights(data = data1,
                                   treated_group = treated_group,
                                   untreated_group = untreated_group)

      estimate1 <- stats::lm(outcome ~ treated + factor(time) + factor(id),
                             data = data1)$coefficients[2]

      two_by_twos[i, "estimate"] <- estimate1
      two_by_twos[i, "weight"] <- weight1
    }
  } else if (length(control_vars > 0)) {
    # Controled ----------------------------------------------------------------
    # Predict Treatment and calulate demeaned residuals
    control_formula <- update(formula,
                              paste0("treated ~ . + factor(time) + factor(id) -",
                                     treated_var))
    data <- calculate_ds(data, control_formula)
    # Ca
    Sigma <- calculate_Sigma(data)
    # one_minus_Sigma <- calculate_one_minus_Sigma(data)
    beta_hat_w <- calculate_beta_hat_w(data)

    N <- nrow(data)
    V_bd <- var(data$d_kt_til)*(N - 1)/N

    for (i in 1:nrow(two_by_twos)) {
      treated_group <- two_by_twos[i, "treated"]
      untreated_group <- two_by_twos[i, "untreated"]
      data1 <- data[data$treat_time %in% c(treated_group, untreated_group), ]
      skl <- calculate_weights_controled(data1, treated_group, untreated_group,
                                         V_bd)
      fit_2x2 <- lm(outcome ~ treated + factor(id) + factor(time), data = data1)
      beta_2x2_kl <- fit_2x2$coefficients["treated"]

      
      two_by_twos[i, "weight"] <- skl
    }
  }
}

#' Unpack Variable Names from Formula
unpack_variable_names <- function(formula) {
  outcome_var <- as.character(formula)[2]
  right_side_vars <- as.character(formula)[3]
  right_side_vars <- strsplit(right_side_vars, " \\+ ")[[1]]
  treated_var <- right_side_vars[1]
  control_vars <- right_side_vars[-1]
  r_list <- list(outcome_var = outcome_var, treated_var = treated_var,
                 control_vars = control_vars)
  return(r_list)
}

#' Rename Variables
rename_vars <- function(data, id_var, time_var, outcome_var, treated_var) {
  colnames(data)[which(colnames(data) %in%
                         c(id_var,
                           time_var,
                           outcome_var,
                           treated_var))] <- c("id", "time", "outcome",
                                               "treated")
  return(data)
}

#' Create Grid of Treatment Groups
#'
#' @param data dataset used to create groups - MUST obey naming convention used
#' in `bacon()`. i.e. columns are ["id", "time", "outcome", "treated"]
#' @param return_merged_df Defaults to `FALSE` whether to return merged data
#' as well as grid of treatment groups.
#'
#' @return data.frame describing treatment groups and empty weight and estimate
#' column set to 0.
create_treatment_groups <- function(data, return_merged_df = FALSE){
  df_treat <- data[data$treated == 1, ]
  df_treat <- df_treat[, c("id", "time")]
  df_treat <- stats::aggregate(time ~ id,  data = df_treat, FUN = min)
  colnames(df_treat) <- c("id", "treat_time")
  data <- merge(data, df_treat, by = "id", all.x = T)
  data[is.na(data$treat_time), "treat_time"] <- 99999

  # Check for weakly increasing treatment
  inc <- ifelse(data$treat_time == 99999, 1,
                ifelse(data$time >= data$treat_time & data$treated == 1, 1,
                       ifelse(data$time < data$treat_time & data$treated == 0,
                              1, 0)))
  if (!all(as.logical(inc))) {
    stop("Treatment not weakly increasing with time")
  }
  # create data.frame of all posible 2x2 estimates
  two_by_twos <- expand.grid(unique(data$treat_time),
                             unique(data$treat_time))
  colnames(two_by_twos) <- c("treated", "untreated")
  two_by_twos <- two_by_twos[!(two_by_twos$treated == two_by_twos$untreated), ]
  two_by_twos <- two_by_twos[!(two_by_twos$treated == 99999), ]
  # Remove first period
  two_by_twos <- two_by_twos[!(two_by_twos$treated == min(data$time)), ]
  two_by_twos[, c("estimate", "weight")] <- 0


  # Whether or not to return the merged data too.
  if (return_merged_df == TRUE) {
    return_data <- list("two_by_twos" = two_by_twos,
                        "data" = data)
  } else {
    return_data <- two_by_twos
  }
  return(return_data)
}

subset_data <- function(data, treated_group, untreated_group) {
  data <- data[data$treat_time %in% c(treated_group, untreated_group), ]
  if (treated_group < untreated_group) {
    data <- data[data$time < untreated_group, ]
  } else if (treated_group > untreated_group) {
    data <- data[data$time >= untreated_group, ]
  }
  return(data)
}

#' Calculate Weights for 2x2 Grid
#'
#'  Calculated weights using:
#'  n_u - observations in untreated group,
#'  n_k - observations in earlier treated group,
#'  n_l - observations in later treated group,
#'  D_k - proportion of time the earlier treated group was treated,
#'  D_l - proportion of time the later treated group was treated.
#'
#' @param data a data.frame
#' @param treated_group the identifier of the treated group
#' @param untreated_group the identifier of the untreated group
#'
#' @return Scalar weight for 2x2 grid
calculate_weights <- function(data,
                              treated_group,
                              untreated_group){
  if (untreated_group == 99999) {
    # Treated vs untreated
    n_u <- sum(data$treat_time == untreated_group)
    n_k <- sum(data$treat_time == treated_group)
    n_ku <- n_k / (n_k + n_u)
    D_k <- mean(data[data$treat_time == treated_group, "treated"])
    V_ku <- n_ku * (1 - n_ku) * D_k * (1 - D_k)
    weight1 <- (n_k + n_u) ^ 2 * V_ku
  } else if (treated_group < untreated_group) {
    # early vs late (before late is treated)
    n_k <- sum(data$treat_time == treated_group)
    n_l <- sum(data$treat_time == untreated_group)
    n_kl <- n_k / (n_k + n_l)
    D_k <- mean(data[data$treat_time == treated_group, "treated"])
    D_l <- mean(data[data$treat_time == untreated_group, "treated"])
    V_kl <- n_kl * (1 - n_kl) * (D_k - D_l) / (1 - D_l) * (1 - D_k) / (1 - D_l)
    weight1 <- ( (n_k + n_l) * (1 - D_l) ) ^ 2 * V_kl
  } else if (treated_group > untreated_group) {
    # late vs early (after early is treated)
    n_k <- sum(data$treat_time == untreated_group)
    n_l <- sum(data$treat_time == treated_group)
    n_kl <- n_k / (n_k + n_l)
    D_k <- mean(data[data$treat_time == untreated_group, "treated"])
    D_l <- mean(data[data$treat_time == treated_group, "treated"])
    V_kl <- n_kl * (1 - n_kl) * (D_l / D_k) * (D_k - D_l) / (D_k)
    weight1 <- ( (n_k + n_l) * D_k) ^ 2 * V_kl
  }
  return(weight1)
}

#' Calculate d_it and Variants
calculate_ds <- function(data, control_formula) {
  
  #  Quoting:
  # To see how the controlled DD coefficient is identified first remove unit- and time-
  # means(indicated by tildes) and then estimate a Frisch-Waugh regression that partials
  # 𝑿𝑿�𝒊𝒊𝒂𝒂out of 𝐷𝐷�𝑖𝑖𝑖𝑖
  
  # Think we need to demean here, before hand
  # predict treatment
  fit_treat <- lm(control_formula, data = data)
  # residulals
  data$d_it <- fit_treat$residuals
  # demean residulas
  data$d_i_bar <- ave(data$d_it, data$id)
  data$d_t_bar <- ave(data$d_it, data$time)
  data$d_bar_bar <- mean(data$d_it)
  data$d_it_til <- data$d_it - data$d_i_bar - data$d_t_bar + data$d_bar_bar
  
  data$d_kt_bar <- ave(data$d_it, data$treat_time, data$time)
  data$d_k_bar <- ave(data$d_it, data$treat_time)
  data$d_ikt_til <- data$d_it - data$d_i_bar - (data$d_kt_bar - data$d_k_bar)
  data$d_kt_til <- (data$d_kt_bar - data$d_k_bar) - (data$d_t_bar - data$d_bar_bar)
  
  data$D_it_til <- data$treated - ave(data$treated, data$id) - ave(data$treated, data$time) + mean(data$treated)
  
  return(data)
}

#' Calculate Sigma
calculate_Sigma <- function(data) {
  # TODO Test that within + between = 1
  N <- nrow(data)
  V_wd <- var(data$d_ikt_til)*(N - 1)/N
  V_d <- var(data$d_it_til)*(N - 1)/N
  Sigma <- V_wd/V_d
  return(Sigma)
}

#' Calculate 1 - Sigma
calculate_one_minus_Sigma <- function(data) {
  N <- nrow(data)
  V_db <- var(data$d_kt_til)*(N - 1)/N
  V_d <- var(data$d_it_til)*(N - 1)/N
  one_minus_Sigma <- V_db/V_d
  return(one_minus_Sigma)
}

#' Calculate beta_hat_w
calculate_beta_hat_w <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_ikt_til)*(N - 1)/N
  V_d <- var(data$d_ikt_til)*(N - 1)/N
  beta_hat_w <- C/V_d
  return(beta_hat_w)
}

calculate_beta_hat_b <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_kt_til)*(N - 1)/N
  V_d <- var(data$d_kt_til)*(N - 1)/N
  beta_hat_b <- C/V_d
  return(beta_hat_b)
}

calculate_ps <- function(data, control_formula) {
  fit_treat <- lm(control_formula, data = data)
  # predicted
  data$p_it <- predict(fit_treat)
  data$p_jt_bar <- ave(data$p_it, data$treat_time, data$time)
  data$p_jt_til <- data$p_it_tile - data$p_jt_bar
}

calculate_weights_controled <- function(data, treated_group,
                                        untreated_group, V_db) {
  # TODO test: between 0-1? (ask him), and sum to 1
  # Treated vs untreated
  n_k <- sum(data$treat_time == untreated_group)
  n_l <- sum(data$treat_time == treated_group)
  N <- nrow(data)
  V_bkl <- var(data$d_kt_til)*(N - 1)/N # this is the V(d_kt) but only for these subgroups
  s_kl <- (n_k + n_l)^2 * V_bkl/V_db
  return(s_kl)
}

calculate_VD_kl <- function(data) {
  D_jt_til <- mean(data$treated) - ave(data$treated, data$treat_time)
  fit <- lm(D_jt_til ~ factor(id) + factor(time), data = data)
  resid <- fit$residuals
  N <- nrow(data)
  VD_kl <- var(resid)*(N - 1)/N
  return(VD_kl)
}

calculate_Vp_bkl <- function(data) {
  fit <- lm(data$p_jt_til ~ factor(id) + factor(time), 
            data = data)
  resid <- fit$residulals
  N <-nrow(data)
  Vp_bkl <- var(resid)*(N - 1)/N
}

calculate_beta_hat_22_kl <- function(data) {
  fit <- lm(outcome ~ treated + factor(id) + factor(time), data = data)
  beta_hat_22_kl <- fit$coefficients["treated"]
  return(beta_hat_22_kl)
}

calculate_beta_hat_p_bkl <- function(data) {
  data$y_jt_bar <- ave(data$outcome, data$treat_time, data$time)
  fit <- lm(y_jt_bar ~ p_jt_bar + factor(id) + factor(time), data = data)
  beta_hat_p_bkl <- fit$coefficients["p_jt_bar"]
  return(beta_hat_p_bkl)
}
