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
    
    library(ggplot2)
    return(two_by_twos)
  } else if (length(control_vars > 0)) {

    # Controled ----------------------------------------------------------------
    # Predict Treatment
    control_formula <- paste(control_vars, collapse = " + ")
    control_formula <- as.formula(paste("treated ~", control_formula))

    data <- calculate_ds(data, control_formula)
    Sigma <- calculate_Sigma(data)
    one_minus_Sigma <- calculate_one_minus_Sigma(data)

    Sigma + one_minus_Sigma


    beta_hat_w <- calculate_beta_hat_w(data)
    N <- nrow(data)
    V_bd <- var(data$d_kt_til)*(N - 1)/N

    for (i in 1:nrow(two_by_twos)) {
      treated_group <- two_by_twos[i, "treated"]
      untreated_group <- two_by_twos[i, "untreated"]
      data1 <- data[data$treat_time %in% c(treated_group, untreated_group), ]
      if (treated_group < untreated_group) {
        data1 <- data1[data1$time < untreated_group, ]
      } else if (treated_group > untreated_group) {
        data1 <- data1[data1$time >= untreated_group, ]
      }

      skl <- calculate_weights_controled(data1, treated_group, untreated_group,
                                         V_bd)
    }
  }
}

#' Unpack variable names from formula
#' 
#' @param fomula
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





calculate_ds <- function(data, control_formula) {
  fit_treat <- lm(control_formula, data = data)
  data$d_it <- fit_treat$residuals

  # OLD WAY OF D_IT_TIL
  # dm_control_formula <- update(control_formula, . ~ . + 0 + factor(time) + factor(id))
  # fit_treat_dm <- lm(dm_control_formula, data = data)
  # data$d_it_til <- fit_treat_dm$residuals

  data$d_it_til <- data$d_it - data$d_i_bar - data$d_t_bar + data$d_bar_bar

  data$d_i_bar <- ave(data$d_it, data$id)
  data$d_t_bar <- ave(data$d_it, data$time)
  data$d_bar_bar <- mean(data$d_it)
  data$d_kt_bar <- ave(data$d_it, data$treat_time, data$time)
  data$d_k_bar <- ave(data$d_it, data$treat_time)
  data$d_ikt_til <- data$d_it - data$d_i_bar - (data$d_kt_bar - data$d_k_bar)
  data$d_kt_til <- (data$d_kt_bar - data$d_k_bar) - (data$d_t_bar - data$d_bar_bar)

  return(data)
}



calculate_Sigma <- function(data) {
  # TODO Test that within + between = 1
  N <- ncol(data)
  V_dw <- var(data$d_ikt_til)*(N - 1)/N
  V_d <- var(data$d_it_til)*(N - 1)/N
  Sigma <- V_wd/V_d
  return(Sigma)
}

calculate_one_minus_Sigma <- function(data) {
  N <- ncol(data)
  V_db <- var(data$d_kt_til)*(N - 1)/N
  V_d <- var(data$d_it_til)*(N - 1)/N
  one_minus_Sigma <- V_db/V_d
  return(one_minus_Sigma)
}

calculate_beta_hat_w <- function(data) {
  # TODO test equation 25
  beta_hat_w <- cov(data$outcome, data$d_ikt_til)
  V_d <- var(data$d_ikt_til)*(N - 1)/N
  return(beta_hat_w)
}


calculate_weights_controled <- function(data, treated_group,
                                        untreated_group, V_db) {
  # TODO test: between 0-1? (ask him), and sum to 1
  if (untreated_group == 99999) {
    # Treated vs untreated
    n_u <- sum(data$treat_time == untreated_group)
    n_k <- sum(data$treat_time == treated_group)
    V_bkl <- var(data$d_kt_til) # this is the V(d_kt) but only for these subgroups
    s_kl <- (n_u + n_k)^2 * V_bkl/V_db

  } else if (treated_group < untreated_group) {
    # early vs late (before late is treated)
    n_k <- sum(data$treat_time == treated_group)
    n_l <- sum(data$treat_time == untreated_group)
    V_bkl <- var(data$d_kt_til) # this is the V(d_kt) but only for these subgroups
    s_kl <- (n_k + n_l)^2 * V_bkl/V_db
  } else if (treated_group > untreated_group) {
    # late vs early (after early is treated)
    n_k <- sum(data$treat_time == untreated_group)
    n_l <- sum(data$treat_time == treated_group)
    V_bkl <- var(data$d_kt_til) # this is the V(d_kt) but only for these subgroups
    s_kl <- (n_k + n_l)^2 * V_bkl/V_db
  }
  return(s_kl)
}







# Function to predict treatment status with covariates
predict_treatment <- function(formula, data) {
  fit <- lm(control_formula, data = data)
  predict(fit)
}

calc_witin_var <- function(resid, groups) {

}
