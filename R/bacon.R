#' Goodman-Bacon Decomposition
#'
#' bacon() is a function that performs the Goodman-Bacon decomposition for
#'  differences-in-differences with variation in treatment timing (with or
#'  without time-varying covariates).
#'
#' @param formula an object of class "formula": a symbolic
#'  representation of the model to be fitted. Must be  of the form y ~ D + controls,
#'  where y is the outcome variable,  D is the binary
#'  treatment indicator, and `controls` can be any additional control variables. Do not
#'  include the fixed effects in the formula. If using `.` notation must be of 
#'  the form y ~ D + . - FE1 - FE2
#' @param data a data.frame containing the variables in the model.
#' @param id_var character, the name of id variable for units.
#' @param time_var character, the name of time variable.
#' @param quietly logical, default = FALSE, if set to TRUE then bacon() does not
#'  print the summary of estimates/weights by type (e.g. Treated vs Untreated)
#'
#' @return If control variables are included in the formula, then an object of
#'  class "list" with three elements:
#'  \item{Omega}{a number between 0 and 1, the weight of the within timing group
#'   coefficient}
#'  \item{beta_hat_w}{a number, the within timing group coefficient}
#'  \item{two_by_twos}{a data.frame with the covariate adjusted 2x2 estimates
#'   and weights}
#'
#' If not control variables are included then only the two_by_twos data.frame
#'  is returned. 
#'
#' @examples
#' # Castle Doctrine (Uncontrolled)
#' df_bacon <- bacon(l_homicide ~ post,
#'                   data = bacondecomp::castle,
#'                   id_var = "state",
#'                   time_var = "year")
#'
#' # Castle Doctrine (Controlled)
#' ret_bacon <- bacon(l_homicide ~ post + l_pop + l_income,
#'                    data = bacondecomp::castle,
#'                    id_var = "state",
#'                    time_var = "year")
#'
#' @import stats
#'
#' @export
bacon <- function(formula,
                  data,
                  id_var,
                  time_var, 
                  quietly = F) {
  # Evaluate formula in data environment
  formula <- formula(terms
                     (formula, data = data))
  
  # Unpack variable names and rename variables
  vars <- unpack_variable_names(formula)
  outcome_var <- vars$outcome_var
  treated_var <- vars$treated_var
  control_vars <- vars$control_vars
  data <- rename_vars(data, id_var, time_var, outcome_var, treated_var)
  
  # Check for a balanced panel
  bal <- aggregate(time ~ id,  data = data, FUN = length)
  balanced <- ifelse(all(bal$time == bal$time[1]), 1, 0)
  if (!balanced) {
    stop("Unbalanced Panel")
  }
  
  # Check for NA observations
  nas <- sum(is.na(data[, c("id", "time", "outcome", "treated")]))
  if (length(control_vars > 0)) {
    control_formula <- update(
      formula,
      paste0("treated ~ . - 1 - ", treated_var)
    )
    mm_control <- model.matrix(control_formula, data = data)
    nas_control <- 1 - (nrow(mm_control) == nrow(data))
    nas <- nas + nas_control
  }
  if (nas > 0) {
    stop("NA observations")
  }

  # Create 2x2 grid of treatment groups
  treatment_group_calc <- create_treatment_groups(data, control_vars,
                                                  return_merged_df = TRUE)
  two_by_twos <- treatment_group_calc$two_by_twos
  data <- treatment_group_calc$data

  # Uncontrolled ---------------------------------------------------------------
  if (length(control_vars) == 0) {
    # Iterate through treatment group dyads
    for (i in 1:nrow(two_by_twos)) {
      treated_group <- two_by_twos[i, "treated"]
      untreated_group <- two_by_twos[i, "untreated"]

      data1 <- subset_data(data, treated_group, untreated_group)

      # Calculate estimate and weight
      weight <- calculate_weights(data = data1,
                                  treated_group = treated_group,
                                  untreated_group = untreated_group)
      estimate <- lm(outcome ~ treated + factor(time) + factor(id),
                     data = data1)$coefficients[2]

      two_by_twos[i, "estimate"] <- estimate
      two_by_twos[i, "weight"] <- weight
    }
    
    # Rescale weights to sum to 1
    two_by_twos <- scale_weights(two_by_twos)
    
    if (quietly == F) {
      print_summary(two_by_twos)
    }
    
    return(two_by_twos)
    
  } else if (length(control_vars) > 0) {
    # Controled ----------------------------------------------------------------
    # Predict Treatment and calulate demeaned residuals
    control_formula <- update(
      formula,
      paste0("treated ~ . + factor(time) + factor(id) -", treated_var)
    )
  
    data <- run_fwl(data, control_formula)

    # Calculate within treatment group estimate and its weight
    Omega <- calculate_Omega(data)
    beta_hat_w <- calculate_beta_hat_w(data)

    # Collapse controls and predicted treatment to treatment group/year level
    r_collapse_x_p <- collapse_x_p(data, control_formula)
    data <- r_collapse_x_p$data
    g_control_formula <- r_collapse_x_p$g_control_formula

    # Iterate through treatment group dyads
    for (i in 1:nrow(two_by_twos)) {
      treated_group <- two_by_twos[i, "treated"]
      untreated_group <- two_by_twos[i, "untreated"]
      data1 <- data[data$treat_time %in% c(treated_group, untreated_group), ]

      # Calculate between group estimate and weight
      weight_est <- calc_controlled_beta_weights(data1, g_control_formula)
      s_kl <- weight_est$s_kl
      beta_hat_d_bkl <- weight_est$beta_hat_d_bkl

      two_by_twos[i, "weight"] <- s_kl
      two_by_twos[i, "estimate"] <- beta_hat_d_bkl
    }

    # Rescale weights to sum to 1
    two_by_twos <- scale_weights(two_by_twos)
    
    if (quietly == F) {
      print_summary(two_by_twos)
    }
    
    r_list <- list("beta_hat_w" = beta_hat_w,
                   "Omega" = Omega,
                   "two_by_twos" = two_by_twos)
    return(r_list)
  }
}


#' Unpack Variable Names from Formula
#'
#' @param formula formula call from bacon()
#'
#' @return A list with 3 elements:
#' \item{outcome_var}{name of the outcome variable from the LHS of formula}
#' \item{treated_var}{name of the treatment indicator variable, the first
#'  variable from the RHS of formula}
#' \item{control_vars}{All variables on RHS of formula except the treated_var.
#'  Only need to know if length(control_vars) > 0}
#'
#' @noRd
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
#'
#' Rename id, time, outcome, and treatment indicator variables
#'
#' @param data a data.frame
#' @param id_var character, name of id variable
#' @param time_var character, name of time variable
#' @param outcome_var character, name of outcome variable
#' @param treated_var character, name of binary treatment variable
#'
#' @return A data.frame with renamed columns
#'
#' @noRd
rename_vars <- function(data, id_var, time_var, outcome_var, treated_var) {
  colnames(data)[colnames(data) == id_var] <- "id"
  colnames(data)[colnames(data) == time_var] <- "time"
  colnames(data)[colnames(data) == outcome_var] <- "outcome"
  colnames(data)[colnames(data) == treated_var] <- "treated"
  return(data)
}


#' Create Grid of Treatment Groups
#'
#' @param data a data.frame used to create groups - MUST obey naming convention
#' used in `bacon()`. i.e. columns are ["id", "time", "outcome", "treated"]
#'
#' @param return_merged_df Defaults to `FALSE` whether to return merged data
#' as well as grid of treatment groups.
#' @param control_vars list of control variables
#'
#' @return If return_merge_df is `TRUE` then a list with two elements:
#' \item{two_by_twos}{data.frame with all treatment group dyads and NA
#'  weight/estimate columns}
#' \item{data}{updated data.frame with the new variable treat_time}
#'
#' If return_merge_df is `FALSE` then only the two_by_twos data.frame is
#'  returned
#'
#' @noRd
create_treatment_groups <- function(data, control_vars,
                                    return_merged_df = FALSE) {
  df_treat <- data[data$treated == 1, ]
  df_treat <- df_treat[, c("id", "time")]
  df_treat <- aggregate(time ~ id,  data = df_treat, FUN = min)
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

  if (length(control_vars) == 0) {
    # Create data.frame of all posible 2x2 estimates. Dyads may appear twice as
    # treatment groups can play the roll of both earlier and later treated
    two_by_twos <- expand.grid(unique(data$treat_time),
                               unique(data$treat_time))
    colnames(two_by_twos) <- c("treated", "untreated")
    two_by_twos <- two_by_twos[!(two_by_twos$treated ==
                                 two_by_twos$untreated), ]
    two_by_twos <- two_by_twos[!(two_by_twos$treated == 99999), ]
    # Remove first period
    two_by_twos <- two_by_twos[!(two_by_twos$treated == min(data$time)), ]
    two_by_twos[, c("estimate", "weight")] <- NA
    # Classify estimate "type"
    two_by_twos[, "type"] <- ifelse(two_by_twos$untreated == 99999,
                                    "Treated vs Untreated",
                                    ifelse(two_by_twos$untreated == 
                                             min(data$time), 
                                           "Later vs Always Treated", 
                                           ifelse(two_by_twos$treated >
                                                    two_by_twos$untreated,
                                                  "Later vs Earlier Treated",
                                                  "Earlier vs Later Treated")))
    
  } else if (length(control_vars) > 0) {
    # In the controlled decomposition, each dyad only appears once becasue we
    # do not make the distinction between earlier vs later treated
    two_by_twos <- data.frame()
    for (i in unique(data$treat_time[data$treat_time != 99999])) {
      for (j in unique(data$treat_time)) {
        if (j > i) {
          two_by_twos1 <- data.frame(treated = i, untreated = j, weight = NA,
                                     estimate = NA)
          two_by_twos <- rbind(two_by_twos, two_by_twos1)
        }
      }
    }
    two_by_twos[, "type"] <- ifelse(two_by_twos$untreated == 99999,
                                    "Treated vs Untreated",
                                    ifelse(two_by_twos$treated == 
                                             min(data$time) | 
                                             two_by_twos$untreated == 
                                             min(data$time), 
                                           "Later vs Always Treated", 
                                           "Both Treated"))
  }

  # Whether or not to return the merged data too.
  if (return_merged_df == TRUE) {
    return_data <- list("two_by_twos" = two_by_twos,
                        "data" = data)
  } else {
    return_data <- two_by_twos
  }
  return(return_data)
}


#' Subset Data
#'
#' In the uncontrolled decomposition, subset_data() subsets the data to only
#'  two treatment groups, and subsets the time to either before the later group
#'  has been treated or after the earlier group has been treated.
#'
#' @param data a data.frame
#' @param treated_group integer, initial treatment time of group acting as
#'  treated
#' @param untreated_group interger, initial treatment time of group acting as
#'  untreated
#'
#' @return subsetted data.frame
#'
#' @noRd
subset_data <- function(data, treated_group, untreated_group) {
  data <- data[data$treat_time %in% c(treated_group, untreated_group), ]
  if (treated_group < untreated_group) {
    data <- data[data$time < untreated_group, ]
  } else if (treated_group > untreated_group) {
    data <- data[data$time >= untreated_group, ]
  }
  return(data)
}


#' Calculate Weights for 2x2 Grid (uncontrolled)
#'
#'  Calculated weights for the uncontrolled decomposition using:
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
#'
#' @noRd
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


#' Scale Weights
#'
#' scale_weights scales the two_by_two weights so that they sum to 1
#'
#' @param two_by_twos data.frame containing 2x2 weights
#'
#' @return two_by_twos data.frame with rescaled weight column
#'
#' @noRd
scale_weights <- function(two_by_twos) {
  sum_weight <- sum(two_by_twos$weight)
  two_by_twos$weight <- two_by_twos$weight / sum_weight
  return(two_by_twos)
}

#' Run Frisch-Waugh-Lowell Regression
#'
#' Predicts treatment using time-varing covariates and calculates predicted
#'  values/residuals
#'
#' @param data a data.frame
#' @param control_formula a fomula
#'
#' @return updated data.frame with predictions/residuals (and their
#'  transformations) from the FWL regession
#'
#' @noRd
run_fwl <- function(data, control_formula) {
  fit_fwl <- lm(control_formula, data = data)
  data$p <- predict(fit_fwl)
  data$d <- fit_fwl$residuals
  data$d_it <- data$d

  # demean residulas
  data$d_i_bar <- ave(data$d, data$id)
  data$d_t_bar <- ave(data$d_it, data$time)
  data$d_bar_bar <- mean(data$d_it)
  data$d_it_til <- data$d_it - data$d_i_bar - data$d_t_bar + data$d_bar_bar

  data$d_kt_bar <- ave(data$d_it, data$treat_time, data$time)
  data$d_k_bar <- ave(data$d_it, data$treat_time)
  data$d_ikt_til <- data$d_it - data$d_i_bar - (data$d_kt_bar - data$d_k_bar)
  data$d_kt_til <- (data$d_kt_bar - data$d_k_bar) -
                   (data$d_t_bar - data$d_bar_bar)
  return(data)
}


#' Calculate Omega
#'
#' Calculate the weight for the within estimate in the controlled decompositon
#'
#' @param data a data frame with the columns d_ikt_til and d_it_til
#'
#' return Omega
#'
#' @noRd
calculate_Omega <- function(data) {
  # TODO Test that within + between = 1
  N <- nrow(data)
  Vd_w <- var(data$d_ikt_til) * (N - 1) / N
  V_d <- var(data$d_it_til) * (N - 1) / N
  Omega <- Vd_w / V_d
  return(Omega)
}

#' Calculate 1 - Omega
#'
#' Calculate the weight for the between estimates. This is used for testing
#'  purposes only.
#'
#' @param data a data frame with the columns d_kt_til and d_it_til
#'
#' @return 1 - Omega
#'
#' @noRd
calculate_one_minus_Omega <- function(data) {
  N <- nrow(data)
  Vd_b <- var(data$d_kt_til) * (N - 1) / N
  V_d <- var(data$d_it_til) * (N - 1) / N
  one_minus_Omega <- Vd_b / V_d
  return(one_minus_Omega)
}


#' Calculate Within Estimate
#'
#' @param data a data.frame with the columns outcome and d_ikt_til
#'
#' @return beta_hat_w, the within estimate
#'
#' @noRd
calculate_beta_hat_w <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_ikt_til) * (N - 1) / N
  Vd_w <- var(data$d_ikt_til) * (N - 1) / N
  beta_hat_w <- C / Vd_w
  return(beta_hat_w)
}


#' Calculate Between Estimate
#'
#' Calculate the overall between estimate (not decomposed by treatement group
#'  dyads). This is used for testing purposes only.
#'
#' @param data a data.frame with the columns outcome and d_kt_til
#'
#' @return beta_hat_w, the between estimate
#'
#' @noRd
calculate_beta_hat_b <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_kt_til) * (N - 1) / N
  Vd_b <- var(data$d_kt_til) * (N - 1) / N
  beta_hat_b <- C / Vd_b
  return(beta_hat_b)
}


#' Collapse Control Variable and Predicted Treatment
#'
#' Collapse controll variables and predicted treatment to treatment group/time
#'  level.
#'
#' @param data a data.frame
#' @param control_formula, a formula
#'
#' @return a list of two elements:
#' \item{data}{updated data frame with new treatment group/time level variables}
#' \item{g_control_formula}{the control formula but with group level variable
#'  names}
#'
#' @noRd
collapse_x_p <- function(data, control_formula) {
  # Group level Xs
  f1 <- update(control_formula, . ~ . - factor(time) - factor(id) - 1)
  control_data <- data.frame(model.matrix(f1, data = data))

  my_ave <- function(x, data) {
    ave(x, data$treat_time, data$time)
  }

  g_control_data <- data.frame(sapply(control_data, my_ave, data))
  colnames(g_control_data) <- paste0("g_", colnames(g_control_data))

  # Group level p
  data$g_p <- ave(data$p, data$treat_time, data$time)
  data <- cbind(data, g_control_data)

  g_control_formula <- as.formula(paste("~", paste(colnames(g_control_data),
                                                   collapse = " + ")))

  r_list <- list(data = data, g_control_formula = g_control_formula)

  return(r_list)
}


#' Calculate Treatment Variance
#'
#' @param data a data.frame
#'
#' @return a list of two items:
#' \item{data}{a data.frame with the new demeaned treatment variable}
#' \item{VD}{the varince of the demeaned treatment variable}
#'
#' @noRd
calc_VD <- function(data) {
  fit_D <- lm(treated ~ factor(id) + factor(time), data = data)
  data$Dtilde <- fit_D$residuals
  N <- nrow(data)
  VD <- var(data$Dtilde) * (N - 1) / N
  r_list <- list(data = data, VD = VD)
  return(r_list)
}


#' Partial Out FE
#'
#' Partial out time and unit fixed effects from group level control and
#'  and predicted treatment variables
#'
#' @param data, a data.frame with the variables in the g_control_formula and
#'  p
#' @param g_control_formula control formula with group level variable names
#'
#' @return updated data.frame with new partialled out group level variables
#'
#' @noRd
partial_group_x <- function(data, g_control_formula) {
  g_vars <- unlist(strsplit(as.character(g_control_formula)[2], " \\+ "))
  for (v in g_vars) {
    g_formula <- as.formula(paste0(v, " ~ factor(id) + factor(time)"))
    fit_g <- lm(g_formula, data = data)
    data[, paste0("p_", v)] <- fit_g$residuals
  }
  return(data)
}


#' Calculate pgjtile and Rsq
#'
#' @param data, a data.frame with the variables in g_control_formula and
#'  Dtilde
#' @param g_control_formula control formula with group level variable names
#'
#' @return a list with two elements:
#' \item{data}{updated data.frame with new variable pgjtilde}
#' \item{Rsq}{Rsq from the regression to generate pgjtilde}
#'
#' @noRd
calc_pgjtile <- function(data, g_control_formula) {
  g_vars <- unlist(strsplit(as.character(g_control_formula)[2], " \\+ "))
  p_g_vars <- paste0("p_", g_vars)

  p_g_control_formula <- formula(paste("Dtilde ~", paste(p_g_vars,
                                                         collapse = " + ")))

  fit_pgj <- lm(p_g_control_formula, data = data)
  data$pgjtilde <- predict(fit_pgj)
  Rsq <- summary(fit_pgj)$r.squared
  r_list <- list(data = data, Rsq = Rsq)
  return(r_list)
}


#' Calculate Vdp
#'
#' @param data a data.frame with the variables g_p, id, time, and ptilde
#'
#' @return list with two elements
#' \item{data}{data.frame with new variables ptilde and dp}
#' \item{Vdp}{Variance of dp}
#'
#' @noRd
calc_Vdp <- function(data) {
  N <- nrow(data)
  fit_p <- lm(g_p ~ factor(id) + factor(time), data = data)
  data$ptilde <- fit_p$residuals
  data$dp <- data$pgjtilde - data$ptilde
  Vdp <- var(data$dp) * (N - 1) / N
  r_list <- list(data = data, Vdp = Vdp)
  return(r_list)
}


#' Calculate BD
#'
#' @param data a data.frame with all variables in the g_control_formula,
#'  outcome, treated, id, and time
#' @param g_control_formula control formula with group level variable names
#' @return BD
#'
#' @noRd
calc_BD <- function(data, g_control_formula) {
  BD_formula <- update(g_control_formula,
                       outcome ~ treated + . + factor(id) + factor(time))
  fit_BD <- lm(BD_formula, data = data)
  BD <- fit_BD$coefficients["treated"]
  return(BD)
}


#' Calculate Bb
#'
#' @param data a data.frame with the variables outcome, dp, time, and id
#'
#' @return Bb
#'
#' @noRd
calc_Bb <- function(data) {
  fit_Bb <- lm(outcome ~ dp + factor(time) + factor(id),
               data = data)
  Bb <- fit_Bb$coefficients["dp"]
  return(Bb)
}


#' Calculate Dyad Covariate Adjusted Between Beta Hat
#'
#' @param Rsq R squared from the regression in calc_pgjtilde
#' @param VD variance of the demeaned treatment indicator (in the dyad)
#' @param BD coeffieicient on treatment from a regression of the outcome on
#'  treatment and group level control variables with unite and time FE
#' @param Vdp variance of dp (pgjtilde - ptilde)
#' @param Bb coefficient on dp from a regression of the outcome on dp with unit
#'  and time fixed effects
#'
#' @return dyad's covariate adjusted between estimate
#'
#' @noRd
calculate_beta_hat_d_bkl <- function(Rsq, VD, BD, Vdp, Bb) {
  beta_hat_d_bkl <- ( (1 - Rsq) * VD * BD + Vdp * Bb) / ( (1 - Rsq) * VD + Vdp)
  return(beta_hat_d_bkl)
}

#' Calculate Dyad Weight (Covariate Adjusted)
#'
#' @param N observations in dyad
#' @param Rsq R squared from the regression in calc_pgjtilde
#' @param VD variance of the demeaned treatment indicator (in the dyad)
#' @param Vdp variance of dp (pgjtilde - ptilde)
#'
#' @return the weight given to the dyad's estimate
#'
#' @noRd
calculate_s_kl <- function(N, Rsq, VD, Vdp) {
  s_kl <- N ^ 2 * ( (1 - Rsq) * VD + Vdp)
  return(s_kl)
}


#' Calculate Covariate Adjusted Betas and Weights
#'
#' Wrapper function for finding the dyad beta/weight in the controled
#'  decomposition
#'
#' @param data a data.frame
#' @param g_control_formula formula with group/timelevel control variables
#'
#' @return a list with two elements
#' \item{s_kl}{the weight given to the dyad kl's estimate}
#' \item{beta_hat_d_bkl}{the covariate adjusted 2x2 estimate for the dyad kl}
#'
#' @noRd
calc_controlled_beta_weights <- function(data, g_control_formula) {
  r_calc_VD <- calc_VD(data)
  VD <- r_calc_VD$VD
  data <- r_calc_VD$data

  data <- partial_group_x(data, g_control_formula)

  r_calc_pgjtilde <- calc_pgjtile(data, g_control_formula)
  Rsq <- r_calc_pgjtilde$Rsq
  data <- r_calc_pgjtilde$data

  r_calc_Vdp <- calc_Vdp(data)
  Vdp <- r_calc_Vdp$Vdp
  data <- r_calc_Vdp$data

  BD <- calc_BD(data, g_control_formula)
  Bb <- calc_Bb(data)

  N <- nrow(data)

  s_kl <- calculate_s_kl(N, Rsq, VD, Vdp)
  beta_hat_d_bkl <- calculate_beta_hat_d_bkl(Rsq, VD, BD, Vdp, Bb)

  r_list <- list(s_kl = s_kl, beta_hat_d_bkl = beta_hat_d_bkl)
  return(r_list)
}

print_summary <- function(two_by_twos, return_df = FALSE) {
  sum_df <- aggregate(weight ~ type,
                      data = two_by_twos,
                      FUN = sum)
  sum_df$avg_est <- c(by(two_by_twos, two_by_twos$type,
                         function(x) round(weighted.mean(x$estimate, x$weight), 5)))
  sum_df$weight <- round(sum_df$weight, 5)
  print(sum_df)
  
  if (return_df == TRUE) {
    return(sum_df)
  }
}
