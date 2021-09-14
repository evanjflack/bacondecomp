#' Goodman-Bacon Decomposition
#'
#' \code{bacon()} is a function that performs the Goodman-Bacon decomposition for
#'  differences-in-differences with variation in treatment timing (with or
#'  without time-varying covariates).
#'
#' @param formula an object of class "formula". Must be  of the form y ~ D + controls,
#'  where y is the outcome variable,  D is the binary treatment indicator, and
#'  `controls` can be any additional control variables. Do not include the
#'  fixed effects in the formula.
#' @param data a data.frame containing the variables in the model.
#' @param id_var character, the name of id variable for units.
#' @param time_var character, the name of time variable.
#' @param quietly logical. If TRUE then bacon() does not
#'  print the summary of estimates/weights by type (e.g. Treated vs Untreated).
#'  Default is TRUE. You can also use \code{bacon_summary()} on the result to view this.
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
#'   data = bacondecomp::castle, id_var = "state", time_var = "year"
#' )
#'
#' # Castle Doctrine (Controlled)
#' ret_bacon <- bacon(l_homicide ~ post + l_pop + l_income,
#'   data = bacondecomp::castle, id_var = "state", time_var = "year"
#' )
#'
#' @import stats fixest
#'
#' @export
bacon <- function(formula, data, id_var, time_var, quietly = TRUE) {

  # Evaluate formula in data environment
  formula <- formula(terms(formula, data = data))
  
  # Unpack variable names and rename variables
  vars <- unpack_variable_names(formula)
  outcome_var <- vars$outcome_var
  treated_var <- vars$treated_var
  control_vars <- vars$control_vars
  data <- rename_vars(data, id_var, time_var, outcome_var, treated_var)

  # Check for NA observations
  nas <- sum(is.na(data[, c("id", "time", "outcome", "treated")]))
  if (length(control_vars > 0)) {
    control_formula <- as.formula(
      paste("~ ", paste(control_vars, collapse = " + "))
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
    return_merged_df = TRUE
  )
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
      weight <- calculate_weights(
        data = data1,
        treated_group = treated_group,
        untreated_group = untreated_group
      )
      estimate <- fixest::feols(outcome ~ treated | time + id,
        data = data1
      )$coefficients[1]

      two_by_twos[i, "estimate"] <- estimate
      two_by_twos[i, "weight"] <- weight
    }

    # Rescale weights to sum to 1
    two_by_twos <- scale_weights(two_by_twos)

    if (quietly == F) {
      bacon_summary(two_by_twos)
    }

    return(two_by_twos)
    # Controled ------------------------------------------------------------------
  } else if (length(control_vars) > 0) {
    # Predict Treatment and calulate demeaned residuals
    control_formula <- as.formula(
      paste("treated ~ ", paste(control_vars, collapse = " + "), " | time + id")
    )

    data <- run_fwl(data, control_formula)

    # Calculate within treatment group estimate and its weight
    Omega <- calculate_Omega(data)
    beta_hat_w <- calculate_beta_hat_w(data)

    # Collapse controls and predicted treatment to treatment group/year level
    r_collapse_x_p <- collapse_x_p(data, control_vars)
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
      bacon_summary(two_by_twos)
    }

    r_list <- list(
      "beta_hat_w" = beta_hat_w,
      "Omega" = Omega,
      "two_by_twos" = two_by_twos
    )
    return(r_list)
  }
}


#' Summary of Goodman-Bacon Decomposition
#' 
#' Uses the two-by-two output produced by 
#' \code{bacondecomp::bacon} to produce 
#' average 2x2 estimate and total weight for the following three comparisons: 
#' Earlier vs. Later (Good), Treated vs. Untreated (Good), and 
#' Later vs. Earlier (Bad).
#' 
#' @param two_by_twos Data.frame produced by \code{bacondecomp::bacon}.
#' @param return_df Logical. If TRUE, then the summary data.frame is returned.
#'   Default is False.
#'   
#' @export
bacon_summary <- function(two_by_twos, return_df = FALSE) {
  sum_df <- aggregate(weight ~ type,
                      data = two_by_twos,
                      FUN = sum
  )
  sum_df$avg_est <- c(by(
    two_by_twos, two_by_twos$type,
    function(x) round(weighted.mean(x$estimate, x$weight), 5)
  ))
  sum_df$weight <- round(sum_df$weight, 5)
  print(sum_df)
  
  if (return_df == TRUE) {
    return(sum_df)
  }
}


unpack_variable_names <- function(formula) {
  outcome_var <- as.character(formula)[2]
  right_side_vars <- as.character(formula)[3]
  right_side_vars <- strsplit(right_side_vars, " \\+ ")[[1]]
  treated_var <- right_side_vars[1]
  control_vars <- right_side_vars[-1]
  r_list <- list(
    outcome_var = outcome_var, treated_var = treated_var,
    control_vars = control_vars
  )
  return(r_list)
}

rename_vars <- function(data, id_var, time_var, outcome_var, treated_var) {
  colnames(data)[colnames(data) == id_var] <- "id"
  colnames(data)[colnames(data) == time_var] <- "time"
  colnames(data)[colnames(data) == outcome_var] <- "outcome"
  colnames(data)[colnames(data) == treated_var] <- "treated"
  return(data)
}

create_treatment_groups <- function(data, control_vars, return_merged_df = FALSE) {
  df_treat <- data[data$treated == 1, ]
  df_treat <- df_treat[, c("id", "time")]
  df_treat <- aggregate(time ~ id, data = df_treat, FUN = min)
  colnames(df_treat) <- c("id", "treat_time")
  data <- merge(data, df_treat, by = "id", all.x = T)
  data[is.na(data$treat_time), "treat_time"] <- Inf

  # Check for weakly increasing treatment
  inc <- ifelse(data$treat_time == Inf, 1,
    ifelse(data$time >= data$treat_time & data$treated == 1, 1,
      ifelse(data$time < data$treat_time & data$treated == 0,
        1, 0
      )
    )
  )

  if (!all(as.logical(inc))) {
    stop("Treatment not weakly increasing with time")
  }

  if (length(control_vars) == 0) {
    # Create data.frame of all posible 2x2 estimates. Dyads may appear twice as
    # treatment groups can play the roll of both earlier and later treated
    two_by_twos <- expand.grid(
      unique(data$treat_time),
      unique(data$treat_time)
    )
    colnames(two_by_twos) <- c("treated", "untreated")
    two_by_twos <- two_by_twos[!(two_by_twos$treated ==
      two_by_twos$untreated), ]
    two_by_twos <- two_by_twos[!(two_by_twos$treated == Inf), ]
    # Remove first period
    two_by_twos <- two_by_twos[!(two_by_twos$treated == min(data$time)), ]
    two_by_twos[, c("estimate", "weight")] <- NA
    # Classify estimate "type"
    two_by_twos[, "type"] <- ifelse(two_by_twos$untreated == Inf,
      "Treated vs Untreated",
      ifelse(two_by_twos$untreated ==
        min(data$time),
      "Later vs Always Treated",
      ifelse(two_by_twos$treated >
        two_by_twos$untreated,
      "Later vs Earlier Treated",
      "Earlier vs Later Treated"
      )
      )
    )
  } else if (length(control_vars) > 0) {
    # In the controlled decomposition, each dyad only appears once becasue we
    # do not make the distinction between earlier vs later treated
    two_by_twos <- data.frame()
    for (i in unique(data$treat_time[data$treat_time != Inf])) {
      for (j in unique(data$treat_time)) {
        if (j > i) {
          two_by_twos1 <- data.frame(
            treated = i, untreated = j, weight = NA,
            estimate = NA
          )
          two_by_twos <- rbind(two_by_twos, two_by_twos1)
        }
      }
    }
    two_by_twos[, "type"] <- ifelse(two_by_twos$untreated == Inf,
      "Treated vs Untreated",
      ifelse(two_by_twos$treated ==
        min(data$time) |
        two_by_twos$untreated ==
          min(data$time),
      "Later vs Always Treated",
      "Both Treated"
      )
    )
  }

  # Whether or not to return the merged data too.
  if (return_merged_df == TRUE) {
    return_data <- list(
      "two_by_twos" = two_by_twos,
      "data" = data
    )
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
calculate_weights <- function(data, treated_group, untreated_group) {
  if (untreated_group == 99999) {
    # Treated vs untreated
    n_u <- sum(data$treat_time == untreated_group)
    n_k <- sum(data$treat_time == treated_group)
    n_ku <- n_k / (n_k + n_u)
    D_k <- mean(data[data$treat_time == treated_group, "treated"])
    V_ku <- n_ku * (1 - n_ku) * D_k * (1 - D_k)
    weight1 <- (n_k + n_u)^2 * V_ku
  } else if (treated_group < untreated_group) {
    # early vs late (before late is treated)
    n_k <- sum(data$treat_time == treated_group)
    n_l <- sum(data$treat_time == untreated_group)
    n_kl <- n_k / (n_k + n_l)
    D_k <- mean(data[data$treat_time == treated_group, "treated"])
    D_l <- mean(data[data$treat_time == untreated_group, "treated"])
    V_kl <- n_kl * (1 - n_kl) * (D_k - D_l) / (1 - D_l) * (1 - D_k) / (1 - D_l)
    weight1 <- ((n_k + n_l) * (1 - D_l))^2 * V_kl
  } else if (treated_group > untreated_group) {
    # late vs early (after early is treated)
    n_k <- sum(data$treat_time == untreated_group)
    n_l <- sum(data$treat_time == treated_group)
    n_kl <- n_k / (n_k + n_l)
    D_k <- mean(data[data$treat_time == untreated_group, "treated"])
    D_l <- mean(data[data$treat_time == treated_group, "treated"])
    V_kl <- n_kl * (1 - n_kl) * (D_l / D_k) * (D_k - D_l) / (D_k)
    weight1 <- ((n_k + n_l) * D_k)^2 * V_kl
  }
  return(weight1)
}
scale_weights <- function(two_by_twos) {
  sum_weight <- sum(two_by_twos$weight)
  two_by_twos$weight <- two_by_twos$weight / sum_weight
  return(two_by_twos)
}
run_fwl <- function(data, control_formula) {
  fit_fwl <- fixest::feols(control_formula, data = data)
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
calculate_Omega <- function(data) {
  # TODO Test that within + between = 1
  N <- nrow(data)
  Vd_w <- var(data$d_ikt_til) * (N - 1) / N
  V_d <- var(data$d_it_til) * (N - 1) / N
  Omega <- Vd_w / V_d
  return(Omega)
}
calculate_one_minus_Omega <- function(data) {
  N <- nrow(data)
  Vd_b <- var(data$d_kt_til) * (N - 1) / N
  V_d <- var(data$d_it_til) * (N - 1) / N
  one_minus_Omega <- Vd_b / V_d
  return(one_minus_Omega)
}
calculate_beta_hat_w <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_ikt_til) * (N - 1) / N
  Vd_w <- var(data$d_ikt_til) * (N - 1) / N
  beta_hat_w <- C / Vd_w
  return(beta_hat_w)
}
calculate_beta_hat_b <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_kt_til) * (N - 1) / N
  Vd_b <- var(data$d_kt_til) * (N - 1) / N
  beta_hat_b <- C / Vd_b
  return(beta_hat_b)
}
collapse_x_p <- function(data, control_vars) {
  # Group level Xs
  control_formula <- as.formula(
    paste("~ ", paste(control_vars, collapse = " + "))
  )
  control_data <- data.frame(model.matrix(control_formula, data = data))

  my_ave <- function(x, data) {
    ave(x, data$treat_time, data$time)
  }

  g_control_data <- data.frame(sapply(control_data, my_ave, data))
  colnames(g_control_data) <- paste0("g_", colnames(g_control_data))

  # Group level p
  data$g_p <- ave(data$p, data$treat_time, data$time)
  data <- cbind(data, g_control_data)

  g_control_formula <- as.formula(
    paste("~", paste(colnames(g_control_data), collapse = " + "))
  )

  r_list <- list(data = data, g_control_formula = g_control_formula)

  return(r_list)
}
calc_VD <- function(data) {
  fit_D <- fixest::feols(treated ~ 1 | time + id, data = data)
  data$Dtilde <- fit_D$residuals
  N <- nrow(data)
  VD <- var(data$Dtilde) * (N - 1) / N
  r_list <- list(data = data, VD = VD)
  return(r_list)
}
partial_group_x <- function(data, g_control_formula) {
  g_vars <- unlist(strsplit(as.character(g_control_formula)[2], " \\+ "))
  for (v in g_vars) {
    # Intercept doesn't work with fixest
    if(length(grep("Intercept", v)) > 0) {
      data[, paste0("p_", v)] <- 0
    } else {
      g_formula <- as.formula(paste0(v, " ~ 1 | time + id"))
      fit_g <- fixest::feols(g_formula, data = data)
      data[, paste0("p_", v)] <- fit_g$residuals
    }
  }
  return(data)
}
calc_pgjtile <- function(data, g_control_formula) {
  g_vars <- unlist(strsplit(as.character(g_control_formula)[2], " \\+ "))
  p_g_vars <- paste0("p_", g_vars)
  
  # Don't include p_g_Intercept
  p_g_vars = p_g_vars[-grep("Intercept", p_g_vars)]

  p_g_control_formula <- as.formula(
    paste("Dtilde ~", paste(p_g_vars, collapse = " + "))
  )

  fit_pgj <- fixest::feols(p_g_control_formula, data = data)
  data$pgjtilde <- predict(fit_pgj)
  Rsq <- fixest::r2(fit_pgj, type = "r2")
  r_list <- list(data = data, Rsq = Rsq)
  return(r_list)
}
calc_Vdp <- function(data) {
  N <- nrow(data)
  fit_p <- fixest::feols(g_p ~ 1 | time + id, data = data)
  data$ptilde <- fit_p$residuals
  data$dp <- data$pgjtilde - data$ptilde
  Vdp <- var(data$dp) * (N - 1) / N
  r_list <- list(data = data, Vdp = Vdp)
  return(r_list)
}
calc_BD <- function(data, g_control_formula) {
  g_vars <- all.vars(g_control_formula)
  # Remove intercept
  g_vars = g_vars[-grep("Intercept", g_vars)]
  
  BD_formula <- as.formula(
    paste0("outcome ~ treated + ", paste(g_vars, collapse = " + "), " | id + time")
  )

  fit_BD <- fixest::feols(BD_formula, data = data)
  BD <- fit_BD$coefficients[[1]]
  return(BD)
}
calc_Bb <- function(data) {
  fit_Bb <- fixest::feols(outcome ~ dp | time + id,
    data = data
  )
  Bb <- fit_Bb$coefficients["dp"]
  return(Bb)
}
calculate_beta_hat_d_bkl <- function(Rsq, VD, BD, Vdp, Bb) {
  beta_hat_d_bkl <- ((1 - Rsq) * VD * BD + Vdp * Bb) / ((1 - Rsq) * VD + Vdp)
  return(beta_hat_d_bkl)
}
calculate_s_kl <- function(N, Rsq, VD, Vdp) {
  s_kl <- N^2 * ((1 - Rsq) * VD + Vdp)
  return(s_kl)
}
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
