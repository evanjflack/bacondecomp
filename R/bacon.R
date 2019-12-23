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
#' @export
bacon <- function(formula,
                  data,
                  id_var,
                  time_var) {
  
  # Unpack variable names and rename variables
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
  
  # Create 2x2 grid of treatment groups
  treatment_group_calc <- create_treatment_groups(data, control_vars, 
                                                  return_merged_df = TRUE)
  two_by_twos <- treatment_group_calc$two_by_twos
  data <- treatment_group_calc$data
  
  # First period in the panel
  first_period <- min(data$time)
  
  # Uncontrolled ---------------------------------------------------------------
  if (length(control_vars) == 0) {
    for (i in 1:nrow(two_by_twos)) {
      treated_group <- two_by_twos[i, "treated"]
      untreated_group <- two_by_twos[i, "untreated"]
      
      data1 <- subset_data(data, treated_group, untreated_group)
      
      
      
    }
    
    two_by_twos <- scale_weights(two_by_twos)
    
    return(two_by_twos)
    
  } else if (length(control_vars) > 0) {
    # Controled ----------------------------------------------------------------
    # Predict Treatment and calulate demeaned residuals
    control_formula <- update(formula,
                              paste0("treated ~ . + factor(time) + factor(id) -",
                                     treated_var))
    
    data <- calculate_dps(data, control_formula)
    Sigma <- calculate_Sigma(data)
    beta_hat_w <- calculate_beta_hat_w(data)
    
    N <- nrow(data)
    Vd_b <- var(data$d_kt_til)*(N - 1)/N
    
    data <- run_fwl(data, control_formula)
    data <- collapse_x_p(data, control_vars)
    
    for (i in 1:nrow(two_by_twos)) {
      treated_group <- two_by_twos[i, "treated"]
      untreated_group <- two_by_twos[i, "untreated"]
      data1 <- data[data$treat_time %in% c(treated_group, untreated_group), ]
      
      r_calc_VD <- calc_VD(data1)
      VD <- r_calc_VD$VD
      data1 <- r_calc_VD$data
      
      data1 <- partial_group_x(data1, control_vars)
      
      r_calc_pgjtilde <- calc_pgjtile(data1)
      Rsq <- r_calc_pgjtilde$Rsq
      data1 <- r_calc_pgjtilde$data
      
      r_calc_Vdp <- calc_Vdp(data1)
      Vdp <- r_calc_Vdp$Vdp
      data1 <- r_calc_Vdp$data
      
      BD <- calc_BD(data1)
      Bb <- calc_Bb(data1)
      
      # weight
      finals <- N^2*((1 - Rsq)*VD + Vdp)
      Beta <- ((1 - Rsq)*VD*BD + Vdp*Bb)/((1- Rsq)*VD + Vdp)
      
    
      two_by_twos[i, "weight"] <- finals
      two_by_twos[i, "estimate"] <- Beta
    }
    
    # two_by_twos <- scale_weights(two_by_twos)
    
    r_list <- list("beta_hat_w" = beta_hat_w,
                   "Sigma" = Sigma,
                   "two_by_twos" = two_by_twos)
    return(r_list)
  }
}


# Frisch-Waugh-Locell Regression
# Predict Treatment and calulate demeaned residuals
run_fwl <- function(data, control_formula) {
  fit_fwl <- lm(control_formula, data = data)
  data$p <- predict(fit_fwl)
  data$d <- fit_fwl$residuals
  return(data)
}

collapse_x_p <- function(data, control_vars) {
  # Collapse X's and p to group/year level
  for (v in c("p", control_vars)) {
    print(v)
    data[, paste0("g_", v)] <-  ave(data[, v], data$treat_time, data$time)
  }
  return(data)
}

calc_VD <- function(data) {
  fit_D <- lm(treated ~ factor(id) + factor(time), data = data)
  data$Dtilde <- fit_D$residuals
  N <- nrow(data)
  VD <- var(data$Dtilde)*(N - 1)/N
  r_list <- list(data = data, VD = VD)
  return(r_list)
}

partial_group_x <- function(data, control_vars) {
  g_vars <- paste0("g_", control_vars)
  for (v in g_vars) {
    g_formula <- as.formula(paste0(v, " ~ factor(id) + factor(time)"))
    fit_g <- lm(g_formula, data = data)
    data[, paste0("p_", v)] <- fit_g$residuals
  }
  return(data)
}

calc_pgjtile <- function(data) {
  fit_pgj <- lm(Dtilde ~ g_l_income + g_l_pop, data = data)
  data$pgjtilde <- predict(fit_pgj)
  Rsq <- summary(fit_pgj)$r.squared
  r_list <- list(data = data, Rsq = Rsq)
  return(r_list)
}

calc_Vdp <- function(data) {
  fit_p <- lm(g_p ~ factor(id) + factor(time), data = data)
  data$ptilde <- predict(fit_p)
  data$dp <- data$pgjtilde - data$ptilde
  Vdp <- var(data$dp)*(N - 1)/N
  r_list <- list(data = data, Vdp = Vdp)
  return(r_list)
}

calc_BD <- function(data) {
  fit_BD <- lm(outcome ~ treated + g_l_income + g_l_pop + factor(time) + factor(id), 
               data = data)
  BD <- fit_BD$coefficients["treated"]
  return(BD)
}

calc_Bb <- function(data) {
  fit_Bb <- lm(outcome ~ dp + factor(time) + factor(id), 
               data = data)
  Bb <- fit_Bb$coefficients["dp"]
  return(Bb)
}

#' Unpack Variable Names from Formula
#' 
#' @param formula 
#' 
#' @return a list with 3 elements: outcome_var, treated_var, control_vars
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
#' @param data a data.frame
#' @param id_var character
#' @param time_var character
#' @param outcome_var character
#' @param treated_var character
#' 
#' @return data.frame with renmaed columns
rename_vars <- function(data, id_var, time_var, outcome_var, treated_var) {
  colnames(data)[colnames(data) == id_var] <- "id"
  colnames(data)[colnames(data) == time_var] <- "time"
  colnames(data)[colnames(data) == outcome_var] <- "outcome" 
  colnames(data)[colnames(data) == treated_var] <- "treated"
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
create_treatment_groups <- function(data, control_vars, return_merged_df = FALSE){
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
  
  if (length(control_vars) == 0) {
    # create data.frame of all posible 2x2 estimates
    two_by_twos <- expand.grid(unique(data$treat_time),
                               unique(data$treat_time))
    colnames(two_by_twos) <- c("treated", "untreated")
    two_by_twos <- two_by_twos[!(two_by_twos$treated == two_by_twos$untreated), ]
    two_by_twos <- two_by_twos[!(two_by_twos$treated == 99999), ]
    # Remove first period
    two_by_twos <- two_by_twos[!(two_by_twos$treated == min(data$time)), ]
    two_by_twos[, c("estimate", "weight")] <- 0
  } else if (length(control_vars) > 0) {
    two_by_twos <- data.frame()
    for (i in unique(data$treat_time[data$treat_time != 99999])) {
      for (j in unique(data$treat_time)) {
        if (j > i) {
          two_by_twos1 <- data.frame(treated = i, untreated = j, weight = 0, estimate = 0)
          two_by_twos <- rbind(two_by_twos, two_by_twos1)
        }
      }
    }
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
#' @param data a data.frame
#' @param treated_group integer
#' @param untreated_group interger
#' 
#' @return subsetted data.frame
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

#' Scale 2x2 Weights
#' 
#' @param two_by_twos
#' 
#' @return two_by_twos
scale_weights <- function(two_by_twos) {
  sum_weight <- sum(two_by_twos$weight)
  two_by_twos$weight <- two_by_twos$weight/sum_weight
  return(two_by_twos)
}

#' Calculate d_it and Variants
#' 
#' @param data 
#' @param control_formula
#' 
#' @return data frame with new "p" variables
calculate_dps <- function(data, control_formula) {
  
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
  
  # Demean the 
  data$D_it_til <- data$treated - ave(data$treated, data$id) - ave(data$treated, data$time) + mean(data$treated)
  data$D_jt_til <- data$D_it_til - ave(data$D_it_til, data$treat_time)
  
  data$p_it_til <- predict(fit_treat)
  data$p_jt_til <- data$p_it_til - ave(data$p_it_til, data$treat_time)
  
  return(data)
}


#' Calculate Sigma
#' 
#' @param data
#' 
#' @return Sigma
calculate_Sigma <- function(data) {
  # TODO Test that within + between = 1
  N <- nrow(data)
  Vd_w <- var(data$d_ikt_til)*(N - 1)/N
  V_d <- var(data$d_it_til)*(N - 1)/N
  Sigma <- Vd_w/V_d
  return(Sigma)
}

#' Calculate 1 - Sigma
#' 
#' @param data 
#' 
#' @return one_minus_Sigma
calculate_one_minus_Sigma <- function(data) {
  N <- nrow(data)
  Vd_b <- var(data$d_kt_til)*(N - 1)/N
  V_d <- var(data$d_it_til)*(N - 1)/N
  one_minus_Sigma <- Vd_b/V_d
  return(one_minus_Sigma)
}

#' Calculate beta_hat_w
calculate_beta_hat_w <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_ikt_til)*(N - 1)/N
  Vd_w <- var(data$d_ikt_til)*(N - 1)/N
  beta_hat_w <- C/Vd_w
  return(beta_hat_w)
}

#' Calculate between beta hat
#' 
#' @param data
#' 
#' @return beta_hat_b
calculate_beta_hat_b <- function(data) {
  N <- nrow(data)
  C <- cov(data$outcome, data$d_kt_til)*(N - 1)/N
  Vd_b <- var(data$d_kt_til)*(N - 1)/N
  beta_hat_b <- C/Vd_b
  return(beta_hat_b)
}

#' Calculate Weights in Controlled Decomposition
#' 
#' @param data
#' @param treated_group integer
#' @param untreated_group integer
#' @param V_db numberurn s_kl
calculate_weights_controled <- function(data, treated_group,
                                        untreated_group, Vd_b) {
  # TODO test: between 0-1? (ask him), and sum to 1
  # Treated vs untreated
  n_k <- sum(data$treat_time == treated_group)
  n_l <- sum(data$treat_time == untreated_group)
  N <- nrow(data)
  Vd_bkl <- var(data$d_kt_til)*(N - 1)/N # this is the V(d_kt) but only for these subgroups
  s_kl <- (n_k + n_l)^2 * Vd_bkl/Vd_b
  return(s_kl)
}

calculate_beta_hat_d_bkl <- function(data) {
  N <- nrow(data)
  Vd_bkl <- var(data$d_kt_til)*(N - 1)/N
  C1 <- cov(data$outcome, data$D_jt_til - data$p_jt_til)
  beta_hat_d_bkl <- C1/Vd_bkl
  
  return(beta_hat_d_bkl)
}

# calculate_VD_kl <- function(data, treated_group, untreated_group) {
#   n_k <- sum(data$treat_time == treated_group)
#   n_l <- sum(data$treat_time == untreated_group)
#   
#   D_bar_k <- mean(data[data$treat_time == treated_group, "treated"])
#   D_bar_l <- mean(data[data$treat_time == untreated_group, "treated"])
#   
#   sum_val <- 0
#   for (t in unique(data$time)) {
#     D_bar_kt <- mean(data[data$treat_time == treated_group & data$time == t, "treated"])
#     D_bar_lt <- mean(data[data$treat_time == untreated_group & data$time == t, "treated"])
#     val <- ((D_bar_kt - D_bar_k) - (D_bar_lt - D_bar_l))^2
#     sum_val = sum_val + val
#   }
#   capT <- length(unique(data$time))
#   adj <- (n_k*n_l)/((n_k + n_l)^2)
#   VD_kl <- adj*(1/capT)*sum_val
#   return(VD_kl)
# }

# "Regression" way
# calculate_VD_kl <- function(data) {
#   data$D_jt_til <- data$D_it - ave(data$D_it_til, data$treat_time) - ave(data$D_it_til, data$time)
#   fit <- lm(D_jt_til ~ factor(id) + factor(treat_time), data = data)
#   resid <- fit$residuals
#   N <- nrow(data)
#   VD_kl <- var(resid)*(N - 1)/N
#   return(VD_kl)
# }

# calculate_Vp_bkl <- function(data, treated_group, untreated_group) {
#   n_k <- sum(data$treat_time == treated_group)
#   n_l <- sum(data$treat_time == untreated_group)
#   
#   p_bar_k <- mean(data[data$treat_time == treated_group, "p_it"])
#   p_bar_l <- mean(data[data$treat_time == untreated_group, "p_it"])
#   
#   sum_val <- 0
#   for (t in unique(data$time)) {
#     p_bar_kt <- mean(data[data$treat_time == treated_group & data$time == t, "p_it"])
#     p_bar_lt <- mean(data[data$treat_time == untreated_group & data$time == t, "p_it"])
#     val <- ((p_bar_kt - p_bar_k) - (p_bar_lt - p_bar_l))^2
#     sum_val = sum_val + val
#   }
#   capT <- length(unique(data$time))
#   adj <- (n_k*n_l)/((n_k + n_l)^2)
#   Vp_bkl <- adj*(1/capT)*sum_val
#   return(Vp_bkl)
# }

# "by hand" way
# calculate_Vp_bkl <- function(data) {
#   fit <- lm(data$p_jt_til ~ factor(id) + factor(treat_time), 
#             data = data)
#   resid <- fit$residuals
#   N <- nrow(data)
#   Vp_bkl <- var(resid)*(N - 1)/N
# }

# calculate_beta_hat_22_kl <- function(data) {
#   fit <- lm(outcome ~ treated + factor(id) + factor(time), data = data)
#   beta_hat_22_kl <- fit$coefficients["treated"]
#   return(beta_hat_22_kl)
# }

# calculate_beta_hat_p_bkl <- function(data) {
#   data$y_jt_bar <- ave(data$outcome, data$treat_time, data$time)
#   fit <- lm(y_jt_bar ~ p_jt_bar + factor(id) + factor(time), data = data)
#   beta_hat_p_bkl <- fit$coefficients["p_jt_bar"]
#   return(beta_hat_p_bkl)
# }

# calculate_beta_hat_d_bkl <- function(data, treated_group, untreated_group) {
#   VD_kl <- calculate_VD_kl(data, treated_group, untreated_group)
#   beta_hat_22_kl <- calculate_beta_hat_22_kl(data)
#   Vp_bkl <- calculate_Vp_bkl(data, treated_group, untreated_group)
#   beta_hat_p_bkl <- calculate_beta_hat_p_bkl(data)
#   
#   N <- nrow(data)
#   Vd_bkl <- var(data$d_kt_til)*(N - 1)/N
#   beta_hat_d_bkl <- (VD_kl*beta_hat_22_kl - Vp_bkl*beta_hat_p_bkl)/Vd_bkl
# }