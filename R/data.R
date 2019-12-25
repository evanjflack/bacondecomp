#' Aggregated data from Goodman (In Press)
#'
#' A data set containing state/year level data on an educational reform and
#'  future income. This is an aggregated version of the data used by Goodman
#'  (2019, JOLE) to estimate the effect of compulsory high school math
#'  coursework on future earnings.
#'
#' @format A data.frame with 520 observations and 5 variables
#' \describe{
#'   \item{state}{The state (unit of analysis).}
#'   \item{class}{The high school class (time).}
#'   \item{reform_math}{Indicator for whether the reform was in place for the state/class.}
#'   \item{reformyr_math}{The year the math reform was first implemented for the state. Set to NA if never implemented.}
#'   \item{incearn_ln}{Natural log of future income.}
#' }
"math_reform"

#' Data from Cheng and Hoekstra (2013, JHR)
#'
#' @format A data.frame with 520 observations and 159 variables
#' \describe{
#'   \item{st}{The state (unit of analysis).}
#'   \item{year}{Calendar year (time).}
#'   \item{l_homicide}{Log of state/year homicide rate}
#'   \item{post}{Indicator whether castle reform has been implemented}
#' }
"castle"


#' Data from Stevenson and Wolfers (2006, QJE)
#'
#' @format A data.frame with 3366 observations and 147 variables
"divorce"
