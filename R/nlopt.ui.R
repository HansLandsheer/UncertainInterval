#' Function for the determination of the population thresholds an uncertain and
#' inconclusive interval for bi-normal distributed test scores.
#'
#' @param UI.Se (default = .55). Desired sensitivity of the test scores within the
#'   uncertain interval. A value <= .5 is not allowed.
#' @param UI.Sp (default = .55). Desired specificity of the test scores within the
#'   uncertain interval. A value <= .5 is not allowed.
#' @param mu0 Population value or estimate of the mean of the test scores of the
#'   persons without the targeted condition.
#' @param sd0 Population value or estimate of the standard deviation of the test
#'   scores of the persons without the targeted condition.
#' @param mu1 Population value or estimate of the mean of the test scores of the
#'   persons with the targeted condition.
#' @param sd1 Population value or estimate of the standard deviation of the test
#'   scores of the persons with the targeted condition.
#' @param intersection Default NULL. If not null, the supplied value is used as
#'   the estimate of the intersection of the two bi-normal distributions.
#'   Otherwise, it is calculated.
#' @param start Default NULL. If not null, the first two values of the supplied
#'   vector are used as the starting values for the \code{nloptr} optimization
#'   function.
#' @param print.level Default is 0. The option print.level controls how much
#'   output is shown during the optimization process. Possible values: 0)
#'   (default)	no output; 1)	show iteration number and value of objective
#'   function; 2)	1 + show value of (in)equalities; 3)	2 + show value of
#'   controls.
#'
#' @details The function can be used to determinate the uncertain interval of
#'   two bi-normal distributions. The Uncertain Interval is defined as an
#'   interval below and above the intersection of the two distributions, with a
#'   sensitivity and specificity below a desired value (default .55).
#'
#'   Only a single intersection is assumed (or a second intersection where the
#'   overlap is negligible).
#'
#'   The function uses an optimization algorithm from the nlopt library
#'   (https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/): the sequential
#'   quadratic programming (SQP) algorithm for nonlinearly constrained
#'   gradient-based optimization (supporting both inequality and equality
#'   constraints), based on the implementation by Dieter Kraft (1988; 1944).
#'
#' @return List of values: \describe{ \item{$status: }{Integer value with the
#'   status of the optimization (0 is success).} \item{$message: }{More
#'   informative message with the status of the optimization} \item{$results:
#'   }{Vector with the following values:} \itemize{ \item{exp.UI.Sp: }{The
#'   population value of the specificity in the Uncertain Interval, given mu0,
#'   sd0, mu1 and sd1. This value should be very near the supplied value of Sp.}
#'   \item{exp.UI.Se: }{The population value of the sensitivity in the Uncertain
#'   Interval, given mu0, sd0, mu1 and sd1. This value should be very near the
#'   supplied value of UI.Se.} \item{mu0: }{The value that has been supplied for
#'   mu0.} \item{sd0: }{The value that has been supplied for sd0.} \item{mu1:
#'   }{The value that has been supplied for mu1.} \item{sd1: }{The value that
#'   has been supplied for sd1.} } \item{$solution: }{Vector with the following
#'   values:} \itemize{ \item{L: }{The population value of the lower threshold
#'   of the Uncertain Interval.} \item{U: }{The population value of the upper
#'   threshold of the Uncertain Interval.} } }
#' @references Dieter Kraft, "A software package for sequential quadratic
#'   programming", Technical Report DFVLR-FB 88-28, Institut für Dynamik der
#'   Flugsysteme, Oberpfaffenhofen, July 1988.
#'
#'   Dieter Kraft, "Algorithm 733: TOMP–Fortran modules for optimal control
#'   calculations," ACM Transactions on Mathematical Software, vol. 20, no. 3,
#'   pp. 262-281 (1994).

#' @export
#' @importFrom rootSolve uniroot.all
#' @importFrom nloptr nloptr
#' @importFrom stats dnorm pnorm sd
#'
#' @examples
#' # A simple test model:
#' nlopt.ui()
#' # Using another bi-normal distribution:
#' nlopt.ui(mu0=0, sd0=1, mu1=1.6, sd1=2)
#'
# UI.Se = .55; UI.Sp = .55; mu0 = 0; sd0 = 1; mu1 = 1; sd1 = 1; intersection = NULL; start=NULL; print.level=0
nlopt.ui <- function(UI.Se = .55, UI.Sp = .55,
                     mu0 = 0, sd0 = 1,
                     mu1 = 1, sd1 = 1,
                     intersection = NULL,
                     start=NULL, print.level=0) {
  if (UI.Se <= .5) stop('Value <= .5 invalid for UI.Se')
  if (UI.Sp <= .5) stop('Value <= .5 invalid for UI.Sp')
  # if (UI.Se > .6) warning('Value > .6 not recommended for UI.Se')
  # if (UI.Sp > .6) warning('Value > .6 not recommended for UI.Sp')

  c01 = UI.Sp / (1 - UI.Sp)
  c11 = UI.Se / (1 - UI.Se)
  if (is.null(intersection)) {
    intersect.binormal1 <- function(mu0, sd0, mu1, sd1) {
        if (sd0==sd1){
          is <- (mu1+mu0)/2
        }else{
          B <- (mu0 / sd0 ^ 2 - mu1 / sd1 ^ 2)
          A <- 0.5 * (1 / sd1 ^ 2 - 1 / sd0 ^ 2)
          C <-
            0.5 * (mu1 ^ 2 / sd1 ^ 2 - mu0 ^ 2 / sd0 ^ 2) - log(sd0 / sd1)

          is = (-B + c(1, -1) * sqrt(B ^ 2 - 4 * A * C)) / (2 * A)
        }
        d = dnorm(is, mu0, sd0) + dnorm(is, mu1, sd1)
        is[order(d)] # tail has highest density
      }

      intersection = intersect.binormal1(mu0, sd0, mu1, sd1)
      if (length(intersection) > 1) {
        intersection=tail(intersection, n=1)
        warning('More than one point of intersection. Highest used.')
      }
  }

  # objective function -(H - L)^2
  eval_f0 <- function(x) {
    return(-(x[2] - x[1]) ^ 2)
    # return(-(x[2] - x[1]) )
  }

  eval_grad_f0 <- function(x) {
    return(c(2 * (x[2] - x[1]),-2 * (x[2] - x[1])))
    # return(c(1,-1))
  }

  # constraint function
  eval_g0 <- function(x) {
    a0 = pnorm(intersection, mu0, sd0) - pnorm(x[1], mu0, sd0) # TN within the uncertain interval
    b0 = pnorm(x[2], mu0, sd0) - pnorm(intersection, mu0, sd0) # FP
    a1 = pnorm(x[2], mu1, sd1) - pnorm(intersection, mu1, sd1) # TP
    b1 = pnorm(intersection, mu1, sd1) - pnorm(x[1], mu1, sd1) # FN

    return(c(a0 - c01 * b0,
             a1 - c11 * b1)) # vector with two constraint values
  }

  # jacobian of constraints x=par0
  eval_jac_g0 <- function(x) {
    return(rbind(c(
      -dnorm(x[1], mu0, sd0) , -c01 * dnorm(x[2], mu0, sd0)
    ),
    c(
      c11 * dnorm(x[1], mu1, sd1), dnorm(x[2], mu1, sd1)
    )))
  }

  if (is.null(start)) par0=c(L=intersection-.5*sd1,
                             H=intersection+.5*sd0) else par0=c(L=start[1], H=start[2])

  res0 <- nloptr( x0=par0,
                  eval_f=eval_f0,
                  eval_grad_f = eval_grad_f0,
                  lb = c(intersection-2*sd1, intersection),
                  ub = c(intersection, intersection+2*sd0),
                  eval_g_ineq = eval_g0,
                  eval_jac_g_ineq = eval_jac_g0,
                  # eval_g_eq = eval_g0,
                  # eval_jac_g_eq = eval_jac_g0,
                  opts = list("algorithm"="NLOPT_LD_SLSQP",
                              "xtol_rel"=1.0e-8,
                              print_level=print.level #,
                              # "check_derivatives" = TRUE,
                              # "check_derivatives_print" = "all"
                              )
  )

  # res0 <- nloptr(
  #   x0 = c(intersection - sd1, intersection + sd0),
  #   eval_f = eval_f0,
  #   eval_grad_f = eval_grad_f0,
  #   eval_g_ineq = eval_g0,
  #   eval_jac_g_ineq = eval_jac_g0,
  #   lb = c(intersection-2*sd1, intersection),
  #   ub = c(intersection, intersection+2*sd0),
  #   opts = list(
  #     "algorithm" = "NLOPT_LD_MMA",
  #     "xtol_rel" = 1.0e-8,
  #     "print_level" = 2 #,
  #     # "check_derivatives" = TRUE,
  #     # "check_derivatives_print" = "all"
  #     )
  #   )

  TN = pnorm(intersection, mu0, sd0) - pnorm(res0$solution[1], mu0, sd0)  # area check UI.Sp: lower area / upper area
  FP = pnorm(res0$solution[2], mu0, sd0) - pnorm(intersection, mu0, sd0)
  TP = pnorm(res0$solution[2], mu1, sd1) - pnorm(intersection, mu1, sd1)  # area check UI.Se: upper area / lower area
  FN = pnorm(intersection, mu1, sd1) - pnorm(res0$solution[1], mu1, sd1)
  res = list()
  res$status = res0$status
  res$message = res0$message
  res$intersection = intersection
  res$results = c(exp.UI.Sp = ifelse((TN > 1e-4),TN / (FP + TN), NA),
                  exp.UI.Se = ifelse(TP > 1e-4, TP / (FN + TP), NA),
                  mu0=unname(mu0), sd0=unname(sd0), mu1=unname(mu1), sd1=unname(sd1))
  res$solution = c(L = res0$solution[1], U = res0$solution[2])

  return(res)
}

