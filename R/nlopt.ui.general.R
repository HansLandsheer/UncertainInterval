#' Function for the determination of the population thresholds an uncertain and
#' inconclusive interval for test scores with a known common distribution.
#'
#' @param UI.Se (default = .55). Desired sensitivity of the test scores within
#'   the uncertain interval. A value <= .5 is not allowed.
#' @param UI.Sp (default = .55). Desired specificity of the test scores within
#'   the uncertain interval. A value <= .5 is not allowed.
#' @param distribution Name of the continuous distribution, exact as used in R
#'   package stats. Equal to density function minus d. For instance when the
#'   density function is 'dnorm', then the distribution is 'norm'.
#' @param parameters.d0 Named vector of population values or estimates of the
#'   parameters of the distribution of the test scores of the persons without
#'   the targeted condition. For instance \code{c(mean = 0, sd = 1)}. This
#'   distribution should have the lower values.
#' @param parameters.d1 Named vector of population values or estimates of the
#'   parameters of the distribution of the test scores of the persons with the
#'   targeted condition. For instance \code{c(mean = 1, sd = 1)}. The test
#'   scores of d1 should have higher values than d0. If not, use -(test scores).
#'   This distribution should have the higher values.
#' @param overlap.interval A vector with a raw estimate of the lower and upper
#'   relevant of the overlap of the two distributions. If NULL, set to quantile
#'   .001 of the distribution of persons with the targeted condition and
#'   quantile .999 of the distribution of persons without the condition. Please
#'   check whether this is a good estimate of the relevant overlap.
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
#'   the two continuous distributions. The Uncertain Interval is defined as an
#'   interval below and above the intersection of the two distributions, with a
#'   sensitivity and specificity below a desired value (default .55).
#'
#'   Important: The test scores of d1 should have higher values than d0. If not,
#'   use -(test scores). The distribution with parameters d1 should have the
#'   higher test scores.
#'
#'   Important: This is a highly complex function, which is less user friendly
#'   and more error prone than other functions. It is included to show that the
#'   technique also works with continuous distributions other than bi-normal.
#'   Use for bi-normal distributions always ui.binormal.
#'
#'   Only a single intersection is assumed (or a second intersection where the
#'   overlap is negligible).
#'
#'   The function uses an optimization algorithm from the nlopt library
#'
#'   (https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/).
#'
#'   It uses the sequential quadratic programming (SQP) algorithm for
#'   nonlinear constrained gradient-based optimization (supporting both
#'   inequality and equality constraints), based on the implementation by Dieter
#'   Kraft (1988; 1944).
#'
#'   N.B. When a normal distribution is expected, the functions
#'   \code{\link{nlopt.ui}} and \code{\link{ui.binormal}} are recommended.
#'
#' @return List of values: \describe{ \item{$status: }{Integer value with the
#'   status of the optimization (0 is success).} \item{$message: }{More
#'   informative message with the status of the optimization} \item{$results:
#'   }{Vector with the following values:} \itemize{ \item{exp.UI.Sp: }{The
#'   population value of the specificity in the Uncertain Interval, given mu0,
#'   sd0, mu1 and sd1. This value should be very near the supplied value of Sp.}
#'   \item{exp.UI.Se: }{The population value of the sensitivity in the Uncertain
#'   Interval, given mu0, sd0, mu1 and sd1. This value should be very near the
#'   supplied value of UI.Se.} \item{vector of parameter values of distribution
#'   d0, that is, the values that have been supplied in \code{parameters.d0}.}
#'   \item{vector of parameter values of distribution d1, that is, the values
#'   that have been supplied in \code{parameters.d1}.}} \item{$solution:
#'   }{Vector with the following values:} \itemize{ \item{L: }{The population
#'   value of the lower threshold of the uncertain interval.} \item{U: }{The
#'   population value of the upper threshold of the uncertain interval.} } }
#' @references Dieter Kraft, "A software package for sequential quadratic
#'   programming", Technical Report DFVLR-FB 88-28, Institut für Dynamik der
#'   Flugsysteme, Oberpfaffenhofen, July 1988.
#'
#'   Dieter Kraft, "Algorithm 733: TOMP–Fortran modules for optimal control
#'   calculations," ACM Transactions on Mathematical Software, vol. 20, no. 3,
#'   pp. 262-281 (1994).
#'
#'   Landsheer, J. A. (2018). The Clinical Relevance of Methods for Handling
#'   Inconclusive Medical Test Results: Quantification of Uncertainty in Medical
#'   Decision-Making and Screening. Diagnostics, 8(2), 32.
#'   https://doi.org/10.3390/diagnostics8020032
#' @export
#' @importFrom nloptr nloptr
#' @importFrom MASS fitdistr
#' @importFrom car qqPlot
#' @importFrom stats uniroot
#'
#'
#' @examples
#' # A simple test model:
#' nlopt.ui.general(UI.Se = .55, UI.Sp = .55,
#'                  distribution = "norm",
#'                  parameters.d0 = c(mean = 0, sd = 1),
#'                  parameters.d1 = c(mean = 1, sd = 1),
#'                  overlap.interval=c(-2,3))
#' # Standard procedure when using a continuous distribution:
#' nlopt.ui.general(parameters.d0 = c(mean = 0, sd = 1),
#'                  parameters.d1 = c(mean = 1.6, sd = 2))
#'
#' # library(MASS)
#' # library(car)
#' # gamma distributed data
#' set.seed(4)
#' d0 = rgamma(100, shape=2, rate=.5)
#' d1 = rgamma(100, shape=7.5, rate=1)

#' # 1. obtain parameters
#' parameters.d0=MASS::fitdistr(d0, 'gamma')$estimate
#' parameters.d1=MASS::fitdistr(d1, 'gamma')$estimate

#' # 2. test if supposed distributions (gamma) is fitting
#' car::qqPlot(d0, distribution='gamma', shape=parameters.d0['shape'])
#' car::qqPlot(d1, distribution='gamma', shape=parameters.d1['shape'])

#' # 3. draw curves and determine overlap
#' curve(dgamma(x, shape=parameters.d0['shape'], rate=parameters.d0['rate']), from=0, to=16)
#' curve(dgamma(x, shape=parameters.d1['shape'], rate=parameters.d1['rate']), from=0, to=16, add=TRUE)
#' overlap.interval=c(1, 15) # ignore intersection at 0; observe large overlap

#' # 4. get empirical AUC
#' simple_auc(d0, d1)
#' # about .65 --> Poor
#' # .90-1 = excellent (A)
#' # .80-.90 = good (B)
#' # .70-.80 = fair (C)
#' # .60-.70 = poor (D)
#' # .50-.60 = fail (F)

#' # 5. Get uncertain interval
#' (res=nlopt.ui.general (UI.Se = .57,
#'                        UI.Sp = .57,
#'                        distribution = 'gamma',
#'                        parameters.d0 = parameters.d0,
#'                        parameters.d1 = parameters.d1,
#'                        overlap.interval,
#'                        intersection = NULL,
#'                        start = NULL,
#'                        print.level = 0))
#' abline(v=c(res$intersection, res$solution))

#' # 6. Assess improvement when diagnosing outside the uncertain interval
#' sel.d0 = d0 < res$solution[1] |  d0 > res$solution[2]
#' sel.d1 = d1 < res$solution[1] |  d1 > res$solution[2]
#' (percentage.selected.d0 = sum(sel.d0) / length(d0))
#' (percentage.selected.d1 = sum(sel.d1) / length(d1))
#' simple_auc(d0[sel.d0], d1[sel.d1])
#' # AUC for selected scores outside the uncertain interval
#' simple_auc(d0[!sel.d0], d1[!sel.d1])
#' # AUC for deselected scores; worst are deselected


#' # weibull distributed data
#' set.seed(4)
#' d0 = rweibull(100, shape=3, scale=50)
#' d1 = rweibull(100, shape=3, scale=70)

#' # 1. obtain parameters
#' parameters.d0=MASS::fitdistr(d0, 'weibull')$estimate
#' parameters.d1=MASS::fitdistr(d1, 'weibull')$estimate

#' # 2. test if supposed distributions (gamma) is fitting
#' car::qqPlot(d0, distribution='weibull', shape=parameters.d0['shape'])
#' car::qqPlot(d1, distribution='weibull', shape=parameters.d1['shape'])

#' # 3. draw curves and determine overlap
#' curve(dweibull(x, shape=parameters.d0['shape'],
#'       scale=parameters.d0['scale']), from=0, to=150)
#' curve(dweibull(x, shape=parameters.d1['shape'],
#'       scale=parameters.d1['scale']), from=0, to=150, add=TRUE)
#' overlap.interval=c(1, 100) # ignore intersection at 0; observe overlap

#' # 4. get empirical AUC
#' simple_auc(d0, d1)
#' # about .65 --> Poor
#' # .90-1 = excellent (A)
#' # .80-.90 = good (B)
#' # .70-.80 = fair (C)
#' # .60-.70 = poor (D)
#' # .50-.60 = fail (F)

#' # 5. Get uncertain interval
#' (res=nlopt.ui.general (UI.Se = .55,
#'                        UI.Sp = .55,
#'                        distribution = 'weibull',
#'                        parameters.d0 = parameters.d0,
#'                        parameters.d1 = parameters.d1,
#'                        overlap.interval,
#'                        intersection = NULL,
#'                        start = NULL,
#'                        print.level = 0))
#' abline(v=c(res$intersection, res$solution))

#' # 6. Assess improvement when diagnosing outside the uncertain interval
#' sel.d0 = d0 < res$solution[1] |  d0 > res$solution[2]
#' sel.d1 = d1 < res$solution[1] |  d1 > res$solution[2]
#' (percentage.selected.d0 = sum(sel.d0) / length(d0))
#' (percentage.selected.d1 = sum(sel.d1) / length(d1))
#' simple_auc(d0[sel.d0], d1[sel.d1])
#' # AUC for selected scores outside the uncertain interval
#' simple_auc(d0[!sel.d0], d1[!sel.d1])
#' # AUC for deselected scores; these scores are almost indistinguishable

# UI.Se = .55; UI.Sp = .55; distribution = 'norm'; parameters.d0 = c(mean = 0, sd = 1);
# parameters.d1 = c(mean = 1, sd = 1); overlap.interval = NULL; intersection = NULL;
# start = NULL; print.level = 0
nlopt.ui.general <- function(UI.Se = .55,
                     UI.Sp = .55,
                     distribution = 'norm',
                     parameters.d0 = c(mean = 0, sd = 1),
                     parameters.d1 = c(mean = 1, sd = 1),
                     overlap.interval = NULL,
                     intersection = NULL,
                     start = NULL,
                     print.level = 0) {

  uniroot.all.copy <- function (f, interval, lower = min(interval), upper = max(interval), 
            tol = .Machine$double.eps^0.2, maxiter = 1000, n = 100, ...) 
  {
    if (!missing(interval) && length(interval) != 2) 
      stop("'interval' must be a vector of length 2")
    if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
        upper) 
      stop("lower < upper  is not fulfilled")
    xseq <- seq(lower, upper, len = n + 1)
    mod <- f(xseq, ...)
    Equi <- xseq[which(mod == 0)]
    ss <- mod[1:n] * mod[2:(n + 1)]
    ii <- which(ss < 0)
    for (i in ii) Equi <- c(Equi, uniroot(f, lower = xseq[i], 
                                          upper = xseq[i + 1], ...)$root)
    return(Equi)
  }
  if (UI.Se <= .5) stop('Value <= .5 invalid for UI.Se')
  if (UI.Sp <= .5) stop('Value <= .5 invalid for UI.Sp')
  # if (UI.Se > .6) warning('Value > .6 not recommended for UI.Se')
  # if (UI.Sp > .6) warning('Value > .6 not recommended for UI.Sp')

  c01 = UI.Sp / (1 - UI.Sp)
  c11 = UI.Se / (1 - UI.Se)

  d0 = get(paste('d', distribution, sep = ""))
  args=formals(d0)
  m = match(names(parameters.d0), names(args))
  if (any(is.na(m)))
    stop(paste(
      "'parameters.d0' specifies names which are not arguments to", deparse(d0)[1]
    ))
  args[m] = parameters.d0
  formals(d0) <- args # d0(2) ; dnorm(2, 0, 1)

  d1 = get(paste('d', distribution, sep = ""))
  m = match(names(parameters.d1), names(args))
  if (any(is.na(m)))
    stop(paste(
      "'parameters.d1' specifies names which are not arguments to", deparse(d1)[1]
    ))
  args[m] = parameters.d1
  formals(d1) <- args # d1(2) ; dnorm(2, 1, 1)

  p0 = get(paste('p', distribution, sep = ""))
  args=formals(p0)
  m = match(names(parameters.d0), names(args))
  if (any(is.na(m)))
    stop(paste(
      "'parameters.d0' specifies names which are not arguments to", deparse(p0)[1]
    ))
  args[m] = parameters.d0
  formals(p0) <- args # p0(2) ; pnorm(2, 0, 1)

  p1 = get(paste('p', distribution, sep = ""))
  m = match(names(parameters.d1), names(args))
  if (any(is.na(m)))
    stop(paste(
      "'parameters.d1' specifies names which are not arguments to", deparse(p1)[1]
    ))
  args[m] = parameters.d1
  formals(p1) <- args # p1(2) ; pnorm(2, 1, 1)

  q0 = get(paste('q', distribution, sep = ""))
  args=formals(q0)
  m = match(names(parameters.d0), names(args))
  if (any(is.na(m)))
    stop(paste(
      "'parameters.d0' specifies names which are not arguments to", deparse(q0)[1]
    ))
  args[m] = parameters.d0
  formals(q0) <- args # p1(2) ; qnorm(.001, 0, 1)

  q1 = get(paste('q', distribution, sep = ""))
  m = match(names(parameters.d1), names(args))
  if (any(is.na(m)))
    stop(paste(
      "'parameters.d1' specifies names which are not arguments to", deparse(q1)[1]
    ))
  args[m] = parameters.d1
  formals(q1) <- args # p1(2) ; qnorm(.001, 1, 1); q1(q=.001)

  if (print.level > 0) {
    cat('function call d0: ', deparse(d0), '\n' )
    cat('function call d1: ', deparse(d1), '\n' )
    cat('function call p0: ', deparse(p0), '\n' )
    cat('function call p1: ', deparse(p1), '\n' )
    cat('function call q0: ', deparse(q0), '\n' )
    cat('function call q1: ', deparse(q1), '\n' )
  }

  if (is.null(overlap.interval)) overlap.interval = c(q1(.001), q0(.999))

  if (is.null(intersection)) {

    f <- function(x) d0(x)-d1(x)

    # library(rootSolve) # overlap.interval=c(-2, 2)
    I = uniroot.all.copy(f, interval = overlap.interval)
    if (length(I) > 1) {
      d = d0(I)
      I = I[which.max(d)]
      warning('More than one point of intersection. Point with highest density used.')
    }
  } else
    I = intersection # abline(v=I)

  # objective function -(H - L)^2
  eval_f0 <- function(x) {
    return(-(x[2] - x[1]) ^ 2)
    # return(-(x[2] - x[1]) )
  }

  eval_grad_f0 <- function(x) {
    return(c(2 * (x[2] - x[1]), -2 * (x[2] - x[1])))
    # return(c(1,-1))
  }

  # constraint function
  eval_g0 <- function(x) {
    a0 = p0(I) - p0(x[1]) # TN within the uncertain interval
    b0 = p0(x[2]) - p0(I) # FP
    a1 = p1(x[2]) - p1(I) # TP
    b1 = p1(I) - p1(x[1]) # FN

    return(c(a0 - c01 * b0,
             a1 - c11 * b1)) # vector with two constraint values
  }

  # jacobian of constraints
  eval_jac_g0 <- function(x) {
    return(rbind(c(
      -d0(x[1]) ,-c01 * d0(x[2])
    ),
    c(
      c11 * d1(x[1]), d1(x[2])
    )))
  }

  if (is.null(start)){
    par0 = c((overlap.interval[1]+I)/2, (overlap.interval[2]+I)/2)
  } else {par0 = c(start[1], start[2])}

  # eval_g0(par0)
  # eval_jac_g0(par0)

  # res0 <- nloptr(
  #   x0 = par0,
  #   eval_f = eval_f0,
  #   eval_g_ineq = eval_g0,
  #   lb = c(overlap.interval[1], I),
  #   ub = c(I, overlap.interval[2]),
  #   opts = list(
  #     "algorithm" = "NLOPT_LN_COBYLA" ,  # "NLOPT_LN_COBYLA" "NLOPT_LD_SLSQP"
  #     "xtol_rel" = 1.0e-8,
  #     "print_level" = print.level #,
  #     # "check_derivatives" = TRUE,
  #     # "check_derivatives_print" = "all"
  #   )
  # )

  # par0=c(L=I-.5*sd1,H=I+.5*sd0)
  # lb = c(I-2*sd1, I)
  # ub = c(I, I+2*sd0)

  res0 <- nloptr(
    x0 = par0,
    lb = c(overlap.interval[1], I),
    ub = c(I, overlap.interval[2]),
    eval_f = eval_f0,
    eval_g_ineq = eval_g0,
    eval_grad_f = eval_grad_f0,
    eval_jac_g_ineq = eval_jac_g0,
    opts = list(
      "algorithm" = "NLOPT_LD_SLSQP" ,  # "NLOPT_LN_COBYLA" "NLOPT_LD_MMA" "NLOPT_LD_SLSQP"
      "xtol_rel" = 1.0e-8,
      "print_level" = print.level #,
      # "check_derivatives" = TRUE,
      # "check_derivatives_print" = "all"
    )
  )

  TN = p0(I) - p0(res0$solution[1])  # area check UI.Sp: lower area / upper area
  FP = p0(res0$solution[2]) - p0(I)
  TP = p1(res0$solution[2]) - p1(I)  # area check UI.Se: upper area / lower area
  FN = p1(I) - p1(res0$solution[1])
  res = list()
  res$status = res0$status
  res$message = res0$message
  res$intersection = I
  res$results = c(
    exp.UI.Sp = ifelse((TN > 1e-4), TN / (FP + TN), NA),
    exp.UI.Se = ifelse(TP > 1e-4, TP / (FN + TP), NA),
    parameters.d0,
    parameters.d1
  )
  res$solution = c(L = res0$solution[1], U = res0$solution[2])

  return(res)
}

