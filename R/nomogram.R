# prob.pre.test=.3; LR=c(44, .022); probs.post.test=c(pos=NULL, neg=NULL)
# prob.pre.test=.3; probs.post.test=c(0.9778, 1-0.9785)); LR=c(NULL, NULL);


#' Fagan's nomogram to show the relationships between the prior probability, the
#' likelihood ratios, sensitivity and specificity, and the posterior
#' probability.
#'
#' Next to plotting Fagan's nomogram, this function also calculates the
#' minimally needed values for specificity and sensitivity to reach desired
#' posttest probabilities (or likelihood ratios) for a grey zone (Coste et al.,
#' 2003, 2006).
#'
#' @param prob.pre.test The prior test probability, with a default value of .5.
#'   Often, (local) prevalence is used.
#' @param probs.post.test A vector of two values that give the desired posttest
#'   probabilities of observing the event in the case of a positive test result
#'   (positive posttest probability: pos), and the posttest probability of
#'   observing the event in the case of a negative test result (negative
#'   posttest probability: neg). When not given, these probabilities are
#'   calculated using the likelihood ratios (LR).
#' @param SeSp A vector of two values that give the desired sensitivity and
#'   specificity. When not given, the Se and Sp values are calculated from the
#'   desired posttest probabilities.
#' @param LR A vector of two values that give the positive likelihood ratio
#'   (sensitivity / (1- specificity)): PLR of observing the event, and the
#'   negative likelihood ratio ((1 - sensitivity) / specificity): NLR of not
#'   observing the event. PLR is a value > 1, NLR is a value between 0 and 1.
#'   When not given, the LR values are calculated from the desired posttest
#'   probabilities.
#' @param plot A Boolean that indicates whether a plot is desired.
#'
#' @return Vector of values: \describe{ \item{$pre: }{The given pre-test
#'   probability.} \item{$min.LRpos: }{The given or calculated minimally
#'   required positive likelihood ratio. If no value is provided, it is
#'   calculated.} \item{$max.LRneg: }{The given or calculated maximally required
#'   negative likelihood ratio. If no value is provided, it is calculated.}
#'   \item{$post.pos: }{The given or calculated positive posttest probability.}
#'   \item{$minSp: }{The minimum value for the specificity, needed to reach the
#'   desired posttest probabilities.} \item{$minSe: }{The minimum value for the
#'   sensitivity, needed to reach the desired posttest probabilities.} }
#' @details Parameter probs.post.test or SeSp or LR must be supplied, the other
#'   two values are calculated. When more than one parameter is given the other
#'   two are ignored. The basis of this function is adapted from package
#'   TeachingDemos.
#' @references { Fagan, T. J. (1975). Nomogram for Bayes theorem. The New
#'   England Journal of Medicine, 293(5), 257-257.
#'
#'   Coste, J., Jourdain, P., & Pouchot, J. (2006). A gray zone assigned to
#'   inconclusive results of quantitative diagnostic tests: application to the
#'   use of brain natriuretic peptide for diagnosis of heart failure in acute
#'   dyspneic patients. Clinical Chemistry, 52(12), 2229-2235.
#'
#'   Coste, J., & Pouchot, J. (2003). A grey zone for quantitative diagnostic
#'   and screening tests. International Journal of Epidemiology, 32(2),
#'   304-313. }
#' @importFrom graphics axis segments text title

#' @export
#'
#' @examples
#' # Show calculated results (first 3 times about the same)
#' (nomogram(prob.pre.test = .10, probs.post.test=c(pos=.70, neg=.001), plot=FALSE))
#' (nomogram(prob.pre.test = .10, SeSp=c(Se=0.991416309, Sp=0.952789700), plot=FALSE))
#' (nomogram(prob.pre.test = .10, LR=c(pos=21, neg=0.0090090091), plot=FALSE))
#' (nomogram(prob.pre.test = .10, SeSp=c(Se=0.99, Sp=0.95), plot=FALSE))
#' # plot only
#' nomogram(prob.pre.test = .10, LR=c(pos=21, neg=0.0090090091))
#' # plot and display precise results
#' (nomogram(prob.pre.test = .10, probs.post.test=c(pos=.70, neg=.001)))
#'
#' # check the influence of different values of prevalence
#' i=1
#' out=matrix(0,nrow = 9, ncol= 7)
#' for (prev in (seq(.1, .9, by=.1))) {
#'   out[i,]=nomogram(prob.pre.test=prev, probs.post.test=c(.95, .05), plot=FALSE)
#'   i=i+1
#' }
#' colnames(out) = names(nomogram(prob.pre.test=prev, probs.post.test=c(.95, .05), plot=FALSE))
#' out
#'

# probs.post.test=c(pos=NULL, neg=.05)
# SeSp = c(.5, .5); prob.pre.test=.5; probs.post.test=c(pos=NULL, neg=NULL);
# LR=c(PLR=NULL, NLR=NULL); plot=FALSE

nomogram <- function(prob.pre.test=.5, probs.post.test=c(pos=NULL, neg=NULL),
                        SeSp=c(Se=NULL, Sp=NULL), LR=c(PLR=NULL, NLR=NULL), plot=T) {

  if (!is.null(probs.post.test[1])) if (probs.post.test[1] < prob.pre.test)
    stop('Invalid assumption. Probs.post.test[1] (pos) must be equal or larger than prob.pre.test.')
  if (!is.null(probs.post.test[2])) if (probs.post.test[2] > prob.pre.test)
    stop('Invalid assumption. Probs.post.test[2] (neg) must be equal or smaller than prob.pre.test.')

  logits <- function(p) log(p/(1 - p))

  # LR = c(LR+, LR-)
  # calc.LR2 = function(prob.pre.test, probs.post.test){
  #   exp(logits(probs.post.test)-logits(prob.pre.test))
  # }
  calc.LR = function(prob.pre.test, probs.post.test){
    (probs.post.test * (1-prob.pre.test)) /
     ((1-probs.post.test) *prob.pre.test)
  }
  # prob.pre.test=.1; probs.post.test=.001;
  # calc.LR(prob.pre.test, probs.post.test);
  # calc.LR2(prob.pre.test, probs.post.test)

  calc.Sp <- function(LR)  (1-LR[1])/(LR[2]-LR[1])
  calc.Se <- function(LR) LR[1]*(1-(1-LR[1])/(LR[2]-LR[1]))

  logits.pre <- logits(prob.pre.test)
  Se = SeSp[1]; Sp=SeSp[2]
  if (!is.null(probs.post.test)){
    LR = calc.LR(prob.pre.test, probs.post.test) # LR positive
    logits.post <- log(LR) + logits.pre
    Se = calc.Se(LR); Sp = calc.Sp(LR)
  } else if (!is.null(SeSp)){
    LR = c(PLR=SeSp[1] / (1- SeSp[2]), NLR=(1 - SeSp[1]) / SeSp[2])
    logits.post <- log(LR) + logits.pre
    probs.post.test <- exp(logits.post)/(1 + exp(logits.post))
  } else  if (!is.null(LR)){
    logits.post <- log(LR) + logits.pre
    probs.post.test <- exp(logits.post)/(1 + exp(logits.post))
    Se = calc.Se(LR); Sp = calc.Sp(LR)
  } else stop("Parameter probs.post.test or SeSP or LR must be supplied.")

  if (any(prob.pre.test > 1) | any(prob.pre.test < 0) |
      any(probs.post.test > 1) | any(probs.post.test < 0) |
      any(LR < 0) | any(is.infinite(LR)) | any(is.nan(LR))) {
      stop("Function call has wrong parameter values.")
  }

  # if (LR[1]!=LR[2]) {Se = calc.Se(LR); Sp = calc.Sp(LR) }

  compl.logit.pre <- logits(1 - prob.pre.test)

    if (plot) {
      opar <- par(no.readonly = T)
      on.exit(par(opar))
      par(mar = c(1.5, 6, 2, 6))
      LR.vec <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2,
                  0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
      prob.vec <- c(0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02,
                    0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                    0.8, 0.9, 0.93, 0.95, 0.97, 0.98, 0.99, 0.993, 0.995,
                    0.997, 0.998, 0.999)
      plot(0, 0, type = "n", ylim = range(logits(prob.vec)), axes = FALSE,
           xlab = "", ylab = "")
      axis(2, rev(logits(prob.vec)), 100 * prob.vec, pos = -1,
           las = 1, cex.axis = 0.7)
      axis(2, rev(logits(prob.vec)), 100 * prob.vec, pos = -1,
           tck = 0.03, labels = FALSE)
      axis(4, logits(prob.vec), 100 * prob.vec, pos = 1, las = 1,
           cex.axis = 0.7)
      axis(4, logits(prob.vec), 100 * prob.vec, pos = 1, tck = 0.03,
           labels = FALSE)
      axis(2, log(LR.vec[1:10])/2, LR.vec[1:10], pos = 0, las = 1,
           cex.axis = 0.7)
      axis(2, log(LR.vec[1:10])/2, LR.vec[1:10], pos = 0, tck = 0.03,
           labels = FALSE)
      axis(4, log(LR.vec[10:19])/2, LR.vec[10:19], pos = 0, las = 1,
           cex.axis = 0.7)
      axis(4, log(LR.vec[10:19])/2, LR.vec[10:19], pos = 0, tck = 0.03,
           labels = FALSE)
      text(0, 4.5, "Likelihood ratio", cex = 1.2)
      segments(-1, compl.logit.pre, 1, logits.post, lwd = 1.5, col = c(3, 2))
      mtext(side = 2, text = "Pretest %", line = 2,
            cex = 1.2)
      mtext(side = 4, text = "Posttest %", line = 2,
            cex = 1.2, las = 3)
      title(main = "Fagan's nomogram")
      text(0, -5.3, paste("Pre test % of disease =",
                          round(100 * prob.pre.test, 2), "% \n",
                          "Likelihood ratios (+;-)=",
                          round(LR[1], 4),';', round(LR[2], 4),"\n", "posttest %. of disease (+;-) =",
                          round(100 * probs.post.test[1], 2), ';', round(100 * probs.post.test[2], 4),"%\n",
                          'Min Se =', round(Se, 2), 'Min Sp =', round(Sp, 2)),
           cex = 0.7) }
    invisible(c(pre = unname(prob.pre.test),
                post.pos= unname(probs.post.test[1]), post.neg = unname(probs.post.test[2]),
                min.LRpos = unname(LR[1]), max.LRneg = unname(LR[2]),
                minSe = unname(Se), minSp = unname(Sp)))

}
