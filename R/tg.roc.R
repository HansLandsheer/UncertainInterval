#' Two-Graphs Receiving Operating Characteristics.
#' @description The function supports the determination and plot of the
#'   sensitivity and specificity against the possible thresholds and shows an
#'   intermediate range of test results that is considered as less accurate.
#'
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the classification by the reference test. The reference standard
#'   must be coded either as 0 (absence of the condition; controls) or 1
#'   (presence of the condition; deviation from controls).
#' @param test The numeric test scores under evaluation. When
#'   \code{mean(test[ref == 0]) > mean(test[ref == 1])} it is assumed that
#'   higher test scores indicate presence of the condition, otherwise that lower
#'   test scores indicate presence of the condition.
#' @param Se.criterion Default = .95. Minimum desired value of Se.
#' @param Sp.criterion Default = .95. Minimum desired value of Sp.
#' @param model Default = 'none'. Model to use, either bi-normal or none
#'   (non-parametric)
#' @param plot Defaults= FALSE. Whether a plot is shown for Se and Sp against
#'   the thresholds.
#' @param position.legend Default: 'left'. Position of the legend. Most used
#'   values: "left", "right".
#' @param cex.legend Default: 1. Relative size of the legend.
#'
#' @details This function implements a non-parametric and a bi-normal model. See
#'   Landsheer(2018) for an evaluative description. When model='none' and the
#'   data have a limited number of values, the upper and lower threshold show
#'   the first values which comply with the criteria.
#'
#'   Warning 1: Whn using test scores where higher test scores indicate presence
#'   of a disease, the whole range of test scores starting at the lowest test
#'   score have perfect Sensitivity (1.00), at the cost of a maximal number of
#'   false positives; the sensitivity is 0.00. When moving to higher test
#'   scores, the value of sensitivity decreases and the value of specificity
#'   increases. Therefore, the lowest test scores are best used for negative
#'   classifications, but these are precisely the test scores with the highest
#'   sensitivity. In TG-ROC the test scores <= the lower limit are interpreted
#'   for negative classifications. However, the whole range of test values >=
#'   lower limit provides the minimal desired positive accuracy (Se.criterion)
#'   (at the cost of a large number of false positives). Similarly, the test
#'   scores >= the upper limit are interpreted for positive classifications,
#'   while the whole range of test values <= upper limit provides the minimal
#'   desired negative accuracy (Sp.criterion) (at the cost of a large number of
#'   false negatives). Of course, this is also true for tests where the lowest
#'   scores indicate the presence of the disease, but only reversed.
#'
#'   Warning 2: The Intermediate range can cover a relatively small part of the
#'   area of overlap between the two distributions. In that case test scores
#'   with relative low number of false classifications are considered as
#'   intermediate.
#'
#'   Please note that the definition of the intermediate interval deviates
#'   substantially from the definition of an uncertain interval.
#'
#'   The TG-ROC (Two Graphs Receiver Operating Characteristics) plot shows the
#'   diminishing values of Se and increasing values of Sp against the possible
#'   thresholds.
#' @return Thresholds for the intermediate zone. Lower threshold < Test scores <
#'   Upper threshold is the intermediate range. The range of test values >=
#'   lower limit provides the desired positive accuracy (Se.criterion), while
#'   the range of test values <= upper limit provides the desired negative
#'   accuracy (Sp.criterion).
#' @references Greiner, M. (1995). Two-graph receiver operating characteristic
#'   (TG-ROC): A Microsoft-EXCEL template for the selection of cut-off values in
#'   diagnostic tests. Journal of Immunological Methods, 185(1), 145-146.
#'
#'   Greiner, M. (1996). Two-graph receiver operating characteristic (TG-ROC):
#'   Update version supports optimisation of cut-off values that minimise
#'   overall misclassification costs. Journal of Immunological Methods, 191(1),
#'   93-94.
#'
#'   Landsheer, J. A. (2018). The Clinical Relevance of Methods for Handling
#'   Inconclusive Medical Test Results: Quantification of Uncertainty in Medical
#'   Decision-Making and Screening. Diagnostics, 8(2), 32.
#'   https://doi.org/10.3390/diagnostics8020032
#'
#' @examples
#' ref = c(rep(0,100), rep(1,100))
#' test = c(rnorm(100, 0, 1), rnorm(100, 1, 1))
#' TG.ROC(ref, test, model='binormal', plot=TRUE)
#' TG.ROC(ref, test, model='none', plot=TRUE)
#' @export
#' @importFrom stats qnorm
# Se.criterion = .9; Sp.criterion = .9; model = 'binormal'
# plot = T; position.legend = 'left'; cex.legend = 1; test=-test
TG.ROC <- function(ref,
                   test,
                   Se.criterion = .9,
                   Sp.criterion = .9,
                   model = c('none', 'binormal'),
                   plot = FALSE,
                   position.legend = 'left',
                   cex.legend = 1) {

  model = match.arg(model) # model='none'

  df = check.data(ref, test, model) # nrow(df)

  normhigher = mean(test[ref==0]) > mean(test[ref==1])

  if (model == 'none') {
    df2 = simple_roc3(df$test[df$ref==0], df$test[df$ref==1]) # head(df2)
    if (normhigher){
      w = which(df2$fpr <= 1-Sp.criterion) # df2$fpr[w]
      wL = w[length(w)]
      wU = which(df2$tpr >= Se.criterion)[1]
    } else {
      w = which(df2$tpr >= Se.criterion) # df2$tpr[w]
      wL = w[length(w)]
      wU = which(df2$fpr <= (1 - Sp.criterion))[1]
    }
    thresholds = c(L = df2$testscores[wL], U = df2$testscores[wU])
    if (plot) {
      plot(
        df2$testscores,
        df2$tpr,
        type = 'l',
        col = 'red',
        xlab = 'thresholds',
        ylab = 'SeSp'
      ) # Se
      # length(df2$test); length(df2$TPR)
      lines(df2$testscores, 1 - df2$fpr, type = 'l', col = 'blue') # Sp
      abline(v = thresholds)
      abline(h = c(Se.criterion,  Sp.criterion))
      legend(position.legend, legend = c('Se', 'Sp'), lty=c(1,1), col=c('red', 'blue'),
             cex=cex.legend)
    }
  }
  if (model == 'binormal') {
    mu0 = mean(df$test[df$ref == 0])
    sd0 = sd(df$test[df$ref == 0])
    mu1 = mean(df$test[df$ref == 1])
    sd1 = sd(df$test[df$ref == 1])
    thresholds = c(L=qnorm(Se.criterion, mu1, sd1, lower.tail = normhigher),
                     U=qnorm(Sp.criterion, mu0, sd0, lower.tail = !normhigher))
    if (normhigher) {
      tmp=thresholds[1]; thresholds[1]=thresholds[2]; thresholds[2]=tmp
    }
    df$test = sort(df$test)
    if (plot) {
      plot(
        df$test,
        pnorm(df$test, mu0, sd0, lower.tail = !normhigher),
        type = 'l',
        col = 'blue',
        xlab = 'thresholds',
        ylab = 'SeSp'
      ) # Sp
      lines(df$test,
            pnorm(df$test, mu1, sd1, lower.tail = normhigher),
            type = 'l',
            col = 'red') # Se
      abline(v = thresholds)
      abline(h = c(Se.criterion,  Sp.criterion))
      legend(position.legend, legend = c('Se', 'Sp'), lty=c(1,1), col=c('red', 'blue'),
             cex=cex.legend)
    }
  }
  return(thresholds)
}



