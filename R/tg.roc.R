#' Two-Graphs Receiving Operating Characteristics.
#' @description The function supports the determination and plot of the
#'   sensitivity and specificity against the possible thresholds and shows an
#'   intermediate range of test results that is considered as less accurate.
#'
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the classification by the reference test. The reference standard
#'   must be coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition).
#' @param test The numeric test scores under evaluation. Higher scores indicate
#'   the presence of the targeted disease. Please use negated values when lower
#'   values indicate the presence of the targeted disease.
#' @param Se.criterion Default = .95. Minimum desired value of Se.
#' @param Sp.criterion Default = .95. Minimum desired value of Sp.
#' @param model Default = 'none'. Model to use, either binormal or none
#'   (non-parametric)
#' @param plot Defaults= FALSE. Whether a plot is shown for Se and Sp against
#'   the thresholds.
#' @param position.legend Default: 'left'. Position of the legend. Most used
#'   values: "left", "right".
#' @param cex.legend Default: 1. Relative size of the legend.
#'
#' @details This function implements a non-parametric and a bi-normal model. See
#'   Landsheer(2018) for an evaluative description.
#'
#'   Warning: Although the test scores <= the lower limit and the test scores >=
#'   the upper limit are interpreted for respectively negative and positive
#'   classifications, the range of test values >= lower limit provides the
#'   desired positive accuracy (Se.citerion), while the range of test values <=
#'   upper limit provides the desired negative accuracy (Sp.citerion). This is
#'   problematic in its contrariness.
#'
#'   Please note that the definition of the intermediate interval deviates from
#'   the definition of an uncertain interval.
#'
#'   The TG-ROC (Two Graphs Receiver Operating Characteristics) plot shows the
#'   diminishing values of Se and increasing values of Sp against the possible
#'   thresholds.
#' @return Thresholds for the intermediate zone. Lower threshold < Test scores <
#'   Upper threshold is the intermediate range. The range of test values >=
#'   lower limit provides the desired positive accuracy (Se.citerion), while the
#'   range of test values <= upper limit provides the desired negative accuracy
#'   (Sp.citerion).
#' @export
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

# Se.criterion = .9; Sp.criterion = .9; model = c('none', 'binormal'); plot = FALSE
TG.ROC <- function(ref,
                   test,
                   Se.criterion = .9,
                   Sp.criterion = .9,
                   model = c('none', 'binormal'),
                   plot = FALSE,
                   position.legend = 'left',
                   cex.legend = 1) {

  mod = match.arg(model) # mod='none'

  df = check.data(ref, test, mod) # nrow(df)

  simple_roc2 <- function(ref, test){
    (tab = table(test, ref)) # head(tab)
    (TPR=rev(cumsum(rev(tab[,2])/sum(tab[,2]))))
    (FPR=rev(cumsum(rev(tab[,1])/sum(tab[,1]))))
    thresholds=as.numeric(rownames(tab))
    data.frame(cbind(thresholds=thresholds, ref0 = tab[,1], ref1=tab[,2], TPR, FPR))
  }

  if (mod == 'none') {
    df2 = simple_roc2(df$ref, df$test) # head(df2); nrow(df2)
    w = which(df2$TPR >= Se.criterion)
    wL = w[length(w)]
    wU = which(df2$FPR <= (1 - Sp.criterion))[1]
    thresholds = c(L = df2$thresholds[wL], U = df2$thresholds[wU])
    if (plot) {
      plot(
        df2$thresholds,
        df2$TPR,
        type = 'l',
        col = 'red',
        xlab = 'thresholds',
        ylab = 'SeSp'
      ) # Se
      # length(df2$test); length(df2$TPR)
      lines(df2$thresholds, 1 - df2$FPR, type = 'l', col = 'blue') # Sp
      abline(v = thresholds)
      abline(h = c(Se.criterion,  Sp.criterion))
      legend(position.legend, legend = c('Se', 'Sp'), lty=c(1,1), col=c('red', 'blue'),
             cex=cex.legend)
    }
  }
  if (mod == 'binormal') {
    mu0 = mean(df$test[df$ref == 0])
    sd0 = sd(df$test[df$ref == 0])
    mu1 = mean(df$test[df$ref == 1])
    sd1 = sd(df$test[df$ref == 1])
    thresholds = c(qnorm(Se.criterion, mu1, sd1, lower.tail = F),
                   # lower threshold
                   qnorm(Sp.criterion, mu0, sd0, lower.tail = T)) # upper threshold)
    df$test = sort(df$test)
    if (plot) {
      plot(
        df$test,
        pnorm(df$test, mu0, sd0),
        type = 'l',
        col = 'blue',
        xlab = 'thresholds',
        ylab = 'SeSp'
      ) # Sp
      lines(df$test,
            1 - pnorm(df$test, mu1, sd1),
            type = 'l',
            col = 'red') # Se
      abline(v = thresholds)
      abline(h = c(Se.criterion,  Sp.criterion))
      legend(position.legend, legend = c('Se', 'Sp'), lty=c(1,1), col=c('red', 'blue'),
             cex=cex.legend)
    }
  }
  if (thresholds[1] > thresholds[2]) thresholds = c(NA, NA)
  return(thresholds)
}



