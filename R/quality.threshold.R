#' Function for the description of the qualities of one or two decision
#' thresholds or threshold.
#'
#' @name quality.threshold
#'
#' @description This function can be used for both dichotomization (single
#'   threshold or cut-point) methods and for trichotomization (two thresholds or
#'   cut-points) methods. In the case of the Uncertain Interval trichotomization
#'   method, it provides descriptive statistics for the test scores outside the
#'   Uncertain Interval. For the TG-ROC trichotomization method it provides the
#'   descriptive statistics for TG-ROC's Valid Ranges.
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the classification by the reference test. The reference standard
#'   must be coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition)
#' @param test The index test or test under evaluation. A column in a dataset or
#'   vector indicating the test results in a continuous scale.
#' @param threshold The decision threshold of a dichotomization method, or the
#'   lower decision threshold of a trichotomization method.
#' @param threshold.upper (default = NULL). The upper decision threshold of a
#'   trichotomization method. When NULL, the test scores are dichotomized and
#'   only threshold is used for the dichotomization.
#' @param model The model to use. Default = 'kernel' for continuous data. For
#'   discrete data 'ordinal' is the better choice.
#'
#' @return{ A list of} \describe{ \item{$table}{The confusion table of {class x
#' ref}, where class is the classification based on the test, when applying the
#' threshold(s). The reference standard (ref) has categories 0 and 1, while the
#' classification based on the test scores (class) has categories 0 and 1 in the
#' case of applying a single threshold (dichotomization), and the categories 0,
#' NA and 1 in the case of trichotomization. In the case of the Uncertain
#' Interval trichotomization method, the row NA shows the count of test scores
#' within the Uncertain Interval. When applying the trichotomization method
#' TG-ROC, the row NA shows the count of the test scores within the Intermediate
#' Range. Table cell {0, 0} shows the True Negatives (TN), cell {0, 1} shows the
#' False Negatives (FN), cell {1, 0} shows the False Positives (FP), and cell
#' {1, 1} shows the True Positives (TP).} \item{$cut}{The values of the
#' threshold(s).} \item{$indices}{A named vector, with the following statistics
#' for the test-scores with classifications 0 or 1: \itemize{ \item{prevalence:
#' }{Proportion of classifiable patients with the targeted condition
#' = (TP+FN)/(TN+FP+FN+TP)} \item{correct.classification.rate (or Accuracy):
#' }{(TP+TN)/(TN+FP+FN+TP)} \item{balance.correct.incorrect : }{(TP+TN)/(FP+FN)}
#' \item{specificity: }{TN/(TN+FN)} \item{sensitivity: }{TP/(TP+FN)}
#' \item{negative.predictive.value: }{TN/(TN+FN)}
#' \item{positive.predictive.value: }{TP/(TN+FN)} \item{SNPV: standardized
#' negative predictive value = specificity / (1- sensitivity + specificity)}
#' \item{SPPV: standardized positive predictive value = sensitivity /
#' (sensitivity +  1 - specificity)} \item{neg.likelihood.ratio:
#' }{(1-sensitivity)/specificity} \item{pos.likelihood.ratio:
#' }{sensitivity/(1-specificity)} \item{concordance: }{The probability that a
#' random chosen patient with the condition is correctly ranked higher than a
#' randomly chosen patient without the condition. Equal to AUC, with for the
#' more certain interval a higher outcome than the overall concordance.}
#'
#' } } }
#'
#' @details The Uncertain Interval is generally defined as an interval below and
#'   above the intersection, where the densities of the two distributions of
#'   patients with and without the targeted impairment are about equal. The
#'   various ui-functions for the estimation of the uncertain interval use a
#'   sensitivity and specificity below a desired value (default .55). Please
#'   refer to the specific function descriptions how the middle section is
#'   defined.
#'
#'   The uncertain area is defined as the scores >= threshold and <=
#'   threshold.upper. When a single threshold is supplied and no uncertain area
#'   is defined, positive classifications (1) are considered for test scores >=
#'   threshold. Please note that the indices are calculated for those who
#'   receive a decision for or against the targeted disease: the test data in
#'   the uncertain interval are ignored. When the lower test scores indicate the
#'   targeted condition, please use the negated test values (-test).
#'
#'   The unstandardized predictive values (negative and positive; NPV and PPV)
#'   present the comparison of the observed frequencies of the two observed
#'   samples, for the evaluated range of test scores. When using a single
#'   cut-point the evaluated range of test scores for the NPV are test scores <
#'   threshold and for PPV test scores >= threshold. When using
#'   trichotomization, the evaluated range is < lower limit (NPV) and > upper
#'   limit (PPV). They provide the exact observed proportions of correctly
#'   classified patients, given the range of test scores.
#'
#'   The standardized predictive values (SNPV and SPPV) present the comparison
#'   of the densities (or relative frequencies) of the two distributions, for
#'   the evaluated range of test scores. When using a single cut-point the
#'   evaluated range of test scores for the SNPV are test scores < threshold and
#'   for SPPV test scores >= threshold. When using trichotomization, the
#'   evaluated range is < lower limit (SNPV) and > upper limit (SPPV). These
#'   predictive values are called standardized, because the two samples are
#'   compared as two independent samples, as if prevalence equals .5.
#'
#'   SNPV and SPPV provide the estimated relative probabilities that a patient
#'   is selected from the population of patients without the targeted condition
#'   or from the population of patients with the targeted condition, given that
#'   the patients test score is in the evaluated range of test scores. Of
#'   course, these estimates are better when the sample sizes are larger. N.B. 1
#'   When negative and predictive values would be calculated for the same range
#'   of test scores, SNPV = 1 - SPPV and SPPV = 1 - SNPV. N.B. 2 SNPV and SPPV
#'   are as independent of prevalence as are specificity and sensitivity and as
#'   are negative and positive Likelihood Ratios.
#'
#' @examples
#' # A simple test
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1))
#' ua = ui.nonpar(ref, test)
#' quality.threshold(ref, test, threshold=ua[1], threshold.upper=ua[2])

#' @export
# threshold=-25; threshold.upper=NULL; model='ordinal'; uncertain = 'ignore';
quality.threshold <- function(ref, test, threshold, threshold.upper=NULL,
                              model = c('kernel', 'binormal', 'ordinal')){
  # .Deprecated('quality.treshold', msg = 'Deprecated. Replaced by quality.mci')

  model <- match.arg(model)
  df=check.data(ref, test, model=model)
  ref=df$ref
  test=df$test
  ia = !is.null(threshold.upper)
  # raw=T
  threshold=unname(unlist(threshold)) # threshold=ua[1]
  threshold.upper=unname(unlist(threshold.upper)) # threshold.upper=NULL
  certain.sel=rep(TRUE, length(ref))
  y.hat=rep(1, length(ref)) # init

  ND0=0
  ND1=0
  if (ia) {
    if (threshold.upper < threshold) {
      temp=threshold; threshold=threshold.upper; threshold.upper=temp
    }
    certain.sel = (test < threshold) | (test > threshold.upper)
    uncertain.obs = sum(!certain.sel)
    # if (!raw){
    #   uncertainty.all = mean((ref-test)^2)
    #   uncertainty.certain.interval = mean((ref[certain.sel]-test[certain.sel])^2)
    #   uncertainty.uncertain.interval = mean((ref[!certain.sel]-test[!certain.sel])^2)
    # }
    y.hat[!certain.sel]= NA
    ND0 = sum(is.na(y.hat) & ref==0)
    ND1 = sum(is.na(y.hat) & ref==1)
  }

  y.hat[test < threshold]=0  # change <=
  # table(y.hat,ref)

  TN = sum(y.hat==0 & ref==0, na.rm=T)
  FN = sum(y.hat==0 & ref==1, na.rm=T)
  FP = sum(y.hat==1 & ref==0, na.rm=T)
  TP = sum(y.hat==1 & ref==1, na.rm=T)

  prevalence=(TP+FN)/(TP+FP+FN+TN)
  sensitivity = TP/(TP+FN)
  specificity = TN/(FP+TN)
  positive.predictive.value=TP/(TP+FP)
  negative.predictive.value=TN/(FN+TN)
  SNPV = specificity / (specificity + 1 - sensitivity)
  SPPV = sensitivity / (sensitivity + 1 - specificity)
  correct.classification.rate=(TP+TN)/(TP+FP+FN+TN)
  balance.correct.incorrect=(TP+TN)/(FP+FN)
  likelihood.ratio.pos = sensitivity /(1-specificity)
  likelihood.ratio.neg = (1-sensitivity)/specificity

  o = outer(test[certain.sel & ref==1], test[certain.sel & ref==0], "-")
  cstat=mean((o>0) + .5*(o==0))

   # threshold = 10; threshold.upper = 13
  if (ia) {
    lowername = '0 (test < threshold.lower)'
    uppername = '1 (test > threshold.upper)'
    t = matrix(data=c(TN, FN, ND0, ND1, FP, TP), ncol=2, byrow=T,
               dimnames=list(class=c(lowername, 'NA', uppername), ref=c('0', '1')))
    cut=c(threshold.lower=threshold, threshold.upper=threshold.upper)
    # if (!raw) {uncertainty = c(uncertain.obs=uncertain.obs, uncertainty.all=uncertainty.all,
    #                            uncertainty.certain.interval=uncertainty.certain.interval,
    #                            uncertainty.uncertain.interval=uncertainty.uncertain.interval)
    # }else{
    #   uncertainty=c(uncertain.obs=uncertain.obs)
    # }
  }else{
    lowername = '0 (test < threshold)'
    uppername = '1 (test >= threshold)'
    t = matrix(data=c(TN, FN, FP, TP), ncol=2, byrow=T,
               dimnames=list(y.hat=c(lowername, uppername), ref=c('0', '1')))
    cut=c(threshold=threshold)
    uncertainty=NA
  }
  return(list(table=addmargins(t), cut=cut,
              # uncertainty = uncertainty,
              indices=c(prevalence=prevalence,
                        correct.classification.rate=correct.classification.rate,
                        balance.correct.incorrect=balance.correct.incorrect,
                        specificity =specificity,
                        sensitivity =sensitivity,
                        negative.predictive.value=negative.predictive.value,
                        positive.predictive.value=positive.predictive.value,
                        SNPV = SNPV,
                        SPPV = SPPV,
                        neg.likelihood.ratio = likelihood.ratio.neg,
                        pos.likelihood.ratio = likelihood.ratio.pos,
                        concordance=cstat)))
}
