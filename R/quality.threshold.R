#' Function for describing the qualities of one or two decision thresholds
#'
#' @name quality.threshold
#'
#' @description This function can be used for both dichotomization methods (single
#'   threshold or cut-point)  and for trichotomization methods (two thresholds or
#'   cut-points). In the case of the Uncertain Interval trichotomization
#'   method, it provides descriptive statistics for the test scores outside the
#'   Uncertain Interval. For the TG-ROC trichotomization method it provides the
#'   descriptive statistics for TG-ROC's Valid Ranges.
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the classification by the reference test. The reference standard
#'   must be coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition).
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
#' @return{ A list of} \describe{ 
#' \item{$table}{The confusion table of {class x ref}, where class is the
#' classification based on the test, when applying the threshold(s). The
#' reference standard (ref) has categories 0 and 1, while the classification
#' based on the test scores (class) has categories 0 and 1 in the case of
#' applying a single threshold (dichotomization), and the categories 0,
#' 'Uncertain / NC' (NC: not classifiable) and 1 in the case of
#' trichotomization. In the case of the Uncertain Interval trichotomization
#' method, the row 'Uncertain / NC' shows the count of test scores within the
#' Uncertain Interval. When applying the trichotomization method TG-ROC, the row
#' 'Uncertain / NC' shows the count of the test scores within the Intermediate
#' Range. Table cell {0, 0} shows the True Negatives (TN), cell {0, 1} shows the
#' False Negatives (FN), cell {1, 0} shows the False Positives (FP), and cell
#' {1, 1} shows the True Positives (TP).}
#' \item{$cut}{The values of the threshold(s).}
#' \item{$indices}{A named vector, with prefix MCI when only the test scores in
#' the more certain intervals are considered. The following statistics are
#' calculated for the test-scores with classifications 0 or 1:
#' \itemize{ 
#' \item{Proportion.True: }{Proportion of true patients of all patients who
#' receive an positive or negative classification: (TP+FN)/(TN+FP+FN+TP). Equal
#' to the sample prevalence in the case of dichotomization when all patients
#' receive a positive or negative classification.}
#' \item{CCR: } {Correct Classification Rate or accuracy of the positive and negative classifications: (TP+TN)/(TN+FP+FN+TP).}
#' \item{balance : }{balance between correct and incorrect classified: (TP+TN)/(FP+FN)}
#' \item{Sp: } {specificity of the positive and negative classifications: TN/(TN+FN).}
#' \item{Se: } {sensitivity of the positive and negative classifications: TP/(TP+FN).}
#' \item{NPV: }{Negative Predictive Value of the negative class: TN/(TN+FN).}
#' \item{PPV: }{Positive Predictive Value of the positive class: TP/(TN+FN).} 
#' \item{SNPV: }{standardized negative predictive value of the negative class.}
#' \item{SPPV: }{standardized positive predictive value of the positive class.} 
#' \item{LR-: }{Negative Likelihood Ratio P(-|D+))/(P(-|D-)) The probability of a person with the
#' condition receiving a negative classification / probability of a person without the
#' condition receiving a negative classification.} 
#' \item{LR+: }{Positive Likelihood Ratio (P(+|D+))/(P(+|D-)) The probability of a person with the
#' condition receiving a positive classification / probability of a person without the
#' condition receiving a positive classification.} 
#' \item{C: }{Concordance, C-Statistic or AUC. The probability that a
#' random chosen patient with the condition is correctly ranked higher than a
#' randomly chosen patient without the condition. Equal to AUC, with for the
#' more certain interval a higher outcome than the overall concordance.}
#' } } }
#'
#' @details The Uncertain Interval is generally defined as an interval around
#'   the intersection, where the densities of the two distributions of patients
#'   with and without the targeted impairment are about equal. The various
#'   ui-functions for the estimation of the uncertain interval use a sensitivity
#'   and specificity below a desired value (default .55). Please refer to the
#'   specific function descriptions how the middle section is defined.
#'
#'   The uncertain area is defined as the scores >= threshold and <=
#'   threshold.upper. When a single threshold is supplied and no uncertain area
#'   is defined, positive classifications (1) are considered for test scores >=
#'   threshold. 
#'   
#'   Please note that the indices are calculated for those who
#'   receive a decision for or against the targeted disease: the test data in
#'   the uncertain interval are ignored. When the lower test scores indicate the
#'   targeted condition, please use the negated test values (-test).
#'
#'   The unstandardized predictive values (negative and positive; NPV and PPV)
#'   present the comparison of the observed frequencies of the two observed
#'   samples, for respectively the negative (0) and positive class (1). 
#'
#'   The standardized predictive values (SNPV and SPPV) present the comparison
#'   of the densities (or relative frequencies) of the two distributions, for
#'   the evaluated range of test scores. These predictive values are called
#'   standardized, because the two samples are compared as two independently
#'   drawn samples, not considering prevalence. It offers the estimated
#'   probability that the person (given the classification) comes from the
#'   population with (1) or from the population without (0) the target disease.
#'
#'   SNPV and SPPV provide the estimated relative probabilities that a patient
#'   is selected from the population of patients without the targeted condition
#'   or from the population of patients with the targeted condition, given that
#'   the patients test score is in the evaluated range of test scores. Of
#'   course, these estimates are better when the sample sizes are larger. 
#'   
#'   N.B. 1 When negative and predictive values would be calculated for the same
#'   range of test scores, NPV = 1 - PPV, SNPV = 1 - SPPV and PPV = 1- NPV, SPPV
#'   = 1 - SNPV.   
#'   N.B. 2 SNPV and SPPV are as independent of the prevalence as specificity 
#'   and sensitivity, as well as negative and positive probability ratios.
#' @seealso \code{\link{UncertainInterval}} for an explanatory glossary of the 
#' different statistics used within this package.
#' @examples
#' # A simple test
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1))
#' ua = ui.nonpar(ref, test)
#' quality.threshold(ref, test, threshold=ua[1], threshold.upper=ua[2])
#' # single threshold
#' quality.threshold(ref, test, threshold=ua[1])

#' @export
# threshold=-25; threshold.upper=NULL; model='kernel'; 
quality.threshold <- function(ref, test, threshold, threshold.upper=NULL,
                              model = c('kernel', 'binormal', 'ordinal')){

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
  correct.classification.rate=(TP+TN)/(TP+FP+FN+TN)
  balance.correct.incorrect=(TP+TN)/(FP+FN)

  o = outer(test[certain.sel & ref==1], test[certain.sel & ref==0], "-")
  cstat=mean((o>0) + .5*(o==0))

  inames = c('Proportion.True', 'CCR', 'balance', 'Sp', 'Se', 'NPV', 
             'PPV', 'SNPV', 'SPPV', 'LR-', 'LR+', 'C')
  if (ia) {
    lowername = '0 (test < threshold.lower)'
    uppername = '1 (test > threshold.upper)'
    t = matrix(data=c(TN, FN, ND0, ND1, FP, TP), ncol=2, byrow=T,
               dimnames=list(class=c(lowername, 'Uncertain / NC', uppername), 
                             ref=c('0', '1')))
    cut=c(threshold.lower=threshold, threshold.upper=threshold.upper)
    pt = addmargins(prop.table(t, margin = 2)	)
    likelihood.ratio.neg = pt[1,2]/pt[1,1]
    likelihood.ratio.pos = pt[3,2]/pt[3,1]
    SNPV = pt[1,1]/pt[1,3]
    SPPV = pt[3,2]/pt[3,3]
  }else{
    lowername = '0 (test < threshold)'
    uppername = '1 (test >= threshold)'
    t = matrix(data=c(TN, FN, FP, TP), ncol=2, byrow=T,
               dimnames=list(y.hat=c(lowername, uppername), ref=c('0', '1')))
    cut=c(threshold=threshold)
    uncertainty=NA
    pt = addmargins(prop.table(t[1:2 ,1:2], margin = 2)	)
    likelihood.ratio.neg = pt[1,2]/pt[1,1]
    likelihood.ratio.pos = pt[2,2]/pt[2,1]
    SNPV = pt[1,1]/pt[1,3]
    SPPV = pt[2,1]/pt[2,3]
  }
  
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
            concordance=cstat)
  
  if (ia){
    names(indices) = paste('MCI.', inames, sep = '')
  }else
    names(indices) = inames
  return(list(table=addmargins(t), cut=cut,
              indices=indices))
}
