#' Function for the description of the qualities of the Uncertain Interval.
#' @name quality.threshold.uncertain
#' @description This function can be used only for trichotomization (double
#'   thresholds or cut-points) methods. In the case of the Uncertain Interval
#'   trichotomization method, it provides descriptive statistics for the test
#'   scores within the Uncertain Interval. For the TG-ROC trichotomization
#'   method it provides the descriptive statistics for TG-ROC's Intermediate
#'   Range.
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the classification by the reference test. The reference standard
#'   must be coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition)
#' @param test The index test or test under evaluation. A column in a dataset or
#'   vector indicating the test results in a continuous scale.
#' @param threshold The lower decision threshold of a trichotomization method.
#' @param threshold.upper The upper decision threshold of a trichotomization
#'   method. Required.
#' @param intersection (default = NULL). When NULL, the intersection is
#'   calculated with \code{get.intersection}, which uses the kernel density
#'   method to obtain the intersection. When another value is assigned to this
#'   parameter, this value is used instead.
#' @param model (default = 'kernel'). The model used defines the intersection.
#'   Default the kernel densities are used with adjust = 1, for ordinal models
#'   adjust = 2 is used. For binormal models the binormal estimate of the
#'   intersection is used. The model defines the intersection, which defines the
#'   output of this function.
#' @return{ A list of} \describe{ \item{intersection}{The value used as estimate
#' of the intersection (that is, the optimal threshold). NOTE: The
#' trichotomization method TG-ROC has no defined position for its Intermediate
#' Range, but usage of the point where Sensitivity=Specificity seems a
#' reasonable choice.} \item{table}{The confusion table of {UI.class x ref} for
#' the Uncertain Interval where the scores are expected to be inconclusive.
#' UI.class is the classification of the UI scores divided by the intersection,
#' 0 (UI scores < intersection and 1 (UI scores >= intersection. This shows
#' therefore the results when compared with applying the intersection is used
#' (that is, the optimal dichotomous threshold). Both the reference standard
#' (ref) and the classification based on the test scores (UI.class) have
#' categories 0 and 1. Table cell {0, 0} shows the True Negatives (TN), cell {0,
#' 1} shows the False Negatives (FN), cell {1, 0} shows the False Positives
#' (FP), and cell {1, 1} shows the True Positives (TP).} \item{cut}{The values
#' of the thresholds.} \item{X2}{Table with the outcomes of three Chi-square
#' tests of the confusion table:} \itemize{ \item{TN.FP: }{Chi-square test of
#' the comparison of TN versus FP.} \item{FN.TP: }{Chi-square test of the
#' comparison of FN versus TP.} \item{overall: }{Chi-square test of all four
#' cells of the table.} } \item{t.test}{Table with t-test results for the
#' comparison of the means. Within the Uncertain Interval, the test scores are
#' compared of individuals without the targeted condition (ref = 0) and
#' individuals with the targeted condition (ref = 1).} \item{indices}{A named
#' vector, with the following statistics for the test-scores within the
#' Uncertain Interval:} \itemize{ \item{prevalence: }{Classifiable patients with
#' the targeted condition / Total sample = (TP+FN)/(TN+FP+FN+TP)}
#' \item{correct.classification.rate (or Accuracy): }{(TP+TN)/(TN+FP+FN+TP)}
#' \item{balance.correct.incorrect: }{(TP+TN)/(FP+FN)} \item{specificity:
#' }{TN/(TN+FN)} \item{sensitivity: }{TP/(TP+FN)}
#' \item{negative.predictive.value: }{TN/(TN+FN)}
#' \item{positive.predictive.value: }{TP/(TN+FN)} \item{SNPV: standardized
#' negative predictive value = specificity / (1- sensitivity + specificity)}
#' \item{SPPV: standardized positive predictive value = sensitivity /
#' (sensitivity +  1 - specificity)} \item{neg.likelihood.ratio:
#' }{(1-sensitivity)/specificity} \item{pos.likelihood.ratio:
#' }{sensitivity/(1-specificity)} \item{concordance: }{The probability that a
#' random chosen patient with the condition is correctly ranked higher than a
#' randomly chosen patient without the condition. Equal to AUC, with for the
#' uncertain interval an expected outcome < .60. (Not equal to a partial AUC.)}
#' } }
#' @details The Uncertain Interval is generally defined as an interval below and
#'   above the intersection, where the densities of the two distributions of
#'   patients with and without the targeted impairment are about equal. The
#'   various functions for the estimation of the uncertain interval use a
#'   sensitivity and specificity below a desired value (default .55). This
#'   function uses the intersection (the optimal dichotomous threshold) to
#'   divide the uncertain interval and provides in this way the indices for the
#'   uncertain interval when the optimal threshold would have been applied.
#'
#'   As a result, it may be expected that Chi-square tests are not significant,
#'   provided that the count of individuals within the Uncertain Interval is not
#'   too large. Most often, the t-test is also non-significant, but as the power
#'   of the t-test is considerably larger than the power of the Chi-square test,
#'   this is less often the case. It is recommended to look at the difference of
#'   the means of the two sub-samples and to visually inspect the
#'   inter-mixedness of the densities of the test scores.
#'
#'   The patients that have test scores within the Uncertain Interval are prone
#'   to be incorrectly classified on the basis of their test result. The results
#'   within the Uncertain Interval differ only slightly for patients with and
#'   without the targeted condition. Patients with slightly lower or higher test
#'   scores too often have the opposite status. They receive the classification
#'   result 'Uncertain'; it is better to apply additional tests or to await
#'   further developments.
#'
#'   When applying the method to the results of a logistic regression, one
#'   should be aware of possible problems concerning the determination of the
#'   intersection. Somewhere in the middle, logistic predictions can have a
#'   range where the distributions have similar densities or have multiple
#'   intersections near to each other. Often, this problem can be approached
#'   effectively by using the linear predictions instead of the logistic
#'   predictions. The linear predictions offer often a far more clear point of
#'   intersection. The solution can then be applied to the prediction values
#'   using the inverse logit of the intersection and the two cut-points. The
#'   logistic predictions and the linear predictions have the same rank
#'   ordering.
#'
#' @export
#'
#' @examples
#' # A simple test model
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,sd=1))
#' ua = ui.nonpar(ref, test)
#' quality.threshold.uncertain(ref, test, ua[1], ua[2])
quality.threshold.uncertain <- function(ref, test,
                                       threshold, threshold.upper, intersection=NULL,
                                       model = c('kernel', 'binormal', 'ordinal')){
  # threshold = threshold.upper = 19; intersection=NULL
  model <- match.arg(model)
  df=check.data(ref, test, model=model)
  ref=df$ref
  test=df$test
  stopifnot(!is.null(threshold.upper))

  # chisq test equal probability of two frequencies
  chisq1 <- function(x){
    E <- sum(x)/(length(x)) # equal frequencies expected
    if (all(E >= 5)) {
      X2=sum((x - E)^ 2 / E)
      return(c(n=x[1], n=x[2], sum=sum(x), X2=round(X2,3), df= 1,
               p=round(pchisq(X2, df=1, lower.tail=F),3))) } else
                 return(c(n1=x[1], n2=x[2], sum=sum(x), X2=NA, df= 1, p=NA))
  }

  # chisq test 2x2 table
  # x=matrix(c(TN,FP,FN,TP), byrow=T, nrow=2)
  chisq2 <- function(x){
    n = sum(x)
    rs <- rowSums(x)
    cs <- colSums(x)
    if(n==0) return(c(n0=cs[1], n1=cs[2], sum=sum(cs), X2=NA, df= 1, p=NA))
    E <- outer(rs, cs, "*") / n
    if (all(E >= 5)) {
      YATES <- min(0.5, abs(x - E))
      X2=sum((abs(x - E) - YATES) ^ 2 / E)
      return(c(n=cs[1], n=cs[2], sum=sum(cs), X2=round(X2,3), df= 1,
               p=round(pchisq(X2, df=1, lower.tail=F),3)))
    } else return(c(n0=cs[1], n1=cs[2], sum=sum(cs), X2=NA, df= 1, p=NA))
  }

  if (threshold.upper < threshold) { temp=threshold; threshold=threshold.upper; threshold.upper=temp}

  threshold=unname(unlist(threshold)) # threshold=ua[1]
  threshold.upper=unname(unlist(threshold.upper)) # threshold.upper=NULL

    if (is.null(intersection)) {
      intersection = tail(get.intersection(ref, test, model=model),1)
      if (length(intersection) > 1) {
        intersection=tail(intersection, n=1)
        warning('More than one point of intersection. Highest used.')
      }
    }
    stopifnot(threshold <= intersection | threshold.upper >= intersection)

  thresholdnew = intersection # center point inside the uncertain interval

  certain.sel = (test < threshold) | (test > threshold.upper)
  # uncertain.obs = sum(!certain.sel)
  # if (!raw) {
  #   uncertainty.all = mean((ref - test) ^ 2)
  #   uncertainty.certain.interval = mean((ref[certain.sel] - test[certain.sel]) ^2)
  #   uncertainty.uncertain.interval = mean((ref[!certain.sel] - test[!certain.sel]) ^2)
  # }

  test.uc = test[!certain.sel] # sum(!certain.sel); length(test.uc)
  ref.uc = ref[!certain.sel] # length(ref.uc)

  y.hat=rep(0, length(ref.uc))
  # only one relevant intersection assumed!

  y.hat[test.uc >= thresholdnew]=1
  #X2=NA
  #X2=try(chisq.test(table(y.hat,ref.uc), simulate.p.value = FALSE), silent=T)
  #(chisq.test(table(y.hat,ref.uc)[,1], simulate.p.value = FALSE))
  #(chisq.test(table(y.hat,ref.uc)[,2], simulate.p.value = FALSE))

  TP = sum(y.hat==1 & ref.uc==1)
  FP = sum(y.hat==1 & ref.uc==0)
  TN = sum(y.hat==0 & ref.uc==0)
  FN = sum(y.hat==0 & ref.uc==1)

  TN.FP=chisq1(c(TN, FP))
  FN.TP=chisq1(c(FN, TP))
  overall=chisq2(matrix(c(TN,FP,FN,TP), byrow=T, nrow=2))

  # test.uc = p0; ref.uc = d0
  t0 = test.uc[ref.uc==0]; n0 = as.numeric(length(t0));
  t1 = test.uc[ref.uc==1]; n1 = as.numeric(length(t1));

  if (n0 > 1 & n1 > 1 & threshold!=threshold.upper) {
    res = t.test(t0, t1)
    t = c(mean.0 = unname(res$estimate[1]),mean.1 = unname(res$estimate[2]),
          res$statistic, res$parameter, p = res$p.value)
  } else {
    t = c(mean.0 = mean(t0),
          mean.1 = mean(t1), t = NA, df = NA, p = NA)
  }

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
  bww = sum(ref.uc==1)/sum(ref.uc==0) # 1/bww

  ta=addmargins(matrix(c(TN,FP,FN,TP),2,2))
  lowername = '0 (threshold.lower <= test < intersection)'
  uppername = '1 (intersection <= test <= threshold.upper)'
  dimnames(ta)=list(UI.class=c(lowername, uppername, 'Sum'), ref= c('0', '1', 'Sum'))

  # if(!raw){
  #   uncertainty = c(uncertain.obs=uncertain.obs, uncertainty.all=uncertainty.all,
  #                   uncertainty.certain.interval=uncertainty.certain.interval,
  #                   uncertainty.uncertain.interval=uncertainty.uncertain.interval)
  # } else{uncertainty = c(uncertain.obs=uncertain.obs)}
  # o = outer(test.uc[ref.uc == 1], test.uc[ref.uc == 0], "-")
  # cstat=mean((o>0) + .5*(o==0))

  r = rank(c(t1,t0))
  cstat = (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n0)


  return(list(intersection=thresholdnew,
              table=ta,
              cut=c(threshold.lower=threshold, threshold.upper=threshold.upper),
              X2=rbind(TN.FP, FN.TP, overall),
              t.test=t,
              # uncertainty=uncertainty,
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
                        balance.with.without = bww,
                        concordance=cstat)) )
}
