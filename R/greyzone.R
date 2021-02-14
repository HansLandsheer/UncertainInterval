#' Function for the determination of a grey zone for ordinal diagnostic and
#' screening tests
#'
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the reference or gold standard. The reference standard must be
#'   coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition).
#' @param test The numeric test scores under evaluation. When
#'   \code{mean(test[ref == 0]) > mean(test[ref == 1])} it is assumed that
#'   higher test scores indicate presence of the condition, otherwise that lower
#'   test scores indicate presence of the condition.
#' @param pretest.prob The pre-test probability to be used. When
#'   NULL, the prevalence found in the sample is used.
#' @param criterion.values The minimum desired values for respectively the
#'   negative and positive post-test probability.
#' @param return.all Default = FALSE. When TRUE the full table of all results
#'   are returned.
#' @param precision Default = 3. Precision used for comparison of the criterion
#'   values and the post-test probabilities.
#'
#' @details This function is proposed by Coste et al. (2003). The current
#'   implementation only handles ordinal test values. This functions uses all
#'   possible test scores as dichotomous thresholds to calculate Se, Sp,
#'   positive and negative likelihood ratios and post-test probabilities. The
#'   likelihood ratios are calculated for the accumulated densities of the test
#'   scores and indicate the levels of seriousness of the disease for all
#'   possible dichotomous thresholds. It uses therefore a cumulative
#'   interpretation of the Likelihood Ratios and posttest probabilities. If a
#'   test has test scores 1 to 5 (with 5 indicating the largest probability of
#'   the disease), Se, positive LR and positive posttest probabilities of the
#'   greyzone function uses dichotomous thresholds that concern  test results >=
#'   1, >= 2, >= 3, >= 4 and >= 5, while Sp, negative LR and negative posttest
#'   probabilities concern test results <= 1, <= 2, <= 3, <= 4 and <= 5. Please
#'   note that in these examples values <= 1 respectively <= 5 concern all
#'   possible test values and have by definition a dichotomous Sensitivity of 1.
#'
#'   Please note that the definition of a grey zone deviates from the definition
#'   of an uncertain interval.
#'
#'   The decision criteria are the required values of post-test probabilities.
#'   This has changed in version 0.7. In earlier versions the criteria was a
#'   value closest to the criterion, which could produce invalid results. These
#'   post-test probabilities of accumulated test scores may require a value over
#'   0.99 or even 0.999 (or under 0.01 or 0.001) to confirm or exclude the
#'   presence of a target disease. Only tests of the highest quality can reach
#'   such criteria. The default criterion values are .05 and .95 for
#'   respectively a negative and positive classification, which may be
#'   sufficient for use by clinicians or Public Health professionals for a first
#'   classification whether a target disease may be present or not (Coste et
#'   al., 2003).
#'
#'   As such the cumulative likelihood ratios differ from the Interval
#'   Likelihood Ratios (see \code{\link{RPV}}), as proposed by Sonis (1999).
#'   These likelihood ratios are calculated for each given interval of test
#'   scores separately and uses their densities. In contrast to the greyzone
#'   method, Interval Likelihood ratios and interval posttest probabilities
#'   concern the separate intervals, that is in this example, the separate score
#'   1 to 5. Interval likelihood ratios assign a specific value to each level of
#'   abnormality, and this value is used to calculate the posttest probabilities
#'   of disease for each given level of a test (Sonis, 1999). These post-test
#'   probabilities differ strongly from the cumulative post-test probabilities
#'   and criterion values can be much lower, especially when diseases are life
#'   threatening and low-cost treatments are available. See Sonis (1999) for
#'   further discussion of the interval interpretation.
#'
#' @return The function returns the lower and upper value of the range of test
#'   scores that are considered 'grey' or inconclusive. Only smaller or larger
#'   values are considered for a decision. When return.all = TRUE the full table
#'   of the results is returned.
#' @seealso \code{\link{RPV}}
#' @export
#' @references { Coste, J., Jourdain, P., & Pouchot, J. (2006). A gray zone
#'   assigned to inconclusive results of quantitative diagnostic tests:
#'   application to the use of brain natriuretic peptide for diagnosis of heart
#'   failure in acute dyspneic patients. Clinical Chemistry, 52(12), 2229-2235.
#'
#'   Coste, J., & Pouchot, J. (2003). A grey zone for quantitative diagnostic
#'   and screening tests. International Journal of Epidemiology, 32(2), 304-313.
#'
#'   Sonis, J. (1999). How to use and interpret interval likelihood ratios.
#'   Family Medicine, 31, 432-437. }
#'
#' @examples
#'  ref=c(rep(0, 250), rep(1, 250))
#'  test = c(rep(1:5, c(90,75,50,34,1)), c(rep(1:5, c(10,25,50,65,100))))
#'  addmargins(table(ref, test))
#'  greyzone(ref, test, ret=TRUE, criterion.values=c(.1, .9))
#'
#'  test = c(rep(14:31, c(0,0,0,0,0,0,3,3,5,7,10,20,30,40,50,24,10,10)),
#'           rep(14:31, c(1,0,0,0,0,0,1,4,4,9, 6,13, 8, 6, 5, 4, 0, 0)))
#'  ref = c(rep(0, 212), rep(1, 61))
#'  barplotMD(ref, test)
#'  addmargins(table(ref, test))
#'  greyzone(ref, test, ret=TRUE, crit=c(.1,.9))

# ref=irondef$ref; test=irondef$test; pretest.prob=.1; criterion.values=c(.7, .1); return.all=TRUE
# pretest.prob=NULL; criterion.values=c(.10, .90); return.all=T; precision=3
greyzone <- function(ref, test, pretest.prob=NULL, 
                     criterion.values=c(.05, .95), return.all=F, precision=3){

  df=check.data(ref, test, model='ordinal')
  
  controlslower = mean(df$test[ref==0]) > mean(df$test[ref==1])

  if (is.null(pretest.prob)) pretest.prob = sum(df$ref==1)/length(df$ref)
  
  rdf = simple_roc3(test[ref==0], test[ref==1])

  totneg = length(df$test[ref==0])
  totpos = length(df$test[ref==1])
  rdf$TN=totneg-rdf$FP
  rdf$FN=totpos-rdf$TP

  rdf$tpr=rdf$TP/totpos          # Sensitivity (fraction true positives)
  rdf$fpr=rdf$FP/totneg  # 1-rdf$FP/totneg         # 1 - specificity (false positives)

  rdf$preodds=pretest.prob/(1-pretest.prob); # pretest odds
  rdf$sp=1-rdf$fpr;
  rdf$se=rdf$tpr;
  rdf$neg.lr=(1-rdf$se)/rdf$sp; # negative likelihood ratio
  rdf$pos.lr=rdf$se/(1-rdf$sp); # positive likelihood ratio
  # By Bayes' rule : posttest(posterior)odds=likelihoodratio * pretest(prior)odds
  rdf$negpostodds=rdf$neg.lr * rdf$preodds ; # negative post odds
  rdf$pospostodds=rdf$pos.lr * rdf$preodds ; # positive post odds
  rdf$negpostprob = ifelse (is.infinite(rdf$negpostodds),1, 
                            rdf$negpostodds / (rdf$negpostodds+1)) # negative post probability
  rdf$pospostprob = ifelse (is.infinite(rdf$pospostodds),1, 
      rdf$pospostodds / (rdf$pospostodds+1)) # positive post probability
  # alternatively
  # rdf$pospostprob = (pretest.prob*rdf$tpr)/(pretest.prob*rdf$tpr+(1-pretest.prob)*rdf$fpr)
  # rdf$negpostprob = 1-((1-pretest.prob)*(1-rdf$fpr))/(pretest.prob*(1-rdf$tpr)+(1-pretest.prob)*(1-rdf$fpr))

  # determine points of intercept with criteria
  # Find points where pospostprob is below negpost.crit.
  # > decreasing negpostprob; < increasing negpostprob
  below <- round(rdf$negpostprob, precision) <= criterion.values[1]
  # Points always intersect when above=TRUE, then FALSE or reverse
  intersect.points.neg <- which(diff(below) != 0)
  if (!controlslower) intersect.points.neg = intersect.points.neg+1
  # Find points where pospostprob is above pospost.crit.
  # > decreasing pospostprob; < increasing negpostprob
  above <- round(rdf$pospostprob, precision) >= criterion.values[2]
  # Points always intersect when above=TRUE, then FALSE or reverse
  intersect.points.pos <- which(diff(above) != 0)
  if (controlslower) intersect.points.pos = intersect.points.pos+1

  ltpos=utpos=NA
  if(length(intersect.points.neg)>0) ltpos = intersect.points.neg # always use last value
  if(length(intersect.points.pos)>0) utpos = intersect.points.pos # always use first value

  if(length(ltpos) > 1 | length(utpos) > 1) warning('Multiple thresholds found.')

  # controlslower; rdf
  if (return.all) {
    return(list(table = rdf, 
                thresholds.grey = c(neg = rdf$testscores[ltpos], 
                                    pos = rdf$testscores[utpos])) )
  } else return(c(neg = rdf$testscores[ltpos], 
                  pos = rdf$testscores[utpos])) 
}


