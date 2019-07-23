#' Calculation of Predictive Values, Standardized Predictive Values, Interval
#' Likelihood Ratios and Posttest Probabilities of intervals or individual test
#' scores of discrete ordinal tests.
#'
#' @export
#' @importFrom zoo rollsum na.fill
#' @description This function can correct for the unreliability of the test. It
#'   also trichotomizes the test results, with an uncertain interval where the
#'   test scores do not allow for an adequate distinction between the two groups
#'   of patients. This function is best applied to large samples with a
#'   sufficient number of patients for each test score.
#' @param ref A vector of two values, ordering 'patients with' > 'patients
#'   without', for instance 0, 1.When using a factor, please check whether the
#'   correct order of the values is used.
#' @param test A vector of ordinal measurement level with a numeric score for
#'   every individual. When using a factor, please check whether the correct
#'   order of the values is used. Further, a warning message is issued
#'   concerning the calculation of the variance of the test.
#' @param reliability (no default) The known reliability of the test, used to
#'   calculate Standard Error of Measurement (SEM). The reliability is expressed
#'   as an applicable correlation coefficient, with values between 0 and 1. A
#'   Pearson's Product Moment correlation or an Intra-Class Coefficient (ICC)
#'   will do. N.B. Setting parameter \code{roll.length} to 1 causes the
#'   reliability of the test to be ignored.
#' @param pretest.prob (default = NULL) value to be used as pre-test
#'   probability. It is used for the calculation of the post-test probabilities.
#'   If pretest.prob = NULL, the sample prevalence is used.
#' @param roll.length (default = NULL) The frame length of the interval of test
#'   scores scores for the calculation of the reliable predictive values. When
#'   NULL, it is calculated as \code{round(SEM)*2+1} (approximately a 68\%
#'   confidence interval around the observed test score in which the true score
#'   is expected). The roll.length needs to be uneven. When the result is even,
#'   you need to choose a value 1 larger or smaller. When roll.length is
#'   supplied, this is used instead of the calculated value. Furthermore,
#'   roll.length has to be >= 1 and uneven. Applicable values are 1, 3, 5, etc.
#'   roll.length = 1 causes the reliability of the test to be ignored.
#' @param extend (default = TRUE) The Reliable Predictive Values cannot be
#'   calculated for the most extreme scores. As the most extreme scores offer
#'   most often least uncertain decisions for or against the disease, the values
#'   that can be calculated are extended. When extend = FALSE, NAs (NOT
#'   Availables) are produced.
#' @param decision.odds (default = 2). The odds for a positive or negative
#'   classification decision. The limit for the Predictive Values =
#'   decision.odds / (decision.odds+1). For a decision, the Predictive Values
#'   needs to be larger than this limit. decision.odds <= 1 causes all test
#'   scores to be used for either positive or negative decisions. NB 1
#'   Decison.odds can be a broken number, such as .55/(1-.55), which defines the
#'   decision limit for predictive values as .55. NB 2 When a test is more
#'   reliable and valid, a higher value for decision.odds can be applied. NB 3
#'   For serious diseases with relatively uncomplicated cures, decision odds can
#'   be smaller than one. In that case, a large number of false positives is
#'   unavoidable and positive decisions are inherently uncertain. See
#'   Sonis(1999) for a discussion.
#' @param decision.use (default = 'standardized.pv'). The probability to be used
#'   for decisions. When 'posttest.probability' is chosen, pt.prob is used for
#'   positive decisions and (1 \- pt.prob) is used for negative decisions. When
#'   'predictive.value' is chosen the rnpv are used for negative decisions and
#'   the rppv for positive decisions. When 'standardized.pv' is chosen, the
#'   standardized positive predictive value is used for positive decisions and
#'   the standardized negative predictive value is used for negative decisions.
#'   N.B. These parameters can be abbreviated as 'post', 'pred' and 'dens'. N.B.
#'   The posttest probability is equal to the positive predictive value when
#'   pre-test probability = sample prevalence. N.B. The posttest probability is
#'   equal to the standardized positive predictive value when pre-test
#'   probability = .5.
#' @param preselected.thresholds (default = c(NULL, NULL)). For use in
#'   comparisons, when preselected.thresholds has valid values these values are
#'   used for the determination of the cut-points. The two cut-points indicate
#'   the limits of the uncertain area. Parameter decision.use is ignored. When
#'   preselected.thresholds[2] > preselected.thresholds[1], the higher scores
#'   are used for positive decisions. positive decisions: test values >
#'   preselected.thresholds[2] negative decisions: test values <
#'   preselected.thresholds[1]) uncertain decisions: test values >=
#'   preselected.thresholds[1]  & <= preselected.thresholds[2]. When
#'   preselected.thresholds[2] < preselected.thresholds[1], the lower scores are
#'   used for positive decisions. positive decisions: test values <=
#'   preselected.thresholds[2] negative decisions: test values >
#'   preselected.thresholds[1]) uncertain decisions: test values >
#'   preselected.thresholds[2]  & <= preselected.thresholds[1]. A single
#'   threshold can be determined with in-between values. For instance c(1.5,
#'   1.3) defines a single threshold of 1 for a descending test.
#' @param digits   the number of digits used in the output.
#'
#' @details This function can be applied to ordinal data. Uncertain test scores
#'   are scores that have about the same density in the two distributions of
#'   patients with and without the targeted condition. This range is typically
#'   found around the optimal cut-point, that is, the point of intersection or
#'   Youden index (Schisterman et al., 2005). The ui-functions use a specified
#'   low value for Se and Sp to find this uncertain interval (default Se = Sp =
#'   .55). Instead, this function uses as a default the decision odds of ordinal
#'   test scores near 1 (default < 2). The limit for the Predictive Values =
#'   decision.odds / (decision.odds+1).
#'
#'   N.B. 1: Sp = Negative Decisions | true.neg.status; Se = Positive Decisions
#'   | true.pos.status. Please note that the values for Se and Sp are
#'   underestimated, as the uncertain test scores are considered here as errors,
#'   which they are not. (Se and Sp are dichotomous indices.). Use
#'   \code{\link{quality.threshold}} and
#'   \code{\link{quality.threshold.uncertain}} for obtaining respectively
#'   quality indices for the test scores when ignoring test scores in the
#'   uncertain interval and the qualitry indices of the test scores within the
#'   uncertain interval.
#'
#'   N.B. 2: For the category Uncertain the odds are (sum of patients with a
#'   positive.status)/(sum of patients with negative.status).
#'
#'   N.B. 3: Set roll.length to 1 to ignore the test reliability and obtain raw
#'   predictive values, likelihood ratios, etc., that are not corrected for the
#'   unreliability of the test.
#'
#'   Raw predictive values compare the frequencies and provide exact sample
#'   values and are most suitable for evaluating the sample results. When
#'   prevalence is low, Positive Predictive Values can be disappointingly low,
#'   even for tests with high Se values. When prevalence is high, Negative
#'   Predictive Values can be low. Reliable Standardized Predictive Values
#'   compare the densities (relative frequencies) and are most suitable for
#'   comparing the two distributions of the scores for patients with and without
#'   the targeted condition.
#'
#'   The predictive values are calculated from the observed frequencies in the
#'   two samples of patients with and without the targeted disease. For a range
#'   of test scores x, if f0(x) and f1(x) are the frequencies of respectively
#'   patients without and with the targeted disease, then the negative
#'   predictive value (NPV) can be defined as: NPV(x) = f0(x) / (f0(x) + f1(x))
#'   and the positive predictive value (PPV) as: PPV(x) = f1(x) / (f0(x) +
#'   f1(x)). The densities for a range of test scores x can be defined d0(x) =
#'   f0(x) / n0 and d1(x) = f1(x) / n1, where n0 and n1 are the number of
#'   observed patients in the two samples.  The standardized negative predictive
#'   value (SNPV) is defined as SNPV(x) = d0(x) / (d0(x) + d1(x)) and the
#'   standardized positive predictive value (SPPV) as SPPV(x) = d1(x) / (d0(x) +
#'   d1(x)). The two distributions are weighed equally, or in other words, the
#'   prevalence is standardized to .5. N.B. The posttest probability is equal to
#'   the positive predictive value when the pretest probability is set to the
#'   sample prevalence, while the standardized positive predictive value is
#'   equal to the posttest probability when the pretest probability is set to
#'   .5.
#'
#'   Reliable estimates of the predictive probabilities correct to a certain
#'   degree for random variations. In test theory this random effect is
#'   estimated with the Standard Error of Measurement (SEM), which is directly
#'   dependent on the reliability of the test: \code{SEM = s * sqrt(1 - r)},
#'   where s is the standard deviation of the test scores and r the estimated
#'   reliability of the test (Crocker & Algina, 1986; Harvill, 1991). The true
#'   score of a patient lies with some probability (roughly 68%) within a range
#'   of +- 1 SEM around the acquired test score. This provides information about
#'   the range of test scores that can be expected due to all kinds of random
#'   circumstances when no real changing agent is effective.
#'
#'   The results show the obtained values for the sample and are not corrected
#'   in any way. The classification 'Uncertain' shows the scores that lead to
#'   odds (d1(x) / d0(x)) that are lower than limit. This indicates that it is
#'   difficult to base classifications on that range of scores. The positive
#'   classifications are less error prone, with realized odds (d1(x) / d0(x)).
#'   These odds are close to 1 and smaller than the decision.odds. The negative
#'   classifications are less error prone than 'Uncertain' (odds = d0(x) /
#'   d1(x)).
#'
#'   The accuracy indices are shown as percentages. Sp = negative
#'   classifications given a true negative status. Se = positive classifications
#'   given a true positive status. NPV = proportion of Negative Classifications
#'   that are correct. PPV = proportion of Positive Classifications that are
#'   correct.
#' @return{ A list of: } \describe{ \item{$parameters: }{ A named vector:
#' \itemize{ \item{pretest.prob: }{provided or calculated pre-test probability.
#' Default, the calculated sample prevalence is used.} \item{sample.prevalence:
#' }{the calculated sample prevalence.} \item{reliability: }{must be provided;
#' ignored when roll.length = 0.} \item{SEM: }{the calculated Standard Error of
#' Measurement (SEM).} \item{roll.length: }{the total length of the range around
#' the test score (2 * SEM + 1).} \item{rel.conf.level: }{the confidence level
#' of the range, given the reliability.} \item{decision odds: }{the provided
#' level of decision certainty.} \item{limit: }{the limit applied to the
#' predictive values for calculating the result.}}} \item{$messages: }{Two
#' messages: 1. The test scores are reported for which reliable predictive
#' values could not be calculated and have been extended from the nearest
#' calculated value, 2. the kind of probabilities that are used for decisions.}
#' \item{$rel.pred.values: }{A table the test scores as columns and with rows:
#' N.B. When roll.length is set to 1, the test reliability is ignored and the
#' outcomes are not corrected for unreliability. \itemize{ \item{rnpv: }{(more)
#' reliable negative predictive value. Fitting for reporting sample results. }
#' \item{rppv: }{(more) reliable positive predictive value.} \item{rndpv:
#' }{(more) reliable standardized negative predictive value.} \item{rpdpv:
#' }{(more) reliable standardized positive predictive value.} \item{rilr:
#' }{(more) reliable interval likelihood ratio.} \item{rpt.odds: }{(more)
#' reliable posttest odds.} \item{rpt.prob: }{(more) reliable posttest
#' probabilities.} } }
#'
#' \item{$result: }{Table of results for the current sample, calculated with the
#' provided parameters. \itemize{ \item{columns: }{Negative Classifications,
#' Uncertain, Positive Classifications.} \item{row scores: }{range of test
#' scores for the three categories.} \item{row total.sample: }{percentage of the
#' total sample.} \item{row correct.decisions: }{percentages of correct negative
#' and positive decisions (NPV and PPV).} \item{row true.neg.status:
#' }{percentage of patients with a true negative status for the 3 categories.}
#' \item{row true.pos.status: }{percentage of patients with a true positive
#' status for the 3 categories.} \item{row realized.odds: }{The odds that are
#' realized in the sample for each of the three categories.} } } }
#' @references Sonis, J. (1999). How to use and interpret interval likelihood
#'   ratios. Family Medicine, 31, 432–437.
#'
#'   Crocker, L., & Algina, J. (1986). Introduction to classical and modern test
#'   theory. Holt, Rinehart and Winston, 6277 Sea Harbor Drive, Orlando, FL
#'   32887 ($44.75).
#'
#'   Harvill, L. M. (1991). Standard error of measurement. Educational
#'   Measurement: Issues and Practice, 10(2), 33–41.

#' @examples
#' set.seed(1)
#' # example of a validation sample
#' ref=c(rep(0,1000), rep(1, 1000))
#' test=round(c(rnorm(1000, 5, 1), rnorm(1000, 8, 2)))
#' # calculated roll.length is invalid. Set to 3. Post test probability equals
#' # Positive Predictive Values. Parameter pretest.prob is set to sample prevalence.
#' RPV(ref, test, reliability = .9, roll.length = 3)
#' # Set roll.length = 1 to ignore test reliability (value of parameter
#' # reliability is ignored, but must be set to some value.)
#' RPV(ref, test, reliability = 0, roll.length = 1)
#' # When pretest.prob is set to .5, the Post-test Probabilities are equal to
#' # the Standardized Positive Predictive Values.
#' RPV(ref, test, pretest.prob = .5, reliability = .9, roll.length = 3)
#'
# pretest.prob = .5; reliability=.9; roll.length = 1; extend = TRUE;
# decision.odds = 3; decision.use = 'LR'; digits=3; preselected.thresholds=NULL
# preselected.thresholds=c(25.9, 25.8)
RPV <- function(ref, test, pretest.prob = NULL, reliability, roll.length = NULL,
                extend = TRUE, decision.odds = 2, decision.use =
                  c('standardized.pv', "posttest.probability", 'LR',
                    "predictive.value"),
                preselected.thresholds=c(NULL, NULL), digits=3){

  # require(zoo)
  df = check.data(ref,test, 'ordinal')
  ref = df$ref
  test = df$test
  decision.use = match.arg(decision.use)


  is.even <- function(x) x %% 2 == 0 # is.even(F)
  percent <- function(x, digits = digits, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  mktxtseq <- function(r) {
    d = abs(diff(c(-Inf, r)))
    w = which(d != 1)
    ranges = list()
    for (i in (1:(length(w) - 1))) {
      ranges[i] = paste(r[w[i]], r[(w[i + 1] - 1)], sep = '-')
    }
    ranges[i + 1] = paste(r[w[i + 1]], r[length(r)], sep = '-')
    paste(unlist(ranges), collapse=' ')
  }

  SEM = sd(test)*sqrt(1-reliability)
  if (is.null(roll.length)) roll.length = round(SEM*2+1) # roll.length=NA
  if (roll.length < 1) stop('Provided roll.length is smaller than 1.')
  if (is.even(roll.length)){
    stop(paste('Calculated or provided roll.length is even. Please define roll.length as',
                 roll.length-1, 'or', roll.length+1))
  }

  # roll.length=1
  if (is.numeric(test)) test = factor(test, ordered=TRUE, levels = min(test):max(test)) # is.ordered(test)
  if (is.numeric(ref)) ref = factor(ref, ordered=T, levels = min(ref):max(ref))
  tt = table(ref, test)
  ts.npv = rollsum(tt[1,], roll.length, fill=NA)/(rollsum(tt[2,]+tt[1,],roll.length, fill=NA))
  ts.ppv = 1- ts.npv
  st0 = sum(tt[1,])
  st1 = sum(tt[2,])
  sample.prevalence = st1/(st0+st1)

  if (is.null(pretest.prob)) pretest.prob = sample.prevalence
  ts.ilr = (rollsum(tt[2,], roll.length, fill=NA)*st0)/(rollsum(tt[1,],roll.length, fill=NA)*st1)
  #
  ts.pdpv =  (rollsum(tt[2,], roll.length, fill=NA)/st1)/
    ((rollsum(tt[2,], roll.length, fill=NA)/st1)+(rollsum(tt[1,],roll.length, fill=NA)/st0))
  ts.ndpv =  (rollsum(tt[1,], roll.length, fill=NA)/st0)/
    ((rollsum(tt[1,], roll.length, fill=NA)/st0)+(rollsum(tt[2,],roll.length, fill=NA)/st1))

  if (is.null(preselected.thresholds)){
    decuse = paste('Decision use = ', decision.use, '.', sep='')
  } else {
    decuse = 'Decision use = preselected.thresholds.'
  }
  if (extend & any(is.na(ts.npv))){
    ext = rbind(paste('Reliable Predictive Values for scores ', paste(names(which(is.na(ts.npv))), collapse=' '), ' have been extended.'),
          unname(decuse))
    ts2.npv=na.fill(ts.npv, fill=c('extend', NA, 'extend'))
    ts2.ppv=na.fill(ts.ppv, fill=c('extend', NA, 'extend'))
    ts2.ndpv=na.fill(ts.ndpv, fill=c('extend', NA, 'extend'))
    ts2.pdpv=na.fill(ts.pdpv, fill=c('extend', NA, 'extend'))
    ts2.ilr=na.fill(ts.ilr, fill=c('extend', NA, 'extend'))
  } else {
    ext = rbind('No extension has been applied.', unname(decuse))
    ts2.npv=ts.npv
    ts2.ppv=ts.ppv
    ts2.ndpv=ts.ndpv
    ts2.pdpv=ts.pdpv
    ts2.ilr=ts.ilr
  }
  pretest.odds = pretest.prob/(1-pretest.prob)
  posttest.odds = pretest.odds * ts2.ilr # posttest.odss = Inf
  posttest.prob = ifelse (is.infinite(posttest.odds) & posttest.odds > 0, 1,
      posttest.odds/(posttest.odds+1))


  rel.conf.level = pnorm((roll.length-1)/2, 0, SEM) - pnorm(-(roll.length-1)/2, 0, SEM)

  if (decision.use == 'LR') {
    limit = decision.odds / pretest.odds   # limit is likelihood ratio
  } else {
    limit = decision.odds/(decision.odds+1) # probability
  }

  arg1 = c(round(c(pretest.prob = pretest.prob, sample.prevalence=sample.prevalence,
                   reliability=reliability, SEM=SEM, roll.length=roll.length,
                 rel.conf.level=rel.conf.level, decision.odds=decision.odds,
                 limit=limit), digits)) # as.numeric(arg1[1])*2

  scores = 1:ncol(tt); names(scores)=colnames(tt)
  if (is.null(preselected.thresholds)){
    if (decision.use == 'predictive.value'){
      pd = which(ts2.ppv > limit)
      if (decision.odds < 1){
        nd = scores[-pd]
        ud = NA
      } else {
        nd = which(ts2.npv > limit)
        ud = which(ts2.ppv <= limit & ts2.npv <= limit)
      }
    } else if (decision.use == 'posttest.probability'){
      pd = which(posttest.prob > limit)
      if (decision.odds < 1){
        nd = scores[-pd]
        ud = NA
      } else {
        nd = which((1-posttest.prob) > limit)
        ud = which(posttest.prob <= limit & ((1-posttest.prob) <= limit)) }
    } else if (decision.use == 'standardized.pv'){
      if (decision.odds < 1){
        pd = which(ts2.pdpv > limit)
        nd = scores[-pd]
        ud = NA
      } else {
        pd = which(ts2.pdpv > limit)
        nd = which(ts2.ndpv > limit)
        ud = which(ts2.ndpv <= limit & (ts2.pdpv <= limit))
      }
    } else if (decision.use == 'LR'){
      if (decision.odds < 1){
        pd = which(ts2.ilr > limit)
        nd = scores[-pd]
        ud = NA
      } else {
        pd = which(ts2.ilr > limit)
        nd = which(ts2.ilr < 1/limit)
        ud = which(ts2.ilr <= limit & (ts2.ilr >= 1/limit))
      }
    }
    } else { # preselected.thresholds = c(25.2,25.1)
      numlevels = as.numeric(levels(test))
      if (preselected.thresholds[2] > preselected.thresholds[1]){ # positive decisions for higher scores
        pd = which(numlevels > preselected.thresholds[2])
        nd = which(numlevels < preselected.thresholds[1])
        ud =  which(numlevels >= preselected.thresholds[1]  & (as.numeric(levels(test)) < preselected.thresholds[2]))
      } else { # positive decisions for lower scores
        pd = which(numlevels < preselected.thresholds[2])
        nd = which(numlevels > preselected.thresholds[1])
        ud =  which(numlevels >= preselected.thresholds[2]  & (as.numeric(levels(test)) <= preselected.thresholds[1]))
      }
      names(pd)=colnames(tt)[pd]
      names(nd)=colnames(tt)[nd]
      names(ud)=colnames(tt)[ud]
    }

  pv = as.matrix(rbind(rnpv=round(ts2.npv, digits), rppv = round(ts2.ppv, digits),
                      rndpv=round(ts2.ndpv, digits), rpdpv=round(ts2.pdpv, digits),
                      rilr=round(ts2.ilr,digits), rpt.odds = round(posttest.odds,digits), rpt.prob = round(posttest.prob, digits)))
  if (roll.length==1) row.names(pv) = c('npv', 'ppv', 'ndpv', 'pdpv', 'ilr', 'pt.odds', 'pt.prob')

  col.title = c('Negative Decisions', 'Uncertain', 'Positive Decisions')
  # scores = c(paste(names(nd), collapse = " "), paste(names(ud),collapse = " "), paste(names(pd), collapse = " "))
  scores = c(mktxtseq(as.numeric(names(nd))),
             mktxtseq(as.numeric(names(ud))),
             mktxtseq(as.numeric(names(pd))))
  n = c(sum(tt[,nd]), sum(tt[,ud]), sum(tt[,pd]))
  total.sample =  percent(c(sum(tt[,nd])/sum(tt), sum(tt[,ud])/sum(tt),sum(tt[,pd])/sum(tt)), 1)
  correct.decisions = percent(c(sum(tt[1,nd])/sum(tt[,nd]), NA, sum(tt[2,pd])/sum(tt[,pd])), 1)
  tns = c(sum(tt[1,nd])/sum(tt[1,]),sum(tt[1,ud])/sum(tt[1, ]), sum(tt[1,pd])/sum(tt[1,]))
  tps = c(sum(tt[2,nd])/sum(tt[2,]),sum(tt[2,ud])/sum(tt[2, ]), sum(tt[2,pd])/sum(tt[2,]))
  true.neg.status = percent(tns, 1)
  true.pos.status = percent(tps, 1)
  ud.odds = sum(tt[2,ud])/sum(tt[1,ud])
  realized.odds = round(c(sum(tt[1,nd])/sum(tt[2,nd]), ud.odds, sum(tt[2,pd])/sum(tt[1,pd])), digits)
  tt2 = as.table(rbind(scores, n, total.sample, correct.decisions, true.neg.status, true.pos.status, realized.odds))
  colnames(tt2) = col.title

  out = list(parameters = arg1, messages = ext, rel.pred.values=pv, result=tt2)
  return(out)
}

# Example with two intersections, i.e. two optimal cut-points
# nobs=1000; set.seed(88);
# test= round(c(rnorm(nobs), rnorm(.5*nobs, 2, 1))); ref= c(rep(0,nobs), rep(1, .5*nobs))
# reliability=.9; roll.length = 5; extend = TRUE; decision.odds = 2
# RPV(ref, test, reliability, roll.length=roll.length, decision.odds = decision.odds)
# (yt = youden.threshold(ref, test))
# quality.threshold(ref, test, yt)
# RPV(ref, test, 0.8, decision.odds = 2)
# RPV(ref, test, 0.8, decision.odds = 3)
# RPV(ref, test, 0.9, decision.odds = 1)
# RPV(ref, test, 0.9, decision.odds = 2)
# RPV(ref, test, 0.9, decision.odds = .7/.3)
# RPV(ref, test, 0.9, roll.length=1, decision.odds = 1)
# RPV(ref, test, 0.9, roll.length=1, decision.odds = 2)
# RPV(ref, test, 0.9, roll.length=1, decision.odds = 3)
# RPV(ref, test, 0.7, roll.length = 5, decision.odds = 2)
# RPV(ref, test, 0.7, roll.length = 1, decision.odds = 1)
# RPV(ref, test, 0.9, roll.length=3, decision.odds = 3)
# RPV(ref, test, 0.9, roll.length=5, decision.odds = 3)
