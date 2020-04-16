#' Function to explore possible uncertain intervals of ordinal test results of
#' individuals with (1) and without (0) the targeted condition.
#' @name ui.ordinal
#' @importFrom grDevices rgb
#' @importFrom graphics legend plot
#' @description This function explores possible uncertain intervals (UI) of the
#'   test results of the two groups. This functions allows for considerable
#'   fine-tuning of the characteristics of the interval of uncertain test
#'   scores, in comparison to other functions for the determination of the
#'   uncertain interval and is intended for tests with a limited number of
#'   ordered values and/or small samples.
#'
#'   This function is intended to be used for tests with 20 or less ordered test
#'   values. The lower range of test scores identifies patients without the
#'   targeted condition (lower More Certain Interval (MCI)), the upper interval
#'   of test scores above the uncertain interval identifies the patients with
#'   the condition (upper MCI). Due to the limited number of distinguishable
#'   scores, the estimations are course. When more than 20 values can be
#'   distinguished, \code{\link{ui.nonpar}} or \code{\link{ui.binormal}} may be
#'   preferred. When a sufficiently large dataset is available, the function
#'   \code{\link{RPV}} may be preferred for the analysis of discrete ordered data.
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the classification by the reference test. The reference standard
#'   must be coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition).
#' @param test The test or predictor under evaluation. A column in a dataset or
#'   vector indicating the test results on an ordinal scale. It is expected that
#'   true patients have higher scores than non-patients. If this is not the
#'   case, the test scores should be negated (test = -(test scores)).
#' @param select.max Selects the candidate thresholds on basis of a desired
#'   property of the More Certain Intervals (MCI). The criteria are: maximum
#'   Se+Sp (default), maximum C (AUC), maximum Accuracy, maximum Sp, maximum Se,
#'   maximum size of MCI. The last alternative 'All' is to choose all possible
#'   details.
#' @param constraints Sets upper constraints for various properties of the
#'   uncertain interval: C-statistic (AUC), Acc (accuracy), lower and upper
#'   limit of the ratio of the proportions with and without the targeted
#'   condition. The default values are C = .57, Acc = .6, lower.ratio = .8,
#'   upper.ratio = 1.25. These values implement the desired uncertainty of the
#'   uncertain interval. The value of C (AUC) is considered the most important
#'   and has the most restrictive default value. For Acc and C, the values
#'   closest to the desired value are found and then all smaller values are
#'   considered. The other two constraints are straightforward lower and upper
#'   limits of the ratio between the number of patients with and without the
#'   targeted disease. If you want to change the values of these constraints, it
#'   is necessary to name all values. C = 1 or Acc = 1 excludes C respectively
#'   accuracy as selection criterion. If no solution is found, the best is
#'   showed together with  a warning message.
#' @param weights (Default = c(1, 1, 1). Vector with weights for the loss
#'   function. weights[1] is the weight of false negatives, weights[2] is the
#'   weight for loss in the uncertain interval (deviations from equal chances to
#'   belong to either distribution), and weights[3] is the weight for false
#'   positives. When a weight is set to a larger value, thresholds are selected
#'   that make the corresponding error smaller while the area grows smaller.
#' @param intersection (Default = NULL). Optional value to de used as value for
#'   the intersection. If no value is supplied, the intersection is calculated
#'   using  \code{get.intersection(ref = ref, test = test, model='ordinal').}
#' @param return.all (Default = FALSE). When TRUE $data.table and
#'   $uncertain.interval are included in the output.
#' @param ... Further parameters that can be transferred to the density
#'   function.
#' @details Due to the limited possibilities of short scales, it is more
#'   difficult to determine a suitable uncertain interval when compared to
#'   longer scales. This problem is aggrevated when samples are small. For any
#'   threshold determination, one needs a large representative sample (200 or
#'   larger). If there are no test scores below the intersection in the
#'   candidate uncertain area, Sp of the Uncertain Interval (UI.Sp) is not
#'   available, while UI.Se equals 1. The essential question is always whether
#'   the patients with the test scores inside the uncertain interval can be
#'   sufficiently distinguished. The candidate intervals are selected on various
#'   properties of the uncertain interval. The defaults are C (AUC) lower than
#'   .6, Acc lower than .6, and the ratio of proportions of persons with /
#'   without the targeted condition between .8 and 1.25. These criteria ensure
#'   that all candidates for the uncertain interval have insufficient accuracy.
#'   The second criterion is the desired property of the More Certain Intervals
#'   (see select.max parameter). The model used is 'ordinal'. This model default
#'   for the adjust parameter send to the density function is 2, but you can
#'   enter another value such as adjust = 1.
#'
#'   Discussion of the first example (please run the code first): Visual
#'   inspection of the mixed densities function \code{\link{plotMD}} shows that
#'   distinguishing patients with and without the targeted condition is almost
#'   impossible for test scores 2, 3 and 4. Sensitivity and Specificity of the
#'   uncertain interval should be not too far from .5. In the first example, the
#'   first interval (3:3) has no lower scores than the intersection (3), and
#'   therefore UI.Sp is not available and UI.Se = 1. The UI.ratio indicates
#'   whether the number of patients with and without the condition is equal in
#'   this interval. For these 110 patients, a diagnosis of uncertainty is
#'   probably the best choice. The second interval (3:4) has an UI.Sp of .22,
#'   which is a large deviation from .5. In this slightly larger interval, the
#'   patients with a test score of 3 have a slightly larger probability to
#'   belong to the group without the condition. UI.Se is .8. UI.ratio is close
#'   to 1, which makes it a feasible candidate. The third interval (2:4) has an
#'   UI.Sp of .35 and an UI.Se of .70 and an UI.ratio still close to one. The
#'   other intervals show either Se or Sp that deviate strongly from .5, which
#'   makes them unsuitable choices. Probably the easiest way to determine the
#'   uncertain interval is the interval with minimum loss. This is interval
#'   (2:4). Dichotomization loss L2 can be defined as the sum of false negatives
#'   and false positives. The Youden threshold minimizes these. The Loss formula
#'   L3 for trichotomization of ordinal test scores is (created by
#'   https://www.codecogs.com/latex/eqneditor.php): \deqn{L_3 =\frac{ \left
#'   (\sum_{i=l}^{u} \left |d0_{i}-d1_{i}  \right | + \sum_{i=u+1}^{h} d1_{i}+
#'   \sum_{i=1}^{l-1}d0_{i}\right )}{N}}{ L3 = 1/N * (sum(abs(d0[u:l] -
#'   d1[u:l])) + sum(d1[1:(l-1)]) + sum(d0[(u+1):h]))}
#'   where \emph{d0}
#'   represents the test scores of the norm group, \emph{d1} represents the test
#'   scores of the targeted patient group, \emph{l} is the lower limit of the
#'   uncertain interval, \emph{u} the upper limit, the first test score is
#'   enumerated 1 and the last test score is enumerated \emph{h}. \emph{N} is
#'   the total number of all persons with test scores.
#'   \itemize{
#'   \item{\eqn{\sum_{i=l}^{u} \left |d0_{i}-d1_{i}  \right |}{sum(abs(d0[u:l] -
#'   d1[u:l])}}{ is the loss in the uncertain interval, that is, the total
#'   deviation from equality.} \item{\eqn{\sum_{i=u+1}^{h}
#'   d1_{i}}{sum(d1[1:(l-1)])}}{ is the loss in the lower More Certain Interval,
#'   that is, the total of False Negatives, the number of patients with the
#'   targeted condition with a test score lower than \emph{l}, and}
#'   \item{\eqn{\sum_{i=u+1}^{h} d0_{i}}{sum(d0[(u+1):h])}}{ is the loss in the
#'   upper More Certain Interval, that is, the total of False Positives, the
#'   number of patients without the targeted condition with a test score higher
#'   than \emph{u}.}}
#'
#'   Loss L is higher when the deviation from equality is higher in the
#'   uncertain area, higher when the number of False Negatives is higher, and
#'   higher when the number of False Positives is higher. The loss of a single
#'   threshold method equals 1 - its Accuracy. In this example, the minimum Loss
#'   is found with interval (2:4). As this agrees with values for UI.C and
#'   UI.ratio that sufficiently indicates the uncertainty of these test scores,
#'   this seems the most suitable choice: the number of patients with test
#'   scores 2 to 4 are almost as likely to come from either population. The
#'   remaining cases outside the uncertain interval (2:4) show high C, Accuracy,
#'   Specificity and Sensitivity.
#'
#' @return List of values: 
#' \describe{ 
#' \item{$Youden}{A vector of statistics
#'   concerning the maximized Youden index:} 
#'   \itemize{ 
#'   \item{max.Youden: }{The
#'   value of the Maximized Youden Index (= max(tpr - fpr)).} 
#'   \item{threshold:
#'   }{The threshold associated with the Maximized Youden Index. Test values >=
#'   threshold indicate the targeted condition.} 
#'   \item{Sp: }{The Specificity of
#'   the test when this threshold is applied.} 
#'   \item{Se: }{The Sensitivity of
#'   the test when this threshold is applied.} 
#'   \item{Acc: }{The Accuracy of the
#'   test when this threshold is applied.} 
#'   \item{Loss: }{min(fnr + fpr) = min(1
#'   - (Se + Sp -1)) = 1 - max(tpr - fpr) lower range ( < threshold): the summed
#'   number of false positives for each test score, divided by the number of
#'   persons that have received that test score. upper range ( >= threshold):
#'   the summed number of false negatives, divided by the number of persons that
#'   have received that test score. The Youden Loss is equal to 1-Youden.index.
#'   } \
#'   item{C: }{Concordance; equals AUROCC (Area Under Receiving Operating
#'   Characteristics Curve or AUC)} } 
#'   \item{$data.table}{A data.frame with the
#'   following columns:} 
#'   \itemize{ 
#'   \item{test: }{The test scores.} 
#'   \item{d0:
#'   }{The frequencies of the test scores of the norm group.} 
#'   \item{d1: }{The
#'   frequencies of the test scores of the group with the targeted condition.}
#'   \item{tot: }{The total frequency of each test scores.} 
#'   \item{TP: }{The
#'   number of True Positives when this test score is used as threshold.}
#'   \item{FP: }{The number of False Positives when this test score is used as
#'   threshold.} 
#'   \item{tpr: }{The true positive rate when this test score is
#'   used as threshold.} 
#'   \item{fpr: }{The false positive rate when this test
#'   score is used as threshold.} 
#'   \item{Y: }{The Youden Index (= tpr - fpr) when
#'   this test score is used as threshold.} } 
#'   \item{$intersection}{The (rounded)
#'   intersection for the distributions of the two groups. Most often, these
#'   distributions have no true point of intersection and the rounded
#'   intersection is an approximation. Often, this equals the Maximized Youden
#'   threshold (see Schisterman 2005). Warning: When a limited range of scores
#'   is available, it is more difficult to estimate the intersection. Different
#'   estimates can easily differ plus minus 1. When using a non-rounded value
#'   (for example 16.1), the effective threshold for the uncertain area is
#'   round(intersection+.5), in the mentioned example: 16.1 becomes 17. }
#'   \item{$uncertain.interval}{Data frame with the statistics of all possible
#'   bounds of the uncertain interval. The columns are the following: }
#'   \itemize{ 
#'   \item{lowerbound: }{Lower bound of the possible uncertain
#'   interval.} 
#'   \item{upperbound: }{Upper bound of the possible uncertain
#'   interval.} 
#'   \item{UI.Sp: }{Specificity of the test scores between and
#'   including the lower and upper boundary. Closer to .5 is 'better', that is,
#'   more uncertain. This estimate is rough and dependent on the intersection
#'   and cannot be recommended as a criterion for a short, ordinal scale. }
#'   \item{UI.Se: }{Sensitivity of the test scores between and including the
#'   lower and upper boundary. Closer to .5 is 'better', that is, more
#'   uncertain. This estimate is rough and dependent on the intersection and
#'   cannot be recommended as a criterion for a short, ordinal scale.}
#'   \item{UI.Acc: }{Accuracy of the test scores between and including the lower
#'   and upper boundary. Closer to .5 is 'better', that is, more uncertain. This
#'   estimate is rough and dependent on the intersection and cannot be
#'   recommended as a criterion for a short, ordinal scale.} 
#'   \item{UI.C:
#'   }{Concordance (AUROC) of the test scores between and including the lower
#'   and upper boundary. Closer to .5 is 'better', that is, more uncertain. Rule
#'   of thumb: <= .6} 
#'   \item{UI.ratio: }{The ratio between the proportion of
#'   patients in the uncertain area with and without the condition. Closer to
#'   one is 'better', that is, more uncertain; 0.8 < UI.ratio < 1.25 as a rule
#'   of fist.} 
#'   \item{UI.n: }{Number of patients with test scores between and
#'   including the lower and upper boundary.} 
#'   \item{MCI.Sp: }{Specificity of the
#'   more certain interval, i.e., the test scores lower than the lower boundary
#'   and higher than the upper boundary.} 
#'   \item{MCI.Se: }{Sensitivity of the
#'   test scores lower than the lower boundary and higher than the upper
#'   boundary.} 
#'   \item{MCI.C: }{Concordance (AUROC) of the test scores outside
#'   the uncertain interval. Closer to .5 is 'better', that is, more uncertain.
#'   Rule of thumb: <= .6} 
#'   \item{MCI.Acc: }{Accuracy of the test scores lower
#'   than the lower boundary and higher than the upper boundary.} \item{MCI.n:
#'   }{Number of patients with test scores lower than the lower boundary and
#'   higher than the upper boundary.} 
#'   \item{Loss: }{Loss of the
#'   trichotomization. The total loss is the sum of the loss of the three areas:
#'   lower MCI: the summed number of false positives for each test score,
#'   divided by the number of persons that have received that test score.
#'   uncertain interval: the sum of the absolute differences in the number of
#'   people in the norm group d0 and the number of persons in the group with the
#'   targeted condition (d1) per test score, divided by the total number of
#'   persons.} upper MCI: the summed number of false negatives, divided by the
#'   number of persons that have received that test score. The Loss can be
#'   compared to the loss of the Youden threshold, provided that the
#'   intersection is equal to the Youden threshold. If necessary, this can be
#'   forced by attributing the value of the Youden threshold to the intersection
#'   parameter. } 
#'   \item{$candidates: }{Candidates with a loss lower than the
#'   Youden loss which might be considered for the Uncertain Interval. The
#'   candidates are selected based on the constraints parameter, that defines
#'   the desired constraints of the uncertain area, and the select.max
#'   parameter, that selects the desired properties of the lower and upper More
#'   Certain Interval. } }
#'
#' @references{ Youden, W. J. (1950). Index for rating diagnostic tests. Cancer,
#' 3(1), 32-35.
#' https://doi.org/10.1002/1097-0142(1950)3:1<32::AID-CNCR2820030106>3.0.CO;2-3
#'
#' Schisterman, E. F., Perkins, N. J., Liu, A., & Bondell, H. (2005). Optimal
#' cut-point and its corresponding Youden Index to discriminate individuals
#' using pooled blood samples. Epidemiology, 73-81.
#'
#' Landsheer, J. A. (2016). Interval of Uncertainty: An alternative approach for
#' the determination of decision thresholds, with an illustrative application
#' for the prediction of prostate cancer. PLOS One.
#'
#' Landsheer, J. A. (2018). The Clinical Relevance of Methods for Handling
#' Inconclusive Medical Test Results: Quantification of Uncertainty in Medical
#' Decision-Making and Screening. Diagnostics, 8(2), 32.
#' https://doi.org/10.3390/diagnostics8020032 }
#' @seealso{  \code{\link{plotMD}} or \code{\link{barplotMD}} for plotting the
#' mixed densities of the test values. \code{\link[stats]{density}} for the
#' parameters of the density function. } \code{\link{ui.nonpar}} or
#' \code{\link{ui.binormal}} can be used when more than 20 values can be
#' distinguished on the ordinal test scale. When a large data set for an ordinal
#' test is available, one might consider \code{\link{RPV}}.
#' @examples
#' # A short test with 5 ordinal values
#' test0     = rep(1:5, times=c(165,14,16,55, 10)) # test results norm group
#' test1     = rep(1:5, times=c( 15,11,13,55,164)) # test results of patients
#' ref = c(rep(0, length(test0)), rep(1, length(test1)))
#' test = c(test0, test1)
#' table(ref, test)
#' plotMD(ref, test, model='ordinal') # visual inspection
#' ui.ordinal(ref, test, select.max='All')
#' # Same solution, but other layout of the results:
#' ui.ordinal(ref, test, select.max=c('MCI.Sp+MCI.Se', 'MCI.C', 'MCI.Acc',
#'                                    'MCI.Se', 'MCI.Sp', 'MCI.n'))
#' # forcing the Youden threshold as intersection gives the same best result.
#' # However, the estimates for ui.Se, ui.Sp and ui.Acc differ:
#' ui.ordinal(ref, test, intersection='Youden', select.max='All')
#'
#' nobs=1000
#' set.seed(6)
#' Z0 <- rnorm(nobs, mean=0)
#' b0=seq(-5, 8, length.out=31)
#' f0=cut(Z0, breaks = b0, labels = c(1:30))
#' x0=as.numeric(levels(f0))[f0]
#' Z1 <- rnorm(nobs, mean=1, sd=1.5)
#' f1=cut(Z1, breaks = b0, labels = c(1:30))
#' x1=as.numeric(levels(f1))[f1]
#' ref=c(rep(0,nobs), rep(1,nobs))
#' test=c(x0,x1)
#' plotMD(ref, test, model='ordinal') # looks like binormal
#' # looks less binormal, but in fact it is a useful approximation:
#' plotMD(ref, test, model='binormal')
#' ui.ordinal(ref, test)
#' ui.binormal(ref, test) # compare application of the bi-normal model

# data(tostbegg2)
# ref=tostbegg2$d; test=tostbegg2$y
# intersection=NULL; weights=c(1,1,1);
# constraints=c( Acc=.6, lower.ratio=.8, upper.ratio=1.25)
# select.max = c('MCI.Sp+MCI.Se', 'MCI.C', 'MCI.Acc', 'MCI.Se', 'MCI.Sp', 'MCI.n')
# select.max="All"

#' @export
ui.ordinal <- function(ref, test,
                       select.max = c('MCI.Sp+MCI.Se', 'MCI.C', 'MCI.Acc', 'MCI.Se', 'MCI.Sp', 'MCI.n', 'All'),
                       constraints=c(C=.57, Acc=.6, lower.ratio=.8, upper.ratio=1.25),
                       weights=c(1,1,1),
                       intersection=NULL,
                       return.all = FALSE, ...){


  # emp.AUC <- function(norm, abnorm) {
  #   o = outer(abnorm, norm, "-")
  #   mean((o > 0) + .5 * (o == 0))
  # }

  find.closest <- function(M, crit){
    M[which(is.na(M))] = Inf # crit+10e10
    mindiff=min(abs(M-crit)) # warning: no non-missing arguments to min; returning Inf
    which((M == crit+mindiff) | (M == crit-mindiff), arr.ind=T)
  }

  AUC <- function(norm, abnorm){
    n1 = as.numeric(length(abnorm));
    n0 = as.numeric(length(norm));
    r = rank(c(abnorm,norm))
    (sum(r[1:n1]) - n1*(n1+1)/2) / (n1*n0)
  }

  tlim = c(C=.57, Acc=.6, lower.ratio=.8, upper.ratio=1.25)
  if (is.null(names(constraints))) {
    tlim[1:length(constraints)] = constraints
  } else {
    if (all(names(constraints) %in% c("C" , "Acc" , "lower.ratio" , "upper.ratio"))){
    tlim[names(constraints)] = constraints[names(constraints)]}
  }
  constraints=tlim
  # if (!(exists(constraints['C'] & exists(constraints['Acc'])))) stop('constraints must be named correctly.')
  select.max <- match.arg(select.max, c('MCI.Sp+MCI.Se', 'MCI.C', 'MCI.Acc','MCI.Se', 'MCI.Sp', 'MCI.n', 'All'), several.ok = TRUE)


  tab=as.matrix(table( test, ref))
  # print(addmargins(tab))
  totpos=sum(tab[,2])                          # The total number of positives (one number)
  totneg=sum(tab[,1])                          # The total number of negatives (one number)

  d=data.frame(test=unique(sort(test)), d0=tab[,'0'], d1=tab[,'1'], row.names=1:nrow(tab))
  d$tot=rowSums(tab)                             # Number of patients w/ each test result
  d$TP=unname(rev(cumsum(rev(tab[,2]))))     # Number of true positives, start with maximum
  d$FP=unname(rev(cumsum(rev(tab[,1]))))     # Number of false positives, start with maximum
  # d$TN=unname(cumsum(tab[,1]))
  d$TN=totneg-d$FP
  d$FN=totpos-d$TP

  d$tpr=d$TP/totpos          # Sensitivity (fraction true positives)
  d$fpr=d$FP/totneg  # 1-d$FP/totneg         # 1 - specificity (false positives)
  # d$fnr=d$FN/totpos  # 1-d$FP/totneg         # 1 - specificity (false positives)
  d$Y = d$tpr-d$fpr # Youden index
  Yt= which.max(d$Y) # Youden threshold as row number
  # Hmisc::find.matches(d$Y, max(d$Y), tol=c(.0001))
  # d$Y[4] > d$Y[3]
  # which( d$Y == max(d$Y) )
  Acc.Y = (sum(d$d0[1:(Yt-1)])+sum(d$d1[(Yt):nrow(d)]))/ (totpos+totneg)
  # Conc = emp.AUC(test[ref==0], test[ref==1])
  Conc = AUC(test[ref==0], test[ref==1])

  lr = nrow(tab)
  t=prop.table(addmargins(table(test, ref),2),2)
  p0 = totneg/(totpos+totneg); p1=totpos/(totpos+totneg) # p0=p1=.5
  # r0 = p0*t[,'0']/(p0*t[,'0']+p1*t[,'1'])
  # R = cbind(r0, r1 = 1-r0)

  # A=Yt
  Loss1 = function(A) {
    lr = nrow(t)
    weights[1]*sum(t[-(A:lr), '1']) +  # FN
      weights[3]*sum(t[A:lr, '0'])   # FP
  }
  Loss.Y = Loss1(Yt) # Yt = rownumber; Loss = 1 - max.Youden

  # is = Yt # is = get.intersection(ref = ref, test = test, model='ordinal')
  # intersection='Youden'
  if (is.null(intersection)) {
    is = get.intersection(ref = ref, test = test, model='ordinal', ...)
  } else if (is.numeric(intersection)) {
    is=intersection
  } else if (intersection=="Youden"){
    is = d$test[Yt]
  }  # closest test value
  isr = round(is[length(is)]) # tail has highest density
  wr = which(d$test==isr) # row number

  dl=nrow(d)
  upperbound = c(wr:dl) # row numbers, not test scores!
  lowerbound = c(wr:1)
  ui= as.data.frame(expand.grid(lowerbound,upperbound))
  colnames(ui) <- c('lowerbound', 'upperbound')

  # x = c(1,5)
  fUI.Sp = function (x){ifelse(x[1] < wr, sum(tab[x[1]:(wr-1),1]), 0)/sum(tab[x[1]:x[2],1])}
  fUI.Se = function (x){sum(tab[wr:x[2],2])/sum(tab[x[1]:x[2],2])}
  fUI.Acc = function(x) {
    (ifelse(x[1] < wr, sum(tab[x[1]:(wr-1),1]), 0)+sum(tab[wr:x[2],2]))/
      sum(tab[(x[1]:x[2]),1]+tab[(x[1]:x[2]),2]) }
  fUI.ratio = function (x){(sum(tab[x[1]:x[2],2])/totpos)/(sum(tab[x[1]:x[2],1])/totneg)}
  # fUI.C = function (x){sel = test>=d$test[x[1]] & test<=d$test[x[2]] # cbind(test,sel)
  #                         emp.AUC(test[sel & ref==0], test[sel & ref==1])} # ref=ref[sel]; test=test[sel]
  fUI.C2 = function (x){sel = test>=d$test[x[1]] & test<=d$test[x[2]] # cbind(test,sel)
      AUC(test[sel & ref==0], test[sel & ref==1])} # ref=ref[sel]; test=test[sel]
  f.UI.n = function(x){ sum(d$tot[x[1]:x[2]]) }
  ui$UI.Sp = apply(ui, 1, fUI.Sp) # x = unlist(ui[1,])
  ui$UI.Se = apply(ui, 1, fUI.Se)
  ui$UI.Acc = apply(ui[,1:2], 1, fUI.Acc)
  ui$UI.ratio = apply(ui[,1:2], 1, fUI.ratio) # x=c(3,3); fUI.ratio(x)
  # ui$UI.C2 = apply(ui[,1:2], 1, fUI.C)
  ui$UI.C = apply(ui[,1:2], 1, fUI.C2)
  ui$UI.n = apply(ui[,1:2], 1, f.UI.n)

  f.Sp = function(x) {
    sum(d$d0[-(x[1]:dl)])/sum(d$d0[-(x[1]:x[2])]) }
  f.Se = function(x) {
    sum(d$d1[-(1:x[2])])/sum(d$d1[-(x[1]:x[2])]) }
  f.Acc = function(x) {
    (sum(d$d0[-(x[1]:dl)])+sum(d$d1[-(1:x[2])]))/
      sum(d$d0[-(x[1]:x[2])]+d$d1[-(x[1]:x[2])]) }
  # f.C = function (x){sel = test>=d$test[x[1]] & test<=d$test[x[2]] # cbind(test,sel)
  #    emp.AUC(test[!sel & ref==0], test[!sel & ref==1])} # ref=ref[sel]; test=test[sel]
  f.C2 = function (x){sel = test>=d$test[x[1]] & test<=d$test[x[2]] # cbind(test,sel)
    AUC(test[!sel & ref==0], test[!sel & ref==1])} # ref=ref[sel]; test=test[sel]

  f.Loss = function(x) {
    weights[1]*sum(t[-(x[1]:lr), '1']) + # FNR
      weights[2]*sum(abs(t[x[1]:x[2], '0'] - t[x[1]:x[2], '1'])) + # errors for the uncertain interval
        weights[3]*sum(t[-(1:x[2]), '0'])  # FPR
  }
  f.MCI.n = function(x){ sum(d$tot[-(x[1]:x[2])]) } # x=c(4,4); f.Loss(x)

  ui$MCI.Sp = apply(ui, 1, f.Sp)
  ui$MCI.Se = apply(ui, 1, f.Se)
  ui$MCI.Acc = apply(ui[,1:2], 1, f.Acc)
  # ui$MCI.C = apply(ui[,1:2], 1, f.C)
  ui$MCI.C = apply(ui[,1:2], 1, f.C2)
  ui$MCI.n = apply(ui[,1:2], 1, f.MCI.n)
  ui$Loss = apply(ui[,1:2], 1, f.Loss)

  # ui$Random.loss= Acc.Y - ((round(ui$MCI.Acc*ui$MCI.n)+
  #                                p1*ui$UI.n)/ (totpos+totneg))

  # replace row numbers by test values
  ui$lowerbound =  d$test[ui$lowerbound]
  ui$upperbound =  d$test[ui$upperbound]


  # sel = ui$Loss <= Loss.Y
  # ui2 = ui[sel,]
  # ui2 = ui2[order(ui2$Loss),]
  C.lim = ui$UI.C[find.closest(ui$UI.C, constraints['C'])][1]
  Acc.lim = ui$UI.Acc[find.closest(ui$UI.Acc, constraints['Acc'])][1]

  sel = ui$UI.C <= C.lim &  ui$UI.Acc <= Acc.lim
  ui2 = ui[sel,]

  sel = ui2$UI.ratio >= constraints['lower.ratio'] &
    ui2$UI.ratio <= constraints['upper.ratio']

  if (sum(sel) >= 1) {ui2 = ui2[sel,]} else warning('No solution found for the ui ratio constraints. \n')
  # select at least one

  if ('All' %in% select.max)
    ui3 = ui2[order(ui2$Loss),]
  else {
    if ('MCI.Sp+MCI.Se' %in% select.max){
      sel = ui$Loss <= Loss.Y
      if (sum(sel) == 0) warning('No solution found with less loss than Youden threshold. \n')
      # select at least one

      ui3 = ui2[which.max(unlist(ui2['MCI.Se']+ui2['MCI.Sp'])),]
      rownames(ui3) <- c('MCI.Sp+MCI.Se')
      select.max = select.max[!select.max %in% 'MCI.Sp+MCI.Se']
    } else { ui3 = ui2[FALSE,] } # empty data.frame
    i=1
    while (select.max[i] %in% c('MCI.C', 'MCI.Se', 'MCI.Sp', 'MCI.Acc', 'MCI.n'))
    {
      temp = nrow(ui3)
      ui3 = rbind(ui3,ui2[which.max(unlist(ui2[select.max[i]])),])
      rownames(ui3)[temp+1] <- select.max[i]
      i=i+1
    }
  }

  # ui2 ; ui3

  if (return.all) { return(
    list(Youden=c(max.Youden=d$Y[Yt], # max Youden index
    threshold=d$test[Yt], # cut-point test values
    Sp=1-d$fpr[Yt], Se=d$tpr[Yt],Acc=Acc.Y, Loss=Loss.Y, C = Conc),
    data.table = d,
    intersection=isr,
    uncertain.interval=ui,
    candidates=ui3) )
  } else   return(
    list(Youden=c(max.Youden=d$Y[Yt], # max Youden index
                  threshold=d$test[Yt], # cut-point test values
                  sp=1-d$fpr[Yt], se=d$tpr[Yt],acc=Acc.Y, loss=Loss.Y, C = Conc),
         intersection=isr,
         candidates=ui3)
  )
}


