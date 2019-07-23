#' Function for the determination of a grey zone for quantitative diagnostic and
#' screening tests
#'
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the reference or gold standard. The reference standard must be
#'   coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition).
#' @param test The ordinal test scores under evaluation. Higher scores indicate
#'   the presence of the targeted disease. Please use negated values when lower
#'   values indicate the presence of the targeted disease.
#' @param prevalence The prevalence or pre-test probability to be used. When
#'   NULL, the prevalence found in the sample is used.
#' @param criterion.values The minimum desired values for respectively the
#'   positive and negative post-test probability.
#' @param return.all Default = FALSE. When TRUE the full table of all results
#'   are returned.
#' @details This function is proposed by Coste et al. (2003). The current
#'   implementation only handles ordinal test values. This functions uses all
#'   possible test scores as dichotomous thresholds to calculate Se, Sp,
#'   positive and negative likelihood ratios and post-test probabilities. The
#'   likelihood ratios are calculated for the cumulated densities of the test
#'   scores and indicate the levels of seriousness of the disease for all possible
#'   dichotomous thresholds. It uses therefore a cumulative interpretation of the
#'   Likelihood Ratios and posttest probabilities. If a test has test scores 1
#'   to 5 (with 5 indicating the largest probability of the disease), Se,
#'   positive LR and positive posttest probabilities of the greyzone function
#'   concern test results >= 1, >= 2, >= 3, >= 4 and >= 5, while Sp, negative LR
#'   and negative posttest probabilities concern test results < 1, < 2, < 3, < 4
#'   and < 5.
#'
#'   Please note that the definition of a grey zone deviates from the definition
#'   of an uncertain interval.
#'
#'   The criterion is a required degree of closeness of post-test probabilities
#'   to 1 or 0. These post-test probabilities of cumulated test scores may
#'   require a value over 0.99 or even 0.999 (or under 0.01 or 0.001) to confirm
#'   or exclude the presence of a target disease. The default criterion values
#'   are .05 and .95 for respectively a negative and positive classification,
#'   which may be sufficient for use by clinicians or Public Health
#'   professionals for a first classification whether a target disease may be
#'   present or not (Coste et al., 2003).
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
#'   further discussion of the interval interprestation.
#'
#' @return The function returns the lower and upper value of the range of test
#'   scores that are considered 'grey' or inconclusive. When return.all = TRUE
#'   the full table of the results is returned.
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
#'  test = c(rep(1:5, c(90,75,50,35,0)), c(rep(1:5, c(10,25,50,65,100))))
#'  table(ref, test)
#'  greyzone(ref, test, ret=TRUE)

# ref=irondef$ref; test=irondef$test; prevalence=.1; criterion.values=c(.7, .1); return.all=TRUE
greyzone <- function(ref, test, prevalence=NULL, criterion.values=c(.05, .95), return.all=F){

  find.closest <- function(M, crit){
    M[which(is.na(M))] = Inf # crit+10e10
    mindiff=min(abs(M-crit)) # warning: no non-missing arguments to min; returning Inf
    which((M == crit+mindiff) | (M == crit-mindiff), arr.ind=T)
  }

  df=check.data(ref, test, model='ordinal')
  if (is.null(prevalence)) prevalence = sum(df$ref==1)/length(df$ref)
  tab=as.matrix(table( df$test, df$ref))
  # print(addmargins(tab))
  totpos=sum(tab[,2])                          # The total number of positives (one number)
  totneg=sum(tab[,1])                          # The total number of negatives (one number)

  rdf=data.frame(thresholds=unique(sort(test)), d0=tab[,'0'], d1=tab[,'1'], row.names=1:nrow(tab))
  rdf$tot=rowSums(tab)                             # Number of patients w/ each threshold result
  rdf$TP=unname(rev(cumsum(rev(tab[,2]))))     # Number of true positives, start with maximum
  rdf$FP=unname(rev(cumsum(rev(tab[,1]))))     # Number of false positives, start with maximum
  # rdf$TN=unname(cumsum(tab[,1]))
  rdf$TN=totneg-rdf$FP
  rdf$FN=totpos-rdf$TP

  rdf$tpr=rdf$TP/totpos          # Sensitivity (fraction true positives)
  rdf$fpr=rdf$FP/totneg  # 1-rdf$FP/totneg         # 1 - specificity (false positives)

  rdf$preodds=prevalence/(1-prevalence); # pretest odds
  rdf$sp=1-rdf$fpr;
  rdf$se=rdf$tpr;
  rdf$neg.lr=(1-rdf$se)/rdf$sp; # negative likelihood ratio
  rdf$pos.lr=rdf$se/(1-rdf$sp); # positive likelihood ratio
  # By Bayes' rule : posttest(posterior)odds=likelihoodratio * pretest(prior)odds
  rdf$negpostodds=rdf$neg.lr * rdf$preodds ; # negative post odds
  rdf$pospostodds=rdf$pos.lr * rdf$preodds ; # poitive post odds
  rdf$negpostprob = rdf$negpostodds / (rdf$negpostodds+1); # negative post probability
  rdf$pospostprob=rdf$pospostodds / (rdf$pospostodds+1); # positive post probability

  # rdf$negpostodds[1:2] / (rdf$negpostodds[1:2]+1)
  # rdf[1:2,]
  # ltpos = rdf$negpostprob < criterion.values[1]
  # utpos = rdf$pospostprob > criterion.values[2]

  ltpos = find.closest(rdf$negpostprob, criterion.values[1])
  utpos = find.closest(rdf$pospostprob, criterion.values[2])

  rdf$threshold[ltpos]
  rdf$threshold[utpos]
  rdf$threshold[which(rdf$negpostprob >= criterion.values[1])]
  rdf$threshold[which(rdf$pospostprob >= criterion.values[2])]

  if(length(ltpos) > 1 | length(utpos) > 1) warning('Multiple thresholds found.')

  if (return.all) {
    return(list(table = rdf, thresholds = c(lt = as.numeric(rownames(tab))[ltpos], ut = as.numeric(rownames(tab))[utpos])))
  } else return(c(lt = as.numeric(rownames(tab))[ltpos], ut = as.numeric(rownames(tab))[utpos]) )
}


