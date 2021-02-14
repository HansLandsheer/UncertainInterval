#' Primitive non-parametric function for calculating false and true positive
#' rate for two comparable samples with ordinal data.
#'
#' @param norm Ordinal data of the norm group (controls).
#' @param abnorm Ordinal data of the abnorm group (patients).
#'
#' @return  List of \itemize{
#' \item $testscores sorted available unique test scores.
#' \item $dich.thresholds Indicates the thresholds and
#'   their interpretation. If \code{mean(norm) > mean(abnorm)}, the thresholds
#'   are <= test score, otherwise the thresholds are >= test score. 
#'   \item $d0
#'   Frequencies of ordinal scores of norm group, from lowest to highest score
#'   \item $d1 Frequencies of ordinal scores of abnorm group, from lowest to
#'   highest score 
#'   \item $TP Cumulative true positive scores. If
#'   \code{mean(norm) > mean(abnorm)} the highest score is the total sum,
#'   otherwise the lowest score is the total sum. 
#'   \item $FP Cumulative false
#'   positive scores. If \code{mean(norm) > mean(abnorm)} the highest score is
#'   the total sum, otherwise the lowest score is the total sum. 
#'   \item $tpr True
#'   positive rates (Sensitivities) for each threshold 
#'   \item $fpr False positive
#'   rates (1 - Specificities) for each threshold }
#' @details This function does not check anything. Argument \code{norm} and
#'   \code{abnorm} must have the correct ordinal data. The thresholds only
#'   concern available test scores and are always ordered from lowest to
#'   highest.
#' @export
#'
#' @examples
#' norm = round(rnorm(100, 3, 1))
#' abnorm= round(rnorm(80, 5, 2))
#' (res=simple_roc3(norm, abnorm))
#' # Plot ROC curve
#' plot(x=res$fpr, y=res$tpr, type='l')
#' abline(a=c(0,0), b=c(1,1))

simple_roc3 <- function(norm, abnorm){
  n0 = length(norm)
  n1 = length(abnorm)
  tab = table(c(norm,abnorm), c(rep(0,n0),rep(1,n1))) 
  if (mean(norm) > mean(abnorm)){
    # highest score is total sum
    TP=unname(cumsum(tab[,2]))    # cumulative sums of abnorm scores
    FP=unname(cumsum(tab[,1]))    # cumsums of norm scores
    thresholds=paste('<=', rownames(tab))
  } else {
    # lowest score is total sum
    TP=unname(rev(cumsum(rev(tab[,2]))))    # cumulative sums of true 1
    FP=unname(rev(cumsum(rev(tab[,1]))))    # cumsums of false 1
    thresholds=paste('>=', rownames(tab))
  }
  data.frame(testscores=as.numeric(row.names(tab)), 
             dich.thresholds=thresholds, d0 = tab[,1], d1=tab[,2], 
                   TP, FP, tpr=TP/n1, fpr=FP/n0, row.names=1:nrow(tab))
}
