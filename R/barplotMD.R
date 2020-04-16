#' Barplot of frequencies, densities or both of the two distributions of
#' patients with and without the targeted condition.
#'
#' @param ref Vector of patient status with two ordered values. The first
#'   indicates the patients without the targeted condition (for instance 0), the
#'   second indicates the patients with the targeted condition (for instance 1).
#'   This order is relevant.
#' @param test Vector of ordinal test scores. The values range from min(test) to
#'   max(test) and need to be ordered. Missing values in between are shown in
#'   the plots as gaps.
#' @param name.test Name used for the title of the x axis.
#' @param fixed.range Default = NULL. If test has numeric values you can set the
#'   values to cover a fixed range. This may enable the comparison of different
#'   samples with truncated test values. Default: use \code{min(test):max(test)}
#'   as the range of values.
#' @param plot Default: 'frequencies'. Which plots to create: 'frequencies',
#'   'densities' or 'both'.
#' @param target.condition Default: 'Target Condition'. Name of target
#'   condition.
#' @param position.legend Default: 'top'. Position of the legend. Most used
#'   values: "topleft", "top", "topright".
#' @param cex.legend Default: 1. Relative size of the legend.
#'
#' @return named matrix with 2 rows and max(test)-min(test)+1 columns that
#'   provide for the position on the x-axis of each of the test values. These
#'   values can be used to draw vertical lines to indicate cutoff scores (see
#'   example).
#' @importFrom graphics barplot
#' @seealso \code{\link{plotMD}}
#' @export
#'
#' @examples
#' data(tostbegg2)
#' barplotMD(ref=tostbegg2$d, test=tostbegg2$y, name='Metastatic Rating', cex=1)
#' x.axis = barplotMD(ref=tostbegg2$d, test=tostbegg2$y, plot='densities',
#' name='Metastatic Rating', cex=1)
#' # Use x.axis to plot vertical line between test score 3 and 4
#' segments(x0=(x.axis[2,4]+x.axis[1,3])/2, y0=0, y1=.4, col='red')
#' # include zero score (in this sample empty)
#' barplotMD(ref=tostbegg2$d, test=tostbegg2$y, fixed.range = c(0, 5),
#'           plot='densities',name='Metastatic Rating', cex=1)
#' op = par(mfrow=c(2,1))
#' barplotMD(ref=tostbegg2$d, test=tostbegg2$y, plot='both',
#'           name='Metastatic Rating', cex.legend=.6, pos='top')
#' par(mfrow=op)

# library(UncertainInterval)
# test = factor(1:3, order=T)
# data(tostbegg2)
# ref=tostbegg2$d; test=tostbegg2$y; name='Metastatic Rating'; cex=1
# plot='frequencies'; target.condition="Target Condition"; position.legend='top'
# cex.legend=1
barplotMD <- function(ref, test, name.test='', fixed.range= c(NULL, NULL),
                      plot=c('frequencies', 'densities', 'both', 'none'),
                      target.condition="Target Condition", position.legend='top',
                      cex.legend=1) {

  df=check.data(ref, test, model='ordinal')
  ref = df$ref
  test = df$test
  plot = match.arg(plot)

  # calculate horizontal postions (x) of barplots
  if (is.numeric(test)) {
    levs = min(test):max(test)
    if (!is.null(fixed.range)) {
      if (is.numeric(fixed.range) & length(fixed.range) == 2) {
        levs = fixed.range[1]:fixed.range[2]
      } else {
        stop("Parameter fixed.range has invalid values.")
      }
    }
    test = factor(test, ordered = TRUE, levels = levs) # is.ordered(test)
  }
  if (is.numeric(ref)) ref = factor(ref, ordered=T, levels = min(ref):max(ref))

  oa.freq = table(ref, test)
  oa.dens = prop.table(oa.freq, margin=1)
  width = 1
  space = c(0, 1)
  space <- space * mean(width)
  NR <- 2
  NC <- ncol(oa.freq) # nr of test scores
  space <- rep.int(c(space[2L], rep.int(space[1L], NR - 1)), NC)
  width <- rep_len(width, nrow(oa.dens))
  delta <- width / 2
  w.r <- cumsum(space + width)
  w.m <- w.r - delta
  w.l <- w.m - delta
  wm = matrix(w.m, ncol = ncol(oa.dens))

  if (plot %in% c('frequencies', 'both')) {
    barplot(
    oa.freq,
    main = "Mixed Frequencies",
    xlab = name.test,
    col = c("white", "black"),
    beside = T
  )
  }
  if (plot %in% c('densities', 'both')) {
    barplot(
    oa.dens,
    main = "Mixed Densities",
    xlab = name.test,
    col = c("white", "black"),
    beside = T
  )
  }
  legend(x=position.legend,
         c(paste(levels(ref)[1], target.condition, "Absent"),
           paste(levels(ref)[2], target.condition, "Present")),
         fill = c("white", "black"), cex=cex.legend)

  invisible(wm)
}

