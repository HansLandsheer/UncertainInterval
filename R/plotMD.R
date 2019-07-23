#' Function to plot the mixed densities of distributions of individuals with (1)
#' and without (0) the targeted condition.
#' @name plotMD
#' @description This plot function shows the densities of the two distributions
#'   and their overlap in a single graph.
#' @param ref The reference standard. A column in a data frame or a vector
#'   indicating the classification by the reference test. The reference standard
#'   must be coded either as 0 (absence of the condition) or 1 (presence of the
#'   condition)
#' @param test The index test or test under evaluation. A column in a dataset or
#'   vector indicating the test results in a continuous scale.
#' @param breaks Breaks used to construct the histograms. Either a single
#'   integer number or a vector containing the actual breaks. In the case of a
#'   vector, the number should cover all available test values. In the case of a
#'   single integer number, this number has to be equal or lower than the
#'   discernable values in the test. For short ordinal scales a vector should be
#'   uses covering all possible test values.
#' @param subtitle Optional subtitle
#' @param position.legend The location can be specified by a single keyword from
#'   the list "topright", "topleft", "top", "right", "bottomright", "bottom",
#'   "bottomleft", "left" and "center". Default is "top.right".
#' @param colspace Use colors, grayscale or only black and white as plot colors.
#'   Default = color.
#' @param model The model used for estimation. Default = 'kernel'. Adapts also
#'   breaks and the call to the density function (parameter adjust). When the
#'   model is obviously wrong, warnings are produced.
#' @param ... passing arguments to the kernel density function, other than
#'   kernel='gaussian' (default).
#' @details The graph shows the densities of the two distributions and their
#'   overlap. Many tests of intermediate quality have a considerable overlap.
#'   Also, the distributions as estimated by the \code{density} function, using
#'   the gaussian kernel is shown. The intersection is indicated by a vertical
#'   line. This graph allows the visual inspection of the two distributions, as
#'   well a visual inspection of the approximation of the \code{density}, based
#'   on the gaussian kernel. When the density estimation is way off, the
#'   standard estimation of the intersection will be incorrect, and another
#'   estimation has to be supplied.
#'
#'   The function \code{plotMD} can also be used for visual inspection of the
#'   Uncertain Interval (see examples). Please note that the sensitivity and
#'   specificity values > .5 (including the default of .55) allows for some
#'   positive bias.
#' @seealso \code{\link{barplotMD}}
#' @return No Value returned.
#' @importFrom grDevices rgb col2rgb gray
#' @importFrom graphics abline hist legend lines mtext plot par
#' @importFrom stats addmargins density pchisq t.test
#' @importFrom utils tail
#' @export
#'
#' @references Landsheer, J. A. (2016). Interval of Uncertainty: An Alternative
#'   Approach for the Determination of Decision Thresholds, with an Illustrative
#'   Application for the Prediction of Prostate Cancer. PloS One, 11(11),
#'   e0166007.
#'
#'   Landsheer, J. A. (2018). The Clinical Relevance of Methods for Handling
#'   Inconclusive Medical Test Results: Quantification of Uncertainty in Medical
#'   Decision-Making and Screening. Diagnostics, 8(2), 32.
#'   https://doi.org/10.3390/diagnostics8020032

#' @examples
#' # A test of intermediate quality
#' set.seed(1)
#' ref=c(rep(0,500), rep(1,500))
#' test=c(rnorm(500,0,1), rnorm(500,1,1.2))
#' plotMD(ref, test)
#' ua = ui.nonpar(ref, test) # with warning message!
#' # Add lines to indicate Uncertain Interval
#' abline(v=ua[1:2])
#' select=(test <= ua[2] & test >= ua[1])
#' # plot the mixed densities for the Uncertain Interval
#' plotMD(ref[select], test[select])
#' plotMD(ref[select], test[select], colspace='gray')
#' plotMD(ref[select], test[select], colspace='BW')
#'
#' # An ordinal test
#' norm     = rep(1:5, times=c(33,6,6,11,2))
#' abnorm   = rep(1:5, times=c(3,2,2,11,33))
#' testres  = c(abnorm,norm)
#' truestat = c(rep(1,length(abnorm)), rep(0,length(norm)))
#' plotMD(ref=truestat, test=testres, model='ordinal')
#'
#' # ordinal test: weak test
#' set.seed(2)
#' nobs=1000
#' Z0 <- rnorm(nobs, mean=0)
#' b0=seq(-5, 5, length.out=31) # range sufficient to cover both z0 and z1
#' f0=cut(Z0, breaks = b0, labels = c(1:30))
#' x0=as.numeric(levels(f0))[f0]
#' Z1 <- rnorm(nobs, mean=.5) # very weak test, not recommended for practical use
#' f1=cut(Z1, breaks = b0, labels = c(1:30))
#' x1=as.numeric(levels(f1))[f1]
#' test=c(x0, x1)
#' ref =c(rep(0, length(x0)), rep(1, length(x1)))
#' (pr=prop.table(table(ref, test)))
#' breaks=c(min(test)-.5, seq(min(test), max(test), by=1)+.5)
#' plotMD(ref, test, model='ordinal')
#' # when model = 'binormal' or 'kernel', default breaks do not work well for
#' # ordinal data, and have to be set by hand
#' plotMD(ref, test, breaks=c(min(test)-.5, seq(min(test), max(test), by=1)+.5),
#'        model='binormal')
#' plotMD(ref, test, breaks=c(min(test)-.5, seq(min(test), max(test), by=1)+.5),
#'        model='kernel')


# colspace='color'; model='binormal'; subtitle=''; position.legend='topright'; colspace='color'; model = 'kernel'
plotMD<-function(ref, test, breaks=20, subtitle='',
                 position.legend='topright',
                 colspace=c('color', 'grayscale', 'BW'),
                 model = c('kernel', 'binormal', 'ordinal'), ...){

  colspace <- match.arg(colspace)
  model <- match.arg(model)
  add.alpha <- function(col, alpha=1){
    x = col2rgb(col)/255
    rgb(x[1], x[2], x[3], alpha=alpha)
  }

  df=check.data(ref, test, model=model)
  # breaks for continuous datasets
  if (length(breaks) == 1  & breaks[1]%%1 ==0 ) {
    low = min(df$test)
    high = max(df$test)
    # p1=mean(df$ref) # getAnywhere(hist.default); methods(hist)
    breaks = seq(low, high, (high - low) / min(breaks, length(unique(df$test))))
  }
  # else {
  #   breaks = breaks
  # }
  if (model=='kernel'){
    d0=density(df$test[df$ref==0], ...)
    d1=density(df$test[df$ref==1], ...)
    intersection=get.intersection(df$ref, df$test, model=model, ...) # raw data always used

  } else if (model=='ordinal'){
    breaks = seq(min(test)-.5, max(test)+.5, by=1)
    if (any(names(list(...))== 'adjust')) {
      d0=density(df$test[df$ref==0], ...)
      d1=density(df$test[df$ref==1], ...)
      intersection=get.intersection(df$ref, df$test, model=model, ...)
    } else {
      d0=density(df$test[df$ref==0], adjust=2, ...)
      d1=density(df$test[df$ref==1], adjust=2, ...)
      intersection=get.intersection(df$ref, df$test, model=model, adjust=2, ...)
    }
  } else if (model=='binormal'){
    mu0 = mean(df$test[df$ref == 0])
    sd0 = sd(df$test[df$ref == 0])
    mu1 = mean(df$test[df$ref == 1])
    sd1 = sd(df$test[df$ref == 1])
    xs <- seq(min(mu0 - 3*sd0, mu1 - 3*sd1), max(mu0 + 3*sd0, mu1 + 3*sd1), .01)
    d0 = dnorm(xs, mu0, sd0)
    d1 = dnorm(xs, mu1, sd1)
    intersection=get.intersection(df$ref, df$test, model=model, ...) # raw data always used
  }

  hy0=hist(df$test[df$ref==0], breaks=breaks, plot=FALSE) # str(hy0)
  hy1=hist(df$test[df$ref==1], breaks=breaks, plot=FALSE)


  if (colspace=='BW'){
    overlap = hy0 # str(hy0) hy0$counts/sum(hy0$counts)
    for(i in 1:length(overlap$density)){
      if(hy0$density[i] > 0 & hy1$density[i] > 0){
        overlap$density[i] = min(hy0$density[i],hy1$density[i])
      } else {
        overlap$density[i] = 0
      }
    }
    plot(hy0, freq=FALSE, main='Mixed Densities',
         col='white',border=T, xlab='Test | Predictor',
         xlim=c(min(df$test), max(df$test)),
         ylim = c(0, max(max(hy0$density), max(hy1$density))) )
    plot(hy1, freq=FALSE, xlim=c(min(df$test), max(df$test)),
         add=T,col='black',border=T, density=15, angle=135)
    plot(overlap, freq=F, xlim=c(min(df$test), max(df$test)), col='black', add=T,
         density=15, angle=45)
    plot(overlap, freq=F, xlim=c(min(df$test), max(df$test)), col='black', add=T,
         density=15, angle=135)
    legend(position.legend, horiz=F, col='black', angle=c(0,45,135), density=c(0,15,15),
           legend=c('0','overlap','1'), bg="transparent", box.col='transparent')
    op = par(bg='transparent')
    on.exit(par(op))
    legend(position.legend, horiz=F, col='black', angle=c(0,135,135), density=c(0,15,15),
           legend=c('0','overlap','1'), bg="transparent", box.col='transparent')
    if (model=='kernel'){
      lines(d0, col='black')
      lines(d1, col='black')
    } else {
      lines(xs, d0, col='black')
      lines(xs, d1, col='black')
    }
  } else {
    r1=0; g1=0; b1=1; a1=1/4
    r2=1; g2=0; b2=0; a2=1/4
    a3 = 1-1*(1-a1)*(1-a2)
    col1 = ifelse(colspace=='color',rgb(r1,g1,b1, a1), gray(0.21*r1+0.71*g1+0.07*b1,a1))
    col2 = ifelse(colspace=='color',rgb(r2,g2,b2, a2), gray(0.21*r2+0.71*g2+0.07*b2,a2))

    plot(hy0, freq=FALSE, main='Mixed Densities',
         col=col1,border=FALSE, xlab='Test | Predictor',
         xlim=c(min(df$test), max(df$test)),
         ylim = c(0, max(max(hy0$density), max(hy1$density))) )
    plot(hy1, freq=FALSE, xlim=c(min(df$test), max(df$test)),
         add=T,col=col2,border=FALSE)
    r3 = r1*a1/a3 + r2*a2*(1-a1)/a3
    g3 = g1*a1/a3 + g2*a2*(1-a1)/a3
    b3 = b1*a1/a3 + b2*a2*(1-a1)/a3
    col3 = ifelse(colspace=='color',rgb(r3,g3,b3, a3), gray(0.21*r3+0.71*g3+0.07*b3,a3))
    legend(position.legend, horiz=F, border=0, bg="transparent", box.col='transparent' ,
           legend=c('0','overlap','1'),fill=c(col1, col3, col2),
           cex=.8)
    if (model=='kernel' | model=='ordinal'){
      lines(d0, col=add.alpha(col1))
      lines(d1, col=add.alpha(col2))
    } else {
      lines(xs, d0, col=add.alpha(col1))
      lines(xs, d1, col=add.alpha(col2))
    }
  }

  abline(v=intersection)

  mtext(subtitle)
}
