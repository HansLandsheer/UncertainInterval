pkgname <- "UncertainInterval"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('UncertainInterval')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("RPV")
### * RPV

flush(stderr()); flush(stdout())

### Name: RPV
### Title: Trichotomization of ordinal test results using predictive values
### Aliases: RPV

### ** Examples

set.seed(1)
# example of a validation sample
ref=c(rep(0,1000), rep(1, 1000))
test=round(c(rnorm(1000, 5, 1), rnorm(1000, 8, 2)))
# calculated roll.length is invalid. Set to 3. Post test probability equals
# Positive Predictive Values. Parameter pretest.prob is set to sample prevalence.
RPV(ref, test, reliability = .9, roll.length = 3)
# Set roll.length = 1 to ignore test reliability (value of parameter
# reliability is ignored, but must be set to some value.)
RPV(ref, test, reliability = 0, roll.length = 1)
# When pretest.prob is set to .5, the Post-test Probabilities are equal to
# the Standardized Positive Predictive Values.
RPV(ref, test, pretest.prob = .5, reliability = .9, roll.length = 3)




cleanEx()
nameEx("TG.ROC")
### * TG.ROC

flush(stderr()); flush(stdout())

### Name: TG.ROC
### Title: Two-Graphs Receiving Operating Characteristics.
### Aliases: TG.ROC

### ** Examples

ref = c(rep(0,100), rep(1,100))
test = c(rnorm(100, 0, 1), rnorm(100, 1, 1))
TG.ROC(ref, test, model='binormal', plot=TRUE)
TG.ROC(ref, test, model='none', plot=TRUE)



cleanEx()
nameEx("barplotMD")
### * barplotMD

flush(stderr()); flush(stdout())

### Name: barplotMD
### Title: Barplot of frequencies, densities or both of the two
###   distributions of patients with and without the targeted condition.
### Aliases: barplotMD

### ** Examples

data(tostbegg2)
barplotMD(ref=tostbegg2$d, test=tostbegg2$y, name='Metastatic Rating', cex=1)
x.axis = barplotMD(ref=tostbegg2$d, test=tostbegg2$y, plot='densities',
name='Metastatic Rating', cex=1)
# Use x.axis to plot vertical line between test score 3 and 4
segments(x0=(x.axis[2,4]+x.axis[1,3])/2, y0=0, y1=.4, col='red')
# include zero score (in this sample empty)
barplotMD(ref=tostbegg2$d, test=tostbegg2$y, fixed.range = c(0, 5),
          plot='densities',name='Metastatic Rating', cex=1)
op = par(mfrow=c(2,1))
barplotMD(ref=tostbegg2$d, test=tostbegg2$y, plot='both',
          name='Metastatic Rating', cex.legend=.6, pos='top')
par(mfrow=op)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("check.data")
### * check.data

flush(stderr()); flush(stdout())

### Name: check.data
### Title: Function to check the dataset of individuals with (1) and
###   without (0) the targeted condition.
### Aliases: check.data

### ** Examples

set.seed(1)
ref=c(rep(0,500), rep(1,500))
test=c(rnorm(500,0,1), rnorm(500,1,1.2))
check.data(ref, test) # model = 'kernel'



cleanEx()
nameEx("get.intersection")
### * get.intersection

flush(stderr()); flush(stdout())

### Name: get.intersection
### Title: get.intersection Obtain the intersection of two distributions
###   using the kernel method
### Aliases: get.intersection

### ** Examples

ref=c(rep(0,500), rep(1,500))
test=c(rnorm(500,0,1), rnorm(500,1,2))
(get.intersection(ref, test)) # two intersections! Generates warning in other functions!



cleanEx()
nameEx("greyzone")
### * greyzone

flush(stderr()); flush(stdout())

### Name: greyzone
### Title: Function for the determination of a grey zone for quantitative
###   diagnostic and screening tests
### Aliases: greyzone

### ** Examples

 ref=c(rep(0, 250), rep(1, 250))
 test = c(rep(1:5, c(90,75,50,35,0)), c(rep(1:5, c(10,25,50,65,100))))
 table(ref, test)
 greyzone(ref, test, ret=TRUE)



cleanEx()
nameEx("nlopt.ui")
### * nlopt.ui

flush(stderr()); flush(stdout())

### Name: nlopt.ui
### Title: Function for the determination of the population thresholds an
###   uncertain and inconclusive interval for bi-normal distributed test
###   scores.
### Aliases: nlopt.ui

### ** Examples

# A simple test model:
nlopt.ui()
# Using another bi-normal distribution:
nlopt.ui(mu0=0, sd0=1, mu1=1.6, sd1=2)




cleanEx()
nameEx("nlopt.ui.general")
### * nlopt.ui.general

flush(stderr()); flush(stdout())

### Name: nlopt.ui.general
### Title: Function for the determination of the population thresholds an
###   uncertain and inconclusive interval for test scores with a known
###   common distribution.
### Aliases: nlopt.ui.general

### ** Examples

# A simple test model:
nlopt.ui.general(Se = .55, Sp = .55,
                 distribution = "norm",
                 parameters.d0 = c(mean = 0, sd = 1),
                 parameters.d1 = c(mean = 1, sd = 1),
                 overlap.interval=c(-2,3))
# Standard procedure when using a continuous distribution:
nlopt.ui.general(parameters.d0 = c(mean = 0, sd = 1),
                 parameters.d1 = c(mean = 1.6, sd = 2))
# Function to calculate the Area under the Receiving Operating Characteristics
# Curve (AUC or C-statistic)
emp.AUC <- function(norm, abnorm) {
  o = outer(abnorm, norm, "-")
  mean((o > 0) + .5 * (o == 0))
}

library(MASS)
library(car)
# gamma distributed data
set.seed(4)
d0 = rgamma(100, shape=2, rate=.5)
d1 = rgamma(100, shape=7.5, rate=1)
# 1. obtain parameters
parameters.d0=fitdistr(d0, 'gamma')$estimate
parameters.d1=fitdistr(d1, 'gamma')$estimate
# 2. test if supposed distributions (gamma) is fitting
qqPlot(d0, distribution='gamma', shape=parameters.d0['shape'])
qqPlot(d1, distribution='gamma', shape=parameters.d1['shape'])
# 3. draw curves and determine overlap
curve(dgamma(x, shape=parameters.d0['shape'], rate=parameters.d0['rate']), from=0, to=16)
curve(dgamma(x, shape=parameters.d1['shape'], rate=parameters.d1['rate']), from=0, to=16, add=TRUE)
overlap.interval=c(1, 15) # ignore intersection at 0; observe large overlap
# 4. get empirical AUC
emp.AUC(d0, d1)
# about .65 --> Poor
# .90-1 = excellent (A)
# .80-.90 = good (B)
# .70-.80 = fair (C)
# .60-.70 = poor (D)
# .50-.60 = fail (F)
# 5. Get uncertain interval
(res=nlopt.ui.general (Se = .57,
                       Sp = .57,
                       distribution = 'gamma',
                       parameters.d0 = parameters.d0,
                       parameters.d1 = parameters.d1,
                       overlap.interval,
                       intersection = NULL,
                       start = NULL,
                       print.level = 0))
abline(v=c(res$intersection, res$solution))
# 6. Assess improvement when diagnosing outside the uncertain interval
sel.d0 = d0 < res$solution[1] |  d0 > res$solution[2]
sel.d1 = d1 < res$solution[1] |  d1 > res$solution[2]
(percentage.selected.d0 = sum(sel.d0) / length(d0))
(percentage.selected.d1 = sum(sel.d1) / length(d1))
emp.AUC(d0[sel.d0], d1[sel.d1])
# AUC for selected scores outside the uncertain interval
emp.AUC(d0[!sel.d0], d1[!sel.d1])
# AUC for deselected scores; worst are deselected
# weibull distributed data
set.seed(4)
d0 = rweibull(100, shape=3, scale=50)
d1 = rweibull(100, shape=3, scale=70)
# 1. obtain parameters
parameters.d0=fitdistr(d0, 'weibull')$estimate
parameters.d1=fitdistr(d1, 'weibull')$estimate
# 2. test if supposed distributions (gamma) is fitting
qqPlot(d0, distribution='weibull', shape=parameters.d0['shape'])
qqPlot(d1, distribution='weibull', shape=parameters.d1['shape'])
# 3. draw curves and determine overlap
curve(dweibull(x, shape=parameters.d0['shape'],
      scale=parameters.d0['scale']), from=0, to=150)
curve(dweibull(x, shape=parameters.d1['shape'],
      scale=parameters.d1['scale']), from=0, to=150, add=TRUE)
overlap.interval=c(1, 100) # ignore intersection at 0; observe overlap
# 4. get empirical AUC
emp.AUC(d0, d1)
# about .65 --> Poor
# .90-1 = excellent (A)
# .80-.90 = good (B)
# .70-.80 = fair (C)
# .60-.70 = poor (D)
# .50-.60 = fail (F)
# 5. Get uncertain interval
(res=nlopt.ui.general (Se = .55,
                       Sp = .55,
                       distribution = 'weibull',
                       parameters.d0 = parameters.d0,
                       parameters.d1 = parameters.d1,
                       overlap.interval,
                       intersection = NULL,
                       start = NULL,
                       print.level = 0))
abline(v=c(res$intersection, res$solution))
# 6. Assess improvement when diagnosing outside the uncertain interval
sel.d0 = d0 < res$solution[1] |  d0 > res$solution[2]
sel.d1 = d1 < res$solution[1] |  d1 > res$solution[2]
(percentage.selected.d0 = sum(sel.d0) / length(d0))
(percentage.selected.d1 = sum(sel.d1) / length(d1))
emp.AUC(d0[sel.d0], d1[sel.d1])
# AUC for selected scores outside the uncertain interval
emp.AUC(d0[!sel.d0], d1[!sel.d1])
# AUC for deselected scores; these scores are almost indistinguishable



cleanEx()
nameEx("nomogram")
### * nomogram

flush(stderr()); flush(stdout())

### Name: nomogram
### Title: Fagan's nomogram to show the relationships between the prior
###   probability, the likelihood ratios, sensitivity and specificity, and
###   the posterior probability.
### Aliases: nomogram

### ** Examples

# Show calculated results (first 3 times about the same)
(nomogram(prob.pre.test = .10, probs.post.test=c(pos=.70, neg=.001), plot=FALSE))
(nomogram(prob.pre.test = .10, SeSp=c(Se=0.991416309, Sp=0.952789700), plot=FALSE))
(nomogram(prob.pre.test = .10, LR=c(pos=21, neg=0.0090090091), plot=FALSE))
(nomogram(prob.pre.test = .10, SeSp=c(Se=0.99, Sp=0.95), plot=FALSE))
# plot only
nomogram(prob.pre.test = .10, LR=c(pos=21, neg=0.0090090091))
# plot and display precise results
(nomogram(prob.pre.test = .10, probs.post.test=c(pos=.70, neg=.001)))

# check the influence of different values of prevalence
i=1
out=matrix(0,nrow = 9, ncol= 7)
for (prev in (seq(.1, .9, by=.1))) {
  out[i,]=nomogram(prob.pre.test=prev, probs.post.test=c(.95, .05), plot=FALSE)
  i=i+1
}
colnames(out) = names(nomogram(prob.pre.test=prev, probs.post.test=c(.95, .05), plot=FALSE))
out




cleanEx()
nameEx("plotMD")
### * plotMD

flush(stderr()); flush(stdout())

### Name: plotMD
### Title: Function to plot the mixed densities of distributions of
###   individuals with (1) and without (0) the targeted condition.
### Aliases: plotMD

### ** Examples

# A test of intermediate quality
set.seed(1)
ref=c(rep(0,500), rep(1,500))
test=c(rnorm(500,0,1), rnorm(500,1,1.2))
plotMD(ref, test)
ua = ui.nonpar(ref, test) # with warning message!
# Add lines to indicate Uncertain Interval
abline(v=ua[1:2])
select=(test <= ua[2] & test >= ua[1])
# plot the mixed densities for the Uncertain Interval
plotMD(ref[select], test[select])
plotMD(ref[select], test[select], colspace='gray')
plotMD(ref[select], test[select], colspace='BW')

# An ordinal test
norm     = rep(1:5, times=c(33,6,6,11,2))
abnorm   = rep(1:5, times=c(3,2,2,11,33))
testres  = c(abnorm,norm)
truestat = c(rep(1,length(abnorm)), rep(0,length(norm)))
plotMD(ref=truestat, test=testres, model='ordinal')

# ordinal test: weak test
set.seed(2)
nobs=1000
Z0 <- rnorm(nobs, mean=0)
b0=seq(-5, 5, length.out=31) # range sufficient to cover both z0 and z1
f0=cut(Z0, breaks = b0, labels = c(1:30))
x0=as.numeric(levels(f0))[f0]
Z1 <- rnorm(nobs, mean=.5) # very weak test, not recommended for practical use
f1=cut(Z1, breaks = b0, labels = c(1:30))
x1=as.numeric(levels(f1))[f1]
test=c(x0, x1)
ref =c(rep(0, length(x0)), rep(1, length(x1)))
(pr=prop.table(table(ref, test)))
breaks=c(min(test)-.5, seq(min(test), max(test), by=1)+.5)
plotMD(ref, test, model='ordinal')
# when model = 'binormal' or 'kernel', default breaks do not work well for
# ordinal data, and have to be set by hand
plotMD(ref, test, breaks=c(min(test)-.5, seq(min(test), max(test), by=1)+.5),
       model='binormal')
plotMD(ref, test, breaks=c(min(test)-.5, seq(min(test), max(test), by=1)+.5),
       model='kernel')



cleanEx()
nameEx("quality.threshold")
### * quality.threshold

flush(stderr()); flush(stdout())

### Name: quality.threshold
### Title: Function for the description of the qualities of one or two
###   decision thresholds or threshold.
### Aliases: quality.threshold

### ** Examples

# A simple test
ref=c(rep(0,500), rep(1,500))
test=c(rnorm(500,0,1), rnorm(500,1,1))
ua = ui.nonpar(ref, test)
quality.threshold(ref, test, threshold=ua[1], threshold.upper=ua[2])



cleanEx()
nameEx("quality.threshold.uncertain")
### * quality.threshold.uncertain

flush(stderr()); flush(stdout())

### Name: quality.threshold.uncertain
### Title: Function for the description of the qualities of the Uncertain
###   Interval.
### Aliases: quality.threshold.uncertain

### ** Examples

# A simple test model
ref=c(rep(0,500), rep(1,500))
test=c(rnorm(500,0,1), rnorm(500,1,sd=1))
ua = ui.nonpar(ref, test)
quality.threshold.uncertain(ref, test, ua[1], ua[2])



cleanEx()
nameEx("synthdata_NACC")
### * synthdata_NACC

flush(stderr()); flush(stdout())

### Name: synthdata_NACC
### Title: synthdata NACC
### Aliases: synthdata_NACC
### Keywords: data

### ** Examples

data(synthdata_NACC) # needs R version 3.5 or later
head(synthdata_NACC) # Show head of the dataset
nrow(synthdata_NACC) # total number of observations
# select part of data for the first measurement
# N.B. ref is not available when it is inconclusive
m1 = synthdata_NACC[!is.na(synthdata_NACC$MOCATOTS.1)
                   & !is.na(synthdata_NACC$ref.1), ]
# preliminary check data for possible missing values
addmargins(table(m1$ref.1, m1$MOCATOTS.1, useNA = 'always'))
# Show the data
barplotMD(m1$ref.1, m1$MOCATOTS.1)

# calculate the difference between the two measurements in days
ddiff = (m1$vdate.2 - m1$vdate.1)
# There is a wide variety !!!
summary(ddiff)
# Estimate the test-retest reliability
library(psych)
ICC(na.omit(cbind(m1$MOCATOTS.1, m1$MOCATOTS.2)))
# Reducing the variety of time between measurements:
timesel = (ddiff >= 335) & (ddiff <= 395)
ICC(na.omit(cbind(m1$MOCATOTS.1[timesel], m1$MOCATOTS.2[timesel])))

# error when using default calculated value for roll.length
# RPV(m1$ref.1, m1$MOCATOTS.1, reliability = .86)
RPV(m1$ref.1, m1$MOCATOTS.1, reliability = .86, roll.length = 5)




cleanEx()
nameEx("ui.binormal")
### * ui.binormal

flush(stderr()); flush(stdout())

### Name: ui.binormal
### Title: Function for the determination of the thresholds of an uncertain
###   interval for bi-normal distributed test scores that are considered as
###   inconclusive.
### Aliases: ui.binormal

### ** Examples

# A simple test model
ref=c(rep(0,500), rep(1,500))
test=c(rnorm(500,0,1), rnorm(500,1,1))
ui.binormal(ref, test)




cleanEx()
nameEx("ui.nonpar")
### * ui.nonpar

flush(stderr()); flush(stdout())

### Name: ui.nonpar
### Title: Function for the determination of an inconclusive interval for
###   continuous test scores
### Aliases: ui.nonpar

### ** Examples

# A simple test model
set.seed(1)
ref=c(rep(0,500), rep(1,500))
test=c(rnorm(500,0,1), rnorm(500,1,1))
ui.nonpar(ref, test, select='limited')

ref = c(rep(0,20), rep(1,20))
test= c(rnorm(20), rnorm(20, mean=1))
ui.nonpar(ref, test)




cleanEx()
nameEx("ui.ordinal")
### * ui.ordinal

flush(stderr()); flush(stdout())

### Name: ui.ordinal
### Title: Function to explore possible uncertain intervals of ordinal test
###   results of individuals with (1) and without (0) the targeted
###   condition.
### Aliases: ui.ordinal

### ** Examples

# A short test with 5 ordinal values
test0     = rep(1:5, times=c(165,14,16,55, 10)) # test results norm group
test1     = rep(1:5, times=c( 15,11,13,55,164)) # test results of patients
ref = c(rep(0, length(test0)), rep(1, length(test1)))
test = c(test0, test1)
table(ref, test)
plotMD(ref, test, model='ordinal') # visual inspection
ui.ordinal(ref, test, select.max='All')
# Same solution, but other layout of the results:
ui.ordinal(ref, test, select.max=c('MCI.Sp+MCI.Se', 'MCI.C', 'MCI.Acc',
                                   'MCI.Se', 'MCI.Sp', 'MCI.n'))
# forcing the Youden threshold as intersection gives the same best result.
# However, the estimates for ui.Se, ui.Sp and ui.Acc differ:
ui.ordinal(ref, test, intersection='Youden', select.max='All')

nobs=1000
set.seed(6)
Z0 <- rnorm(nobs, mean=0)
b0=seq(-5, 8, length.out=31)
f0=cut(Z0, breaks = b0, labels = c(1:30))
x0=as.numeric(levels(f0))[f0]
Z1 <- rnorm(nobs, mean=1, sd=1.5)
f1=cut(Z1, breaks = b0, labels = c(1:30))
x1=as.numeric(levels(f1))[f1]
ref=c(rep(0,nobs), rep(1,nobs))
test=c(x0,x1)
plotMD(ref, test, model='ordinal') # looks like binormal
# looks less binormal, but in fact it is a useful approximation:
plotMD(ref, test, model='binormal')
ui.ordinal(ref, test)
ui.binormal(ref, test) # compare application of the bi-normal model



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
