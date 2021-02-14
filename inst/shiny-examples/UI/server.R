library(shiny)
library(UncertainInterval)


server <- function(input, output, session) {
  
  output$distPlot <- renderPlot({
    # generate pdf's based on input$bins from ui.R
    
    if (input$combineSliders) updateSliderInput(session, "FPR", value = input$FNR)
    
    FPR = input$FPR / 100
    FNR = input$FNR / 100
    acc = ifelse(input$acc==1, .55, .60)
    
    m0=0; sd0=1; 
    is = qnorm(1-FPR, 0, 1) # x value for 15% FP = intersection
    Z = qnorm(FNR) # Z value = 15% FN; Z = (x-mean) / sd
    dens = dnorm(is, 0, 1) # density at the point of intersection
    sd1 = (1/(dens*sqrt(2*pi)))*exp(-.5 * Z^2)
    m1 = is-Z*sd1 
    
    x <- seq(-4, 4+m1, length=200)
    y0 <- dnorm(x, m0, sd0)
    y1 = dnorm(x, m1, sd1)
    
    plot(x, y0, type="l", col='black', xlab='Test score', ylab='Density', 
         main='Test Score Distributions and Uncertain Interval (UI)', 
         ylim=c(0,.5)) 
    lines(x, y1, type="l", col='red')
    
    res = nlopt.ui(
      UI.Se = acc,
      UI.Sp = acc,
      mu0 = m0,
      sd0 = sd0,
      mu1 = m1,
      sd1 = sd1,
      intersection = is,
      start = NULL,
      print.level = 0
    )
    v1 = res$solution[1]
    v2 = res$solution[2]
    lines(x=c(v1,v1), y=c(0,.5),col='black', lty=2)
    lines(x=c(v2,v2), y=c(0,.5),col='black', lty=2)
    lines(x=c(is,is), y=c(0,.5),col='blue', lty=3)
    a = (m1-m0)/sd1
    b = sd0/sd1
    AUC = round(unname(pnorm(a/sqrt(1+b^2))), 4)
    p = pnorm(v2, m1, sd1) - pnorm(v1, m1, sd1)
    q = pnorm(v2, m0, sd0) - pnorm(v1, m0, sd0)
    err.ui = pnorm(is, m1, sd1) - pnorm(v1, m1, sd1) +
      pnorm(v2, m1, sd1) - pnorm(is, m1, sd1)
    legend('topleft', c(paste(' Non-Patients = N(', m0,', ', sd0,')', sep=''),
                        paste('True Patients = N(', round(m1,2),', ', round(sd1,2),')', sep=''),
                        paste('AUC = ', AUC, sep=''),
                        paste('UI.size = ', round((p+q)/2 * 100,0), '%', sep=''),
                        paste('UI.err = ', round((err.ui/(FNR+FPR))*100, 0), '%', sep='')),
           text.col= c('black','red', 'black', 'black'))
    
    legend('topright', c(paste('Se.UI =', round(res$results["exp.UI.Se"], 2)),
                         paste('Sp.Ui =',round(res$results["exp.UI.Sp"], 2)), 'intersection',
                         "borders UI"), lty=c(0, 0,3, 2), col=c('red', 'black', 'blue', 'black'), 
           text.col=c('red', 'black', 'blue', 'black'))
    
    Se.opt = 1-pnorm(is, m1, sd1)
    Sp.opt = pnorm(is, m0, sd0)
    Se.MCI = (1-pnorm(v2, m1, sd1))/(pnorm(v1, m1, sd1)+ 1-pnorm(v2, m1, sd1))
    Sp.MCI = (pnorm(v1, m0, sd0))/(pnorm(v1, m0, sd0)+ 1-pnorm(v2, m0, sd0))
    
    legend('right', c(paste('Se.opt =', round(Se.opt,2)),
                      paste('Se.MCI =', round(Se.MCI,2))), text.col=c('red', 'red'))
    legend('left', c(paste('Sp.opt =', round(Sp.opt,2)),
                     paste('Sp.MCI =', round(Sp.MCI,2))), text.col=c('black', 'black'))
  })
  
  output$overlap = renderText(paste("Overlap Distributions =", input$FNR+input$FPR, '%'))
  
}
