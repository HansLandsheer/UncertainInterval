library(shiny)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  output$distPlot <- renderPlot({
    # generate pdf's based on input$bins from ui.R
    FNR <- input$FNR
    if (input$combineSliders) updateSliderInput(session, "FPR", value = FNR)
    FNR <- input$FNR
    FPR = input$FPR / 100
    FNR = input$FNR / 100
    acc = ifelse(input$acc==1, .9, .95)
    
    m0=0; sd0=1; 
    is = qnorm(1-FPR, 0, 1) # x value for 15% FP = intersection
    Z = qnorm(FNR) # Z value = 15% FN; Z = (x-mean) / sd
    dens = dnorm(is, 0, 1) # density at the point of intersection
    sd1 = (1/(dens*sqrt(2*pi)))*exp(-.5 * Z^2)
    m1 = is-Z*sd1 
    
    x <- seq(-4, 4+m1, length=200)
    y0 <- dnorm(x, m0, sd0)
    y1 = dnorm(x, m1, sd1)
    par(mfrow=c(1,2))
    plot(x, y0, type="l", col='black', xlab='Test score', ylab='Density', 
         main='Probability Density Functions', ylim=c(0,.5))
    lines(x, y1, type="l", col='red')
    v1 = qnorm(acc)
    v2 = qnorm(p=1-acc, mean=m1, sd=sd1)
    lines(x=c(v1,v1), y=c(0,.5),col='black', lty=2)
    lines(x=c(v2,v2), y=c(0,.5),col='red', lty=2)
    lines(x=c(is,is), y=c(0,.5),col='blue', lty=3)
    a = (m1-m0)/sd1
    b = sd0/sd1
    AUC = round(unname(pnorm(a/sqrt(1+b^2))), 4)
    legend('topleft', c(paste('Non-Patients = N(', m0,', ', sd0,')', sep=''),
                        paste('True Patients = N(', round(m1,2),', ', round(sd1,2),')', sep=''),
                        paste('AUC = ', AUC, sep='')),
           text.col= c('black','red', 'black'))
    
    legend('topright', c(paste('Se =', acc),paste('Sp =',acc), 'intersection'), lty=c(2,2,3), col=c('red', 'black', 'blue'))
    
    plot(x, 1-pnorm(x, m1, sd1), type='l', col='red', xlab='Test score', 
         ylab='Probability', main="TG-ROC curves of Sensitivity and Specificity") # Sensitivity
    lines(x, pnorm(x, m0, sd0), type='l', col= 'black') # Specificity
    lines(x=c(v1,v1), y=c(0,1),col='black', lty=2)
    lines(x=c(v2,v2), y=c(0,1),col='red', lty=2)
    abline(h=acc)
    if (v1 >= v2){
      Se.VR = (1-pnorm(v1, m1, sd1))/(pnorm(v2, m1, sd1)+ 1-pnorm(v1, m1, sd1))
      Sp.VR = (pnorm(v2, m0, sd0))/(pnorm(v2, m0, sd0)+ 1-pnorm(v1, m0, sd0))
    } else {
      Se.VR = (1-pnorm(v2, m1, sd1))/(pnorm(v1, m1, sd1)+ 1-pnorm(v2, m1, sd1))
      Sp.VR = (pnorm(v1, m0, sd0))/(pnorm(v1, 0, 1)+ 1-pnorm(v2, m0, sd0))
    }
    legend('left', paste('Se.VR =', round(Se.VR,2)), text.col='red')
    legend('right', paste('Sp.VR =', round(Sp.VR,2)), text.col='black')
    par(mfrow=c(1,1))
  })
}
