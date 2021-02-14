#' Run example of the ROC method, simulating a large variety
#' of tests with varying qualities.
#' 
#' @description Documentation is provided in the ShinyApp itself.
#' @seealso \code{\link{runExample.UI}} \code{\link{runExample.TGROC}} 
#' \code{\link{runExample.GreyZone}} 
#' 

#' @export
#' @importFrom shiny runApp
runExample.ROC <- function() {
  appDir <- system.file("shiny-examples", "ROC", package = "UncertainInterval")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `UncertainInterval` package.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}