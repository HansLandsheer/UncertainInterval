#' Run example of the TG-ROC method, simulating a large variety
#' of tests with varying qualities.
#' 
#' @description Documentation is provided in the ShinyApp itself.
#' @seealso \code{\link{runExample.UI}} \code{\link{runExample.GreyZone}} 
#' \code{\link{runExample.ROC}}

#' @export
#' @importFrom shiny runApp
runExample.TGROC <- function() {
  appDir <- system.file("shiny-examples", "TG-ROC", package = "UncertainInterval")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `UncertainInterval` package.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}