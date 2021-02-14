#' Run example of the UI (UncertainInterval method), simulating a large variety
#' of tests with varying qualities.
#' 
#' @description Documentation is provided in the ShinyApp itself.
#' @seealso \code{\link{runExample.TGROC}} \code{\link{runExample.GreyZone}} 
#' \code{\link{runExample.ROC}}

#' @export
#' @importFrom shiny runApp
runExample.UI <- function() {
  appDir <- system.file("shiny-examples", "UI", package = "UncertainInterval")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `UncertainInterval` package.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}