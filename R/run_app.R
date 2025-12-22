#' @title Run sfi Shiny App
#' @description A function to run the shiny app for sfi package
#' @export
#' @importFrom rmarkdown run
#' @return A shiny app
run_app <- function() {
  app_dir <- system.file("shiny", "app.Rmd", package = "sfi")
  if (app_dir == "") {
    stop("Could not find shiny app directory. Try re-installing `sfi`.", call. = FALSE)
  }
  options(shiny.maxRequestSize = 100000 * 1024^2) # Set max upload size to 100GB
  rmarkdown::run(app_dir)
}
