#' Show Clean Description in RStudio
#' @export
#' @noRd
.rs.describeObject <- function(env, objName, computeSize = TRUE) {
  existing <- get(".rs.describeObject",
                  envir = as.environment('tools:rstudio'))(objName = objName,
                                                           env = env,
                                                           computeSize = computeSize)
  if (inherits(get(objName), "DataFrame")) {
    savedattr <- attr(existing$value, "class")
    existing$value <- "Formal class DataFrame (dplyr-compatible)"
    attr(existing$value, "class") <- savedattr
  }
  existing
}
