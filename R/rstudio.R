#' Show Clean Description in RStudio
#' @export
#' @noRd
.rs.describeObject <- function(env, objName) {
    existing <- get(".rs.describeObject",
        envir = as.environment("tools:rstudio")
    )(objName = objName, env = env)
    if (inherits(get(objName), "DataFrame")) {
        savedattr <- attr(existing$value, "class")
        existing$value <- "Formal class DataFrame (dplyr-compatible)"
        attr(existing$value, "class") <- savedattr
    }
    existing
}
