#' @export
ggplot2::fortify
#' @export
fortify.DataFrame <- function(model, data, ...) {
  convert_with_group(model)
}
