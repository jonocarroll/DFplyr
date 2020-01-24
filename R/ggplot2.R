#' @export
ggplot2::fortify
#' @export
fortify.src_DF <- function(model, data, ...) {
  convert_with_group(model)
}
