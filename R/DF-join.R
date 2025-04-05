check_common_columns <- function(names_x, names_y) {
  by <- intersect(names_x, names_y)

  # from dplyr:::join_by_common
  if (length(by) == 0)
    rlang::abort("`by` must be supplied when `x` and `y` have no common variables.")

  message("Joining with `by = ", deparse(by), "`")

  return(by)
}

join_internal <- function(x, y, by = NULL, ...) {
  if (is.null(by))
    by <- check_common_columns(names(x), names(y))

  grps <- get_group_data(x)
  x <- S4Vectors::merge(x, y, by = by, sort = FALSE, ...)
  set_group_data(x, grps)
}



#' @title Mutating joins
#'
#' @param x a `DataFrame`
#' @param y a `DataFrame` or `data.frame`
#' @param by columns to use for joining the objects. If `NULL`, the function will look for common columns.
#'
#' @return a `DataFrame`
#' @examples
#' library(dplyr)
#' library(S4Vectors)
#' da <- starwars[, c("name", "mass", "species")][1:10, ]
#' db <- starwars[, c("name", "homeworld")]
#'
#' Da <- as(da, "DataFrame")
#' Db <- as(db, "DataFrame")
#'
#' Res_inner <- inner_join(Da, Db[1:3, ])
#'
#' @export
inner_join.DataFrame <- function(x, y, by = NULL)
  join_internal(x = x, y = y, by = by, all = FALSE)

#' @export
left_join.DataFrame <- function(x, y, by = NULL)
  join_internal(x = x, y = y, by = by, all.x = TRUE)

#' @export
right_join.DataFrame <- function(x, y, by = NULL)
  join_internal(x = x, y = y, by = by, all.y = TRUE)

#' @export
full_join.DataFrame <- function(x, y, by = NULL)
  join_internal(x = x, y = y, by = by, all = TRUE)
