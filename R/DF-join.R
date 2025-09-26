check_common_columns <- function(names_x, names_y) {
    by <- intersect(names_x, names_y)

    # from dplyr:::join_by_common
    if (length(by) == 0) {
        rlang::abort(
            "`by` must be supplied when `x` and `y` have no common variables."
        )
    }

    message("Joining with `by = ", deparse(by), "`")

    return(by)
}

is_leftish <- function(...) {
    # does this look like a non-right join?
    args <- list(...)
    if (utils::hasName(args, "all.y") && args$all.y) {
        return(FALSE)
    }
    TRUE
}

join_internal <- function(x, y, by = NULL, ...) {
    use_rownames <- is_leftish(...)
    if (use_rownames) x$...rownames <- rownames(x)
    if (is.null(by)) by <- check_common_columns(names(x), names(y))

    grps <- group_vars(x)

    x <- S4Vectors::merge(x, y, by = by, sort = FALSE, ...)
    if (use_rownames) rownames(x) <- x$...rownames
    if (use_rownames) x$...rownames <- NULL

    grps <- intersect(grps, colnames(x)) # <-- preserve remaining groups
    # if no grouping return object
    if (length(grps) == 0) {
        return(x)
    }
    # else rebuild groups with new DF
    group_by(x, !!!rlang::syms(grps))
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
#' @keywords internal
inner_join.DataFrame <- function(
        x,
        y,
        by = NULL,
        copy = FALSE,
        suffix = c(".x", ".y"),
        ...,
        keep = NULL) {
    join_internal(x = x, y = y, by = by, all = FALSE)
}

#' @export
#' @keywords internal
left_join.DataFrame <- function(
        x,
        y,
        by = NULL,
        copy = FALSE,
        suffix = c(".x", ".y"),
        ...,
        keep = NULL) {
    join_internal(x = x, y = y, by = by, all.x = TRUE)
}

#' @export
#' @keywords internal
right_join.DataFrame <- function(
        x,
        y,
        by = NULL,
        copy = FALSE,
        suffix = c(".x", ".y"),
        ...,
        keep = NULL) {
    join_internal(x = x, y = y, by = by, all.y = TRUE)
}

#' @export
#' @keywords internal
full_join.DataFrame <- function(
        x,
        y,
        by = NULL,
        copy = FALSE,
        suffix = c(".x", ".y"),
        ...,
        keep = NULL) {
    join_internal(x = x, y = y, by = by, all = TRUE)
}
