#' Treat a `S4Vectors::DataFrame` as a `dplyr` data source
#'
#' Add \pkg{dplyr} compatibility to `S4Vectors::DataFrame`
#' for use with a selection of \pkg{dplyr} verbs.
#'
#' @param x A `S4Vectors::DataFrame` object
#'
#' @md
#' @export
#' @examples
#' library(S4Vectors)
#' library(dplyr)
#'
#' DF <- as(mtcars, "DataFrame")
#' d <- add_dplyr_compat(DF)
#' d
#'
#' mutate(d, newvar = cyl + hp)
#'
#' mutate_at(d, vars(starts_with("c")), ~.^2)
#'
#' group_by(d, cyl, am) %>%
#'    tally(gear)
#'
#' count(d, gear, am, cyl)
#'
#' select(d, am, cyl)
#'
#' select(d, am, cyl) %>%
#'    rename(foo = am)
#'
#' arrange(d, desc(hp))
#'
#' rbind(data.frame(mtcars[1, ], row.names = "MyCar"), DF) %>%
#'    add_dplyr_compat() %>%
#'    distinct()
#'
#' filter(d, am == 0)
#'
#' slice(d, 3:6)
add_dplyr_compat <- function(x) {
  if (inherits(x, "DataFrame"))
    dplyr::src(subclass = "DF", x = x)
  else
    x
}


#' Drop support as a `dplyr` backend
#'
#' In case something goes awry, this function can be used
#' to revert bak to just `S4Vectors::DataFrame` without `dplyr`
#' support.
#'
#' @param x A `S4Vectors::DataFrame` object
#'
#' @md
#' @export
#' @examples
#' library(S4Vectors)
#' library(dplyr)
#'
#' DF <- as(mtcars, "DataFrame")
#' d <- add_dplyr_compat(DF)
#' d
#'
#' drop_dplyr_compat(d)
drop_dplyr_compat <- function(x) {
  x$x
}

#' @export
format.src_DF <- function(x, ...) {
  x
}

#' @importFrom S4Vectors show
#' @export
print.src_DF <- function(x, ...) {
  S4Vectors::show(x)
}

setClass("src_DF", representation("VIRTUAL"))
setMethod("show", "src_DF", function(object) .show_src_DF(object))

#' @export
dplyr::filter
#' @importFrom digest digest
#' @export
filter.src_DF <- function(.data, ..., .preserve = FALSE) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  t$rowid <- seq_len(nrow(t))
  tf <- dplyr::filter(t, ..., .preserve = .preserve)
  tDF <- as(tf, "DataFrame")
  rownames(tDF) <- rn[tf$rowid]
  t$rowid <- NULL
  add_dplyr_compat(tDF)
}

#' @export
dplyr::mutate
#' @export
mutate.src_DF <- function(.data, ...) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  tm <- dplyr::mutate(t, ...)
  tDF <- as(tm, "DataFrame")
  rownames(tDF) <- rn
  add_dplyr_compat(tDF)
}

#' @export
dplyr::tbl_vars
#' @export
tbl_vars.src_DF <- function(x) {
  names(x$x)
}

#' @export
dplyr::select
#' @export
select.src_DF <- function(.data, ...) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  tm <- dplyr::select(t, ...)
  tDF <- as(tm, "DataFrame")
  rownames(tDF) <- rn
  add_dplyr_compat(tDF)
}

#' @export
dplyr::rename
#' @export
rename.src_DF <- function(.data, ...) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  tm <- dplyr::rename(t, ...)
  tDF <- as(tm, "DataFrame")
  rownames(tDF) <- rn
  add_dplyr_compat(tDF)
}

#' @export
dplyr::count
#' @export
count.src_DF <- function(.data, ...) {
  t <- convert_with_group(.data)
  dplyr::count(t, ...)
}

#' @export
dplyr::tally
#' @export
tally.src_DF <- function(x, wt = NULL, sort = FALSE, name = "n") {
  t <- convert_with_group(x)
  dplyr::tally(t, wt = wt, sort = sort, name = name)
}

#' @export
dplyr::summarise
dplyr::summarize
#' @export
summarise.src_DF <- function(.data, ...) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  dplyr::summarise(t, ...)
}

#' @export
dplyr::group_data
#' @export
group_data.src_DF <- function(.data) {
  attr(.data$x, "groups")
}

#' @export
dplyr::group_vars
#' @export
group_vars.src_DF <- function(x) {
  ## dplyr:::group_vars.grouped_df (not exported)
  groups <- group_data(x)[[1]]
  if (is.character(groups)) {
    groups
  }
  else if (is.data.frame(groups)) {
    head(names(groups), -1L)
  }
  else if (is.list(groups)) {
    purrr::map_chr(groups, rlang:::as_string)
  }
}

#' @export
dplyr::group_by
#' @export
group_by.src_DF <- function(.data, ..., add = FALSE, .drop = group_by_drop_default(.data)) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  tm <- dplyr::group_by(t, ..., add = add, .drop = .drop)
  groupdata <- group_data(tm)
  tDF <- as(tm, "DataFrame")
  rownames(tDF) <- rn
  attr(tDF, "groups") <- list(groupdata)
  add_dplyr_compat(tDF)
}

#' @export
dplyr::ungroup
#' @export
ungroup.src_DF <- function(x, ...) {
  attr(x$x, "groups") <- NULL
  attr(x$x@listData, "groups") <- NULL
  x
}

#' @export
dplyr::arrange
#' @export
arrange.src_DF <- function(.data, ...) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  ta <- dplyr::arrange(t, ...)
  tDF <- as(ta, "DataFrame")
  rownames(tDF) <- rn
  add_dplyr_compat(tDF)
}

#' @export
dplyr::distinct
#' @export
distinct.src_DF <- function(.data, ..., .keep_all = FALSE) {
  t <- convert_with_group(.data)
  td <- dplyr::distinct(t, ..., .keep_all = .keep_all)
  tDF <- as(td, "DataFrame")
  add_dplyr_compat(tDF)
}

#' @export
dplyr::pull
#' @export
pull.src_DF <- function(.data, var = -1) {
  var <- tidyselect::vars_pull(names(.data$x), !!enquo(var))
  .data$x[[var]]
}

#' @export
dplyr::slice
#' @export
slice.src_DF <- function(.data, ..., .preserve = FALSE) {
  rn <- rownames(.data$x)
  t <- convert_with_group(.data)
  ts <- dplyr::slice(t, ..., .preserve = .preserve)
  tDF <- as(ts, "DataFrame")
  rownames(tDF) <- rn[...]
  add_dplyr_compat(tDF)
}

convert_with_group <- function(.data) {
  t <- dplyr::tbl_df(.data$x)
  if (!is.null(group_data(.data))) {
    for (gvar in group_vars(.data)) {
      t <- dplyr::group_by(t, !!sym(gvar), add = TRUE)
    }
  }
  t
}
