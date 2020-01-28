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

#' @export
format.DataFrame <- function(x, ...) {
  x
}

S4_cols <- function(x) {
  vapply(x, isS4, logical(1))
}

has_S4 <- function(x) {
  any(S4_cols(x))
}

restore_S4 <- function(old, new, id) {
  ## restore S4 columns
  s4 <- names(which(S4_cols(old)))

  for (s4i in s4) {
    broken_cols <- if (rlang::has_name(new, s4i)) {
      s4i
    } else {
      names(new)[startsWith(names(new), paste0(s4i, "."))]
    }
    new[[broken_cols[1]]] <- old[[s4i]][id]
    names(new) <- replace(names(new), which(names(new) == broken_cols[[1]]), s4i)
    if (length(broken_cols) > 1L) {
      for (bc in broken_cols[2:length(broken_cols)]) {
        new[[bc]] <- NULL
      }
    }
  }

  return(new)
}

#' @export
dplyr::filter
#' @importFrom digest digest
#' @export
filter.DataFrame <- function(.data, ..., .preserve = FALSE) {

  t <- convert_with_group(.data)
  t$rowid <- seq_len(nrow(t))
  tf <- dplyr::filter(t, ..., .preserve = .preserve)
  tDF <- restore_DF(tf, rownames(.data)[tf$rowid])

  if (has_S4(.data)) {
    tDF <- restore_S4(.data, tDF, tf$rowid)
  }
  tDF$rowid <- NULL

  tDF
}


#' @export
dplyr::mutate
#' @export
mutate.DataFrame <- function(.data, ...) {
  t <- convert_with_group(.data)
  tm <- dplyr::mutate(t, ...)
  tDF <- restore_DF(tm, rownames(.data))

  if (has_S4(.data)) {
    tDF <- restore_S4(.data, tDF, seq_len(nrow(tm)))
  }

  tDF
}


#' @export
dplyr::tbl_vars
#' @export
tbl_vars.DataFrame <- function(x) {
  names(x)
}


#' @export
dplyr::select
#' @export
select.DataFrame <- function(.data, ...) {
  t <- convert_with_group(.data)
  ts <- dplyr::select(t, ...)
  restore_DF(ts, rownames(.data))
}


#' @export
dplyr::rename
#' @export
rename.DataFrame <- function(.data, ...) {
  t <- convert_with_group(.data)
  tr <- dplyr::rename(t, ...)
  restore_DF(tr, rownames(.data))
}


#' @export
dplyr::group_by_drop_default
#' @export
group_by_drop_default.DataFrame <- function(.tbl) {
  if (!is.null(group_data(.tbl)) && nrow(group_data(.tbl)) > 1L) {
    tryCatch({
      !identical(attr(group_data(.tbl), ".drop"), FALSE)
    }, error = function(e) {
      TRUE
    })
  } else {
    TRUE
  }
}


#' @export
dplyr::tally
#' @export
tally.DataFrame <- function(x, wt = NULL, sort = FALSE, name = "n") {
  t <- convert_with_group(x)
  wt <- enquo(wt)
  dplyr::tally(t, wt = !!(wt), sort = sort, name = name)
}


#' @export
dplyr::summarise
dplyr::summarize
#' @export
summarise.DataFrame <- function(.data, ...) {
  t <- convert_with_group(.data)
  dplyr::summarise(t, ...)
}



#' @export
dplyr::group_data
#' @export
group_data.DataFrame <- function(.data) {
  group_attr <- attr(.data@listData, "groups")
  if (!is.null(group_attr) && nrow(group_attr) > 1L) {
    group_attr
  } else {
    rows <- list(seq_len(nrow(.data)))
    tibble(`:=`(".rows", rows))
  }
}


#' @export
dplyr::group_vars
#' @export
group_vars.DataFrame <- function(x) {
  groups <- group_data(x)
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
group_by.DataFrame <- function(.data, ..., add = FALSE, .drop = group_by_drop_default(.data)) {
  t <- convert_with_group(.data)
  tm <- dplyr::group_by(t, ..., add = add, .drop = .drop)
  restore_DF(tm, rownames(.data))
}


#' @export
dplyr::ungroup
#' @export
ungroup.DataFrame <- function(x, ...) {
  attr(x@listData, "groups") <- NULL
  x
}


#' @export
dplyr::arrange
#' @export
arrange.DataFrame <- function(.data, ...) {
  t <- convert_with_group(.data)
  ta <- dplyr::arrange(t, ...)
  restore_DF(ta, rownames(.data))
}


#' @export
dplyr::distinct
#' @export
distinct.DataFrame <- function(.data, ..., .keep_all = FALSE) {
  t <- convert_with_group(.data)
  td <- dplyr::distinct(t, ..., .keep_all = .keep_all)
  restore_DF(td, NULL) # no rownames since they can't be determined
}


#' @export
dplyr::pull
#' @export
pull.DataFrame <- function(.data, var = -1) {
  var <- tidyselect::vars_pull(names(.data), !!enquo(var))
  .data[[var]]
}


#' @export
dplyr::slice
#' @export
slice.DataFrame <- function(.data, ..., .preserve = FALSE) {
  t <- convert_with_group(.data)
  ts <- dplyr::slice(t, ..., .preserve = .preserve)
  restore_DF(ts, rownames(.data)[...])
}


convert_with_group <- function(.data) {
  t <- dplyr::tbl_df(.data)
  if (!is.null(group_data(.data)) && nrow(group_data(.data)) > 1L) {
    for (gvar in group_vars(.data)) {
      t <- dplyr::group_by(t, !!sym(gvar), add = TRUE)
    }
  }
  t
}

restore_DF <- function(.data, rn) {
  DF <- as(.data, "DataFrame")
  rownames(DF) <- rn
  DF
}
