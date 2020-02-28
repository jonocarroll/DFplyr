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

#' @importFrom dplyr filter
#' @export
dplyr::filter

#' @importFrom rlang quos eval_tidy quo_squash
#' @importFrom S4Vectors groupInfo subset
#' @export
filter.DataFrame <- function(.data, ..., .preserve = FALSE, quiet = FALSE, ungroup = FALSE) {
  FNS <- lapply(rlang::quos(...), rlang::quo_squash)
  if (length(FNS) > 1L && !quiet)
    message("Note that logial predicates will be applied in order provided.
Refer ?DFplyr")
  if (any(names(FNS) %in% names(.data))) {
    message("Arguments should not be named. Do you mean to use `==`?")
  }
  groupvars <- groupInfo(.data)
  if (length(groupvars) > 0L) {
    fact <- as.list(d[groupvars])
    split_data <- S4Vectors::split(.data, fact, drop = TRUE)
    for (f in seq_along(FNS)) {
      split_data <- lapply(split_data, function(xx) {
        with(xx, subset(xx, rlang::eval_tidy(FNS[[f]])))
      })
    }
    split_data <- S4Vectors::List(split_data)
    names(split_data) <- NULL
    .data <- DataFrame(split_data)
  } else {
    for (f in seq_along(FNS)) {
      .data <- with(.data, S4Vectors::subset(.data, rlang::eval_tidy(FNS[[f]])))
    }
    .data
  }
  if (ungroup) {
    .data <- ungroup(.data)
  }
  .data
}

top_n_rank <- function (n, wt) {
  if (n > 0) {
    min_rank(desc(wt)) <= n
  }
  else {
    min_rank(wt) <= abs(n)
  }
}

#' @importFrom dplyr mutate
#' @export
dplyr::mutate

#' @importFrom rlang quos quo_squash
#' @importFrom S4Vectors groupInfo subset
#' @export
mutate.DataFrame <- function(.data, ..., ungroup = FALSE) {

  FNS <- lapply(rlang::quos(...), rlang::quo_squash)
  groupvars <- S4Vectors::groupInfo(.data)
  if (length(groupvars) > 0L) {
    fact <- as.list(d[groupvars])
    split_data <- S4Vectors::split(.data, fact, drop = TRUE)
    split_data <- lapply(split_data, function(xx) {
      mutate_internal(xx, FNS, rlang::quos(...))
    })
    split_data <- S4Vectors::List(split_data)
    names(split_data) <- NULL
    .data <- BiocGenerics::unlist(split_data)
  } else {
    .data <-  mutate_internal(.data, FNS, rlang::quos(...))
  }
  if (ungroup) {
    .data <- ungroup(.data)
  }
  .data
}

#' @importFrom rlang eval_tidy quo_get_env
mutate_internal <- function(.data, FUNS, quos) {
  ## hack: inject the local env with the scoped data,
  ## excluding data that was already here.
  scope_env <- rlang::quo_get_env(quos[[1]])
  for (n in setdiff(ls(scope_env, all.names = TRUE), c("...", ls(all.names = TRUE)))) {
    assign(n, get(n, scope_env), pos = as.environment(-1L))
  }
  EXPRS <- lapply(names(FUNS), function(x) {
    FUNS_expl <- with(.data, rlang::eval_tidy(FUNS[[x]]))
    FUNS_obj <- with(.data, eval(rlang::eval_tidy(FUNS[[x]])))
    if (!inherits(FUNS_obj, "numeric")) {
      sprintf('%s <- %s', x, paste0(deparse(FUNS_expl), collapse = ""))
    } else {
      sprintf('%s <- c(%s)', x, paste0(FUNS_expl, collapse = ", "))
    }
  })
  S4Vectors::within(.data, eval(parse(text = paste0(unlist(EXPRS), collapse = '\n'))))
}


#' @importFrom dplyr tbl_vars
#' @export
dplyr::tbl_vars

#' @export
tbl_vars.DataFrame <- function(x) {
  names(x)
}

#' @importFrom dplyr select
#' @export
dplyr::select

#' @importFrom rlang quos quo_squash exprs
#' @export
select.DataFrame <- function(.data, ...) {
  dotnames <- names(rlang::exprs(...))
  .data <- S4Vectors::subset(.data,
                  select = unlist(lapply(
                    rlang::quos(...),
                    function(x){rlang::eval_tidy(rlang::quo_squash(x))})
                  ))
  if (any(dotnames != "")) {
    non_empty <- which(dotnames != "")
    for (ne in non_empty) {
      colnames(.data)[ne] <- dotnames[ne]
    }
  }
  .data
}

#' @importFrom dplyr rename
#' @export
dplyr::rename

#' @importFrom rlang quos quo_squash
#' @export
rename.DataFrame <- function(x, ...) {
  FNS <- lapply(rlang::quos(...), rlang::quo_squash)
  ## S4Vectors::rename not imported as it would mask the existing generic
  S4Vectors::rename(x, setNames(names(FNS), unlist(FNS)))
}

#' @importFrom dplyr count
#' @export
dplyr::count

#' @importFrom rlang quos quo_squash
#' @export
count <- function(x, ..., wt = NULL, sort = FALSE, name = "n", .drop = group_by_drop_default(x)) {

  ## count is not a generic, so preserve dplyr functionality
  if (!inherits(x, "DataFrame")) {
    return(dplyr::count(x, ..., wt = !!enquo(wt), sort = sort, name = name, .drop = .drop))
  }

  EXPRS <- lapply(rlang::quos(...), function(x) rlang::quo_squash(x))
  groupvars <- S4Vectors::groupInfo(x)
  groups <- group_combos(x)
  if (length(groupvars) > 0L) {
    split_data <- S4Vectors::split(x, x[groupvars])
    split_data <- lapply(split_data, function(xx) {
      if (nrow(xx) == 0L)
        return(NULL)
      tbl_grp <- S4Vectors::with(xx, do.call(table, EXPRS))
      unique_grp <- unique(xx[, groupvars])
      rownames(unique_grp) <- NULL
      cbind(unique_grp, as(tbl_grp, "data.frame"))
    })
    .data <- do.call(rbind, split_data)
  } else {
    tbl_grp <- with(x, do.call(table, EXPRS))
    .data <- as(tbl_grp, "DataFrame")
    .data <- as.data.frame(.data)
  }
  .data <- ungroup(.data)
  names(.data)[ncol(.data)] <- name
  .data <- .data[.data$n != 0, ]
  .data <- .data[with(.data, do.call(order, rlang::syms(EXPRS))), ]
  rownames(.data) <- seq_len(nrow(.data))
  .data
}

#' @importFrom dplyr group_by_drop_default
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

#' @importFrom dplyr summarise summarize
#' @export
dplyr::summarise
dplyr::summarize

#' @importFrom rlang quos quo_squash eval_tidy
#' @export
summarise.DataFrame <- function(.data, ...) {
  FNS <- lapply(rlang::quos(...), rlang::quo_squash)
  groupvars <- groupInfo(.data)

  if (length(groupvars) > 0L) {
    fact <- as.list(d[groupvars])
    groups <- group_combos(.data)
    split_data <- S4Vectors::split(.data, fact, drop = TRUE)
      split_data <- lapply(seq_len(length(split_data)), function(xx) {
        tbl_grp <- lapply(FNS, function(f) with(split_data[[xx]], rlang::eval_tidy(f)))
        cbind(groups[xx,], as(tbl_grp, "DataFrame"))
      })
    split_data <- S4Vectors::List(split_data)
    names(split_data) <- NULL
    RET <- BiocGenerics::unlist(split_data)



    # split_data <- lapply(seq_len(nrow(groups)), function(xx) {
    #   .data_grp <- .data[groups$.rows[xx][[1]], ,drop = FALSE]
    #   tbl_grp <- lapply(FNS, function(xx) {with(.data_grp, rlang::eval_tidy(xx))})
    #   cbind(groups, as(tbl_grp, "DataFrame"))
    # })
    # RET <- do.call(rbind, split_data)
  } else {
    RET <- as(lapply(FNS, function(xx) {with(.data, eval(xx))}), "DataFrame")
  }
  RET
}

#' @export
summarize.DataFrame <- summarise.DataFrame

#' @importFrom dplyr group_data
#' @export
dplyr::group_data

#' @importFrom tibble tibble
#' @export
group_data.DataFrame <- function(.data) {
  group_attr <- attr(.data@listData, "groups")
  if (!is.null(group_attr) && nrow(group_attr) > 1L) {
    group_attr
  } else {
    rows <- list(seq_len(nrow(.data)))
    tibble::tibble(`:=`(".rows", rows))
  }
}

#' @importFrom dplyr group_vars
#' @export
dplyr::group_vars

#' @importFrom rlang as_string
#' @export
group_vars.DataFrame <- function(x) {
  groups <- group_data(x)
  if (is.character(groups)) {
    groups
  }
  else if (is.data.frame(groups)) {
    head(names(groups), -1L)
  }
}

group_combos <- function(.data, groupvars = NULL) {
  if (is.null(groupvars))
    groupvars <- S4Vectors::groupInfo(.data)
  if (is.null(groupvars))
    return(NULL)
  uniques <- unique(select(.data, !!!syms(unlist(groupvars))))
  uniques <- as.data.frame(uniques)
  rownames(uniques) <- NULL
  uniques[with(uniques, do.call(order, rlang::syms(groupvars))), ]
}

#' @importFrom dplyr group_by
#' @export
dplyr::group_by

#' @importFrom rlang quos quo_squash as_string syms
#' @importFrom tibble as.tibble
#' @export
group_by.DataFrame <- function(.data, ..., add = FALSE, .drop = group_by_drop_default(.data)) {

  groupvars <- vapply(rlang::quos(...),
                      function(x) rlang::as_string(rlang::quo_squash(x)),
                      character(1))
  groupvars <- intersect(groupvars, names(.data))
  if (is.null(groupInfo(.data)) && !add) {
    if (!length(groupvars)) {
      return(.data)
    } else {
      groupInfo(.data) <- unname(groupvars)
    }
  } else if (!is.null(groupInfo) && length(groupvars) && add) {
    groupInfo(.data) <- c(groupInfo(.data), groupvars)
  }
  .data
}

#' @importFrom dplyr ungroup
#' @export
dplyr::ungroup

#' @export
ungroup.DataFrame <- function(x, ...) {
  S4Vectors::groupInfo(x) <- NULL
  x
}

#' @importFrom dplyr arrange
#' @export
dplyr::arrange

#' @export
dplyr::arrange
#' @export
arrange.DataFrame <- function(.data, ...) {

  EXPRS <- lapply(rlang::quos(...), function(x) rlang::quo_squash(x))

  groupvars <- group_vars(.data)
  if (length(groupvars) > 0L) {
    groups <- group_data(.data)
    split_data <- lapply(seq_len(nrow(groups)), function(x) {
      .data_grp <- .data[groups$.rows[x][[1]], ,drop = FALSE]
      .data_grp[with(.data_grp, do.call(order, EXPRS)), ]
    })
    do.call(rbind, split_data)
  } else {
    .data[with(.data, do.call(order, EXPRS)), ]
  }
}

#' @export
desc <- function(x) {
  -xtfrm(x)
}

#' @importFrom dplyr distinct
#' @export
dplyr::distinct

#' @export
dplyr::distinct
#' @export
distinct.DataFrame <- function(.data, ..., .keep_all = FALSE) {
  .data[!BiocGenerics::duplicated(.data), ]
}

#' @importFrom dplyr pull
#' @export
dplyr::pull

#' @export
dplyr::pull
#' @export
pull.DataFrame <- function(.data, var = -1, name = NULL) {
  var <- tidyselect::vars_pull(names(.data), !!enquo(var))
  name <- rlang::enquo(name)
  if (rlang::quo_is_null(name)) {
    return(.data[[var]])
  }
  name <- tidyselect::vars_pull(names(.data), !!name)
  rlang::set_names(.data[[var]], nm = .data[[name]])
}

#' @importFrom dplyr slice
#' @export
dplyr::slice

#' @export
dplyr::slice
#' @export
slice.DataFrame <- function(.data, ..., .preserve = FALSE) {
  .data[..., ]
}

#' @keywords internal
convert_with_group <- function(.data) {
  t <- tibble::as_tibble(.data)
  if (!is.null(group_data(.data)) && nrow(group_data(.data)) > 1L) {
    for (gvar in group_vars(.data)) {
      t <- dplyr::group_by(t, !!sym(gvar), add = TRUE)
    }
  }
  t
}

#' @keywords internal
restore_DF <- function(.data, rn) {
  DF <- as(.data, "DataFrame")
  rownames(DF) <- rn
  DF
}
