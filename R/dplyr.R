utils::globalVariables(c(".data", ":="))

#' @inherit base::format
#' @export
format.DataFrame <- function(x, ...) {
    x
}

#' @inherit dplyr::filter
#' @importFrom rlang quos eval_tidy
#' @export
filter.DataFrame <- function(.data, ..., .preserve = FALSE) {
    FNS <- lapply(rlang::quos(...), rlang::quo_squash)
    groupvars <- group_vars(.data)
    if (length(groupvars) > 0L) {
        groups <- group_data(.data)
        split_data <- lapply(seq_len(nrow(groups)), function(x) {
            .data_grp <- .data[groups$.rows[x][[1]], , drop = FALSE]
            for (f in seq_along(FNS)) {
                .data_grp <- with(
                    .data_grp,
                    base::subset(
                        .data_grp,
                        rlang::eval_tidy(FNS[[f]])
                    )
                )
            }
            .data_grp
        })
        do.call(rbind, split_data)
    } else {
        for (f in seq_along(FNS)) {
            .data <- with(
                .data,
                base::subset(
                    .data,
                    rlang::eval_tidy(FNS[[f]])
                )
            )
        }
        .data
    }
}

#' @inherit dplyr::mutate
#' @importFrom rlang quos quo_squash
#' @export
mutate.DataFrame <- function(.data, ...) {
    FNS <- lapply(rlang::quos(...), rlang::quo_squash)

    groupvars <- group_vars(.data)
    if (length(groupvars) > 0L) {
        groups <- group_data(.data)
        split_data <- lapply(seq_len(nrow(groups)), function(x) {
            .data_grp <- .data[groups$.rows[x][[1]], , drop = FALSE]
            mutate_internal(.data_grp, FNS, rlang::quos(...))
        })
        do.call(rbind, split_data)
    } else {
        mutate_internal(.data, FNS, rlang::quos(...))
    }
}

#' @keywords internal
#' @importFrom rlang eval_tidy quo_get_env
mutate_internal <- function(.data, FUNS, quos) {
    op <- options("useFancyQuotes")
    on.exit(options(op))
    options(useFancyQuotes = FALSE)

    ## hack: inject the local env with the scoped data,
    ## excluding data that was already here.
    ## this is important for capturing RHS variables
    scope_env <- rlang::quo_get_env(quos[[1]])
    for (n in setdiff(
        ls(scope_env, all.names = TRUE),
        c("...", ls(all.names = TRUE))
    )) {
        assign(n, get(n, scope_env), pos = as.environment(-1L))
    }
    EXPRS <- lapply(names(FUNS), function(x) {
        FUNS_expl <- with(.data, rlang::eval_tidy(FUNS[[x]]))
        FUNS_obj <- with(.data, eval(rlang::eval_tidy(FUNS[[x]])))
        if (!inherits(FUNS_obj, "numeric")) {
            sprintf("%s <- %s", x, paste0(deparse(FUNS_expl), collapse = ""))
        } else {
            sprintf("%s <- c(%s)", x, paste0(FUNS_expl, collapse = ", "))
        }
    })
    S4Vectors::within(
        .data,
        eval(
            parse(
                text = paste0(unlist(EXPRS), collapse = "\n")
            )
        )
    )
}

#' @inherit dplyr::tbl_vars
#' @return all variables, with a `groups` attribute when grouped.
#' @export
tbl_vars.DataFrame <- function(x) {
    names(x)
}

#' @inherit dplyr::select
#' @importFrom rlang quos quo_squash exprs eval_tidy
#' @export
select.DataFrame <- function(.data, ...) {
    dotnames <- names(rlang::exprs(...))
    .data <- base::subset(
        .data,
        select = unlist(lapply(
            rlang::quos(...),
            function(x) {
                rlang::eval_tidy(rlang::quo_squash(x))
            }
        ))
    )
    if (any(dotnames != "")) {
        non_empty <- which(dotnames != "")
        for (ne in non_empty) {
            colnames(.data)[ne] <- dotnames[ne]
        }
    }
    .data
}

#' @importFrom rlang quos quo_squash set_names
#' @importFrom S4Vectors rename
#'
.rename <- function(x, ...) {
    FNS <- lapply(rlang::quos(...), rlang::quo_squash)
    methods::getMethod(
        "rename",
        "Vector",
        where = asNamespace("S4Vectors")
    )(x, rlang::set_names(names(FNS), unlist(FNS)))
}

#' Rename Columns of a DataFrame
#' @param x a DataFrame
#' @param ... NSE syntax; `new_name = old_name` to rename selected variables
#' @export
setMethod("rename", "DataFrame", .rename)

#' Rename Columns of a DataFrame (deprecated)
#' @param .data a DataFrame
#' @param ... columns to be renamed with syntax `new = old`
#' @export
rename2 <- function(.data, ...) {
    .Deprecated(
        "rename",
        msg = "DFplyr now properly supports rename with NSE syntax"
    )
}

#' @inherit dplyr::count
#' @importFrom rlang quos quo_squash enquo
#' @export
count.DataFrame <- function(
        x,
        ...,
        wt = NULL,
        sort = FALSE,
        name = "n",
        .drop = group_by_drop_default(x)) {
    if (!inherits(x, "DataFrame")) {
        return(
            dplyr::count(
                x,
                ...,
                wt = !!rlang::enquo(wt),
                sort = sort,
                name = name,
                .drop = .drop
            )
        )
    }

    groupvars <- group_vars(x)
    EXPRS <- lapply(rlang::quos(...), function(x) rlang::quo_squash(x))
    if (length(groupvars) > 0L) {
        groups <- group_data(x)
        if (!length(EXPRS)) {
            RET <- select(
                mutate(groups, n = lengths(.data[[".rows"]])),
                -".rows"
            )
            names(RET)[ncol(RET)] <- name
            RET <- RET[RET[[name]] != 0, ]
            return(methods::as(RET, "DataFrame"))
        }
        split_data <- lapply(seq_len(nrow(groups)), function(xx) {
            .data_grp <- x[groups$.rows[xx][[1]], , drop = FALSE]
            tbl_grp <- with(.data_grp, do.call(table, EXPRS))
            cbind(groups[xx, -ncol(groups)], methods::as(tbl_grp, "DataFrame"))
        })
        RET <- methods::as(do.call(rbind, split_data), "DataFrame")
    } else {
        RET <- methods::as(with(x, do.call(table, EXPRS)), "DataFrame")
    }

    names(RET)[ncol(RET)] <- name
    RET <- RET[RET[[name]] != 0, ]
    RET <- RET[with(RET, do.call(order, EXPRS)), ]

    RET
}


#' @inherit dplyr::group_by_drop_default
#' @export
group_by_drop_default.DataFrame <- function(.tbl) {
    if (!is.null(group_data(.tbl)) && nrow(group_data(.tbl)) > 1L) {
        tryCatch(
            {
                !identical(attr(group_data(.tbl), ".drop"), FALSE)
            },
            error = function(e) {
                TRUE
            }
        )
    } else {
        TRUE
    }
}

#' @inherit dplyr::summarise
#' @importFrom rlang quos quo_squash eval_tidy
#' @export
summarise.DataFrame <- function(.data, ...) {
    FNS <- lapply(rlang::quos(...), rlang::quo_squash)
    groupvars <- group_vars(.data)

    if (length(groupvars) > 0L) {
        groups <- group_data(.data)
        split_data <- lapply(seq_len(nrow(groups)), function(xx) {
            .data_grp <- .data[groups$.rows[xx][[1]], , drop = FALSE]
            tbl_grp <- lapply(FNS, function(xx) {
                with(.data_grp, rlang::eval_tidy(xx))
            })
            cbind(groups[xx, -ncol(groups)], methods::as(tbl_grp, "DataFrame"))
        })
        RET <- do.call(rbind, split_data)
    } else {
        RET <- lapply(FNS, function(xx) {
            with(.data, eval(xx))
        })
    }
    RET <- methods::as(RET, "DataFrame")
    RET
}

#' @inherit dplyr::summarize
#' @export
summarize.DataFrame <- summarise.DataFrame

#' @inherit dplyr::group_data title description
#' @param .data a [S4Vectors::DataFrame()]
#' @return a `data.frame` of group data
#' @export
group_data.DataFrame <- function(.data) {
    # group_attr <- attr(.data@listData, "groups")
    group_attr <- get_group_data(.data)
    if (!is.null(group_attr) && nrow(group_attr) > 1L) {
        group_attr
    } else {
        rows <- list(seq_len(nrow(.data)))
        data.frame(.rows = I(rows))
    }
}

#' @inherit dplyr::group_vars title description
#' @importFrom rlang as_string
#' @param x a [S4Vectors::DataFrame()], likely grouped
#' @return the grouping variables as a character vector
#' @export
group_vars.DataFrame <- function(x) {
    groups <- group_data(x)
    if (is.character(groups)) {
        groups
    } else if (is.data.frame(groups)) {
        head(names(groups), -1L)
    }
}

#' @inherit dplyr::group_by
#' @importFrom rlang quos quo_squash as_string syms
#' @export
group_by.DataFrame <- function(
        .data,
        ...,
        add = FALSE,
        .drop = group_by_drop_default(.data)) {
    if (is.null(group_data(.data)) || nrow(group_data(.data)) == 1L) {
        groupvars <- lapply(
            rlang::quos(...),
            function(x) rlang::as_string(rlang::quo_squash(x))
        )
        for (v in groupvars) {
            if (!utils::hasName(.data, v)) stop("Column '", v, "' not found in data")
        }
        uniques <- unique(select(.data, !!!rlang::syms(unlist(groupvars))))
        flagged <- S4Vectors::merge(
            mutate(.data, rowid = seq_len(nrow(.data))),
            mutate(uniques, flag = seq_len(nrow(uniques))),
            by = unlist(groupvars),
            sort = FALSE
        )
        groups <- split(as.integer(flagged$rowid), flagged$flag)
        uniques <- as.data.frame(uniques)
        uniques$.rows <- unname(groups)
        groupdata <- uniques[
            with(
                uniques,
                do.call(order, rlang::syms(groupvars))
            ),
        ]
        if (!inherits(groupdata, "data.frame")) {
            return(.data)
        }
        rownames(groupdata) <- seq_len(nrow(groupdata))
        .data <- set_group_data(.data, groupdata, .drop)
    }
    .data
}

#' Set and Get Group Data on a DataFrame
#'
#' The location of group data is an internal implemnetation
#' detail, so these get and set methods enable interfacing
#' with that data.
#'
#' @rdname groupdata
#' @name group_data
#' @param x A [S4Vectors::DataFrame] on which to set group data.
#' @param g Group data (a `data.frame`).
#' @param .drop Drop groups formed by factor levels that don't appear in
#'     the data?
#' @return For `set_group_data`, the input `x` with group data set as metadata.
#'     For `get_group_data`, the group data that is set on `x`.
#' @keywords internal
set_group_data <- function(x, g, .drop = group_by_drop_default(x)) {
    if (!is.null(g)) {
        attr(g, ".drop") <- .drop
    }
    metadata(x)$groups <- g
    x
}

#' @rdname groupdata
#' @keywords internal
get_group_data <- function(x) {
    metadata(x)$groups
}

#' @inherit dplyr::ungroup
#' @export
ungroup.DataFrame <- function(x, ...) {
    set_group_data(x, NULL)
}

#' @inherit dplyr::arrange
#' @export
arrange.DataFrame <- function(.data, ...) {
    EXPRS <- lapply(rlang::quos(...), function(x) rlang::quo_squash(x))

    groupvars <- group_vars(.data)
    if (length(groupvars) > 0L) {
        groups <- group_data(.data)
        split_data <- lapply(seq_len(nrow(groups)), function(x) {
            .data_grp <- ungroup(.data)[groups$.rows[x][[1]], , drop = FALSE]
            .data_grp[with(.data_grp, do.call(order, EXPRS)), ]
        })
        group_by(do.call(rbind, split_data), !!!groupvars)
    } else {
        .data[with(.data, do.call(order, EXPRS)), ]
    }
}


#' @inherit dplyr::desc
#' @return the input vector in a format that will be sorted in descending order.
#' @export
desc <- function(x) {
    -xtfrm(x)
}

#' @inherit dplyr::distinct
#' @importFrom BiocGenerics duplicated
#' @export
distinct.DataFrame <- function(.data, ..., .keep_all = FALSE) {
    .data[!BiocGenerics::duplicated(.data), ]
}

#' @inherit dplyr::pull
#' @export
#' @examples
#' d <- S4Vectors::DataFrame(mtcars)
#' pull(d, cyl)
#'
pull.DataFrame <- function(.data, var = -1, name = NULL, ...) {
    var <- tidyselect::vars_pull(names(.data), !!rlang::enquo(var))
    name <- rlang::enquo(name)
    if (rlang::quo_is_null(name)) {
        return(.data[[var]])
    }
    name <- tidyselect::vars_pull(names(.data), !!name)
    rlang::set_names(.data[[var]], nm = .data[[name]])
}

#' @inherit dplyr::slice
#' @export
slice.DataFrame <- function(.data, ..., .preserve = FALSE) {
    groupvars <- group_vars(.data)
    if (length(groupvars) > 0L) {
        groups <- group_data(.data)
        split_data <- lapply(seq_len(nrow(groups)), function(x) {
            .data_grp <- .data[groups$.rows[x][[1]], , drop = FALSE]
            .data_grp[..., ]
        })
        do.call(rbind, split_data)
    } else {
        .data[..., ]
    }
}

# #' @keywords internal
# convert_with_group <- function(.data) {
#     t <- tibble::as_tibble(.data)
#     if (!is.null(group_data(.data)) && nrow(group_data(.data)) > 1L) {
#         for (gvar in group_vars(.data)) {
#             t <- dplyr::group_by(t, !!rlang::sym(gvar), add = TRUE)
#         }
#     }
#     t
# }

# #' @keywords internal
# restore_DF <- function(.data, rn) {
#     DF <- methods::as(.data, "DataFrame")
#     rownames(DF) <- rn
#     DF
# }

#' @inherit dplyr::tally
#' @export
tally.DataFrame <- function(x, wt = NULL, sort = FALSE, name = NULL) {
    name <- .check_n_name(name, group_vars(x))
    out <- .tally_n(x, {{ wt }}, name)
    if (sort) {
        arrange(out, desc(!!rlang::sym(name)))
    } else {
        out
    }
}

#' @importFrom rlang caller_arg caller_env inform
#' @keywords internal
.check_n_name <- function(
        name,
        vars,
        arg = rlang::caller_arg(name),
        call = rlang::caller_env()) {
    if (is.null(name)) {
        name <- .n_name(vars)
        if (name != "n") {
            rlang::inform(c(
                paste0(
                    "Storing counts in `",
                    name,
                    "`,
                    as `n` already present in input"
                ),
                i = "Use `name = \"new_name\"` to pick a new name."
            ))
        }
    }
    name
}

#' @keywords internal
.n_name <- function(x) {
    name <- "n"
    while (name %in% x) {
        name <- paste0("n", name)
    }
    name
}

#' @importFrom rlang enquo quo_get_expr warn quo_is_null expr
#' @keywords internal
.tally_n <- function(x, wt, name) {
    wt <- rlang::enquo(wt)
    if (rlang::is_call(rlang::quo_get_expr(wt), "n", n = 0)) {
        rlang::warn(c(
            "`wt = n()` is deprecated",
            i = "You can now omit the `wt` argument"
        ))
        wt <- rlang::quo(NULL)
    }
    if (rlang::quo_is_null(wt)) {
        count(x, sort = sort, name = name)
    } else {
        n <- rlang::expr(sum(!!wt, na.rm = TRUE))
        summarise(x, `:=`(!!name, !!n))
    }
}

group_intersect <- function(x, new) {
    intersect(group_vars(x), names(new))
}


.s4_subset <- function(x, i, j, ..., drop = TRUE) {
    if (!isTRUEorFALSE(drop)) stop("'drop' must be TRUE or FALSE")
    if (length(list(...)) > 0L) warning("parameters in '...' not supported")

    ## NOTE: matrix-style subsetting by logical matrix not supported.
    list_style_subsetting <- (nargs() - !missing(drop)) < 3L
    if (list_style_subsetting || !missing(j)) {
        if (list_style_subsetting) {
            if (!missing(drop)) {
                warning("'drop' argument ignored by list-style subsetting")
            }
            if (missing(i)) {
                return(x)
            }
            j <- i
        }
        x <- S4Vectors::extractCOLS(x, j)
        if (list_style_subsetting) {
            return(x)
        }
    }
    if (!missing(i)) x <- S4Vectors::extractROWS(x, i)
    if (missing(drop)) {
        # drop by default if only one column left
        drop <- ncol(x) == 1L
    }
    if (drop) {
        ## one column left
        if (ncol(x) == 1L) {
            return(x[[1L]])
        }
        ## one row left
        if (nrow(x) == 1L) {
            return(methods::as(x, "list"))
        }
    }
    x
}


.grp_subset <- function(x, i, j, ..., drop = FALSE) {
    if (is.null(get_group_data(x))) {
        return(.s4_subset(x, i, j, ..., drop = drop))
    }
    out <- ungroup(x)[i, j, drop = drop]
    if (drop) {
        out
    } else {
        groups <- group_intersect(x, out)
        if ((missing(i) || nargs() == 2) && identical(groups, group_vars(x))) {
            out
        } else {
            group_by(out, !!!groups)
        }
    }
}

#' Subset a DataFrame
#'
#' @param x a `DataFrame` to be subset
#' @param i rows to subset
#' @param j columns to subset
#' @param ... other params, passed to regular `S4Vectors` subsetting
#' @param drop drop dimensions?
#' @export
setMethod("[", "DataFrame", .grp_subset)

.grp_bindROWS <- function(x, objects = list()) {
    combined <- methods::getMethod(
        "bindROWS",
        "DataFrame",
        where = asNamespace("S4Vectors")
    )(
        x,
        objects
    )
    if (is.null(group_data(x))) {
        return(combined)
    }
    group_by(ungroup(combined), !!!group_vars(x))
}

#' rbind DataFrames
#' @param x a DataFrame
#' @param objects a list of DataFrames
#' @export
setMethod("bindROWS", "DataFrame", .grp_bindROWS)
