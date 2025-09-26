.show_DF <- function(object) {
    nhead <- S4Vectors::get_showHeadLines()
    ntail <- S4Vectors::get_showTailLines()
    x_nrow <- nrow(object)
    x_ncol <- ncol(object)
    cat(
        S4Vectors::classNameForDisplay(object),
        " with ", x_nrow, " row",
        ifelse(x_nrow == 1L, "", "s"),
        " and ", x_ncol, " column",
        ifelse(x_ncol == 1L, "", "s"),
        "\n",
        sep = ""
    )
    ## add group information
    if (!is.null(group_data(object)) & nrow(group_data(object)) > 1L) {
        cat("Groups: ", toString(group_vars(object)), "\n")
    }
    if (x_nrow != 0L && x_ncol != 0L) {
        x_rownames <- rownames(object)
        if (x_nrow <= nhead + ntail + 1L) {
            m <- S4Vectors::makeNakedCharacterMatrixForDisplay(object)
            if (!is.null(x_rownames)) {
                rownames(m) <- x_rownames
            }
        } else {
            m <- rbind(
                S4Vectors::makeNakedCharacterMatrixForDisplay(
                    head(ungroup(object), nhead)
                ),
                rbind(rep.int("...", x_ncol)),
                S4Vectors::makeNakedCharacterMatrixForDisplay(
                    tail(ungroup(object), ntail)
                )
            )
            rownames(m) <- .make_rownames_for_RectangularData_display(
                x_rownames,
                x_nrow,
                nhead,
                ntail
            )
        }
        m <- rbind(.make_class_info_for_DataFrame_display(object), m)
        print(m, quote = FALSE, right = TRUE)
    }
    invisible(NULL)
}

setMethod("show", "DataFrame", .show_DF)

.make_rownames_for_RectangularData_display <- function(
        x_rownames,
        nrow,
        nhead,
        ntail) {
    p1 <- ifelse(nhead == 0L, 0L, 1L)
    p2 <- ifelse(ntail == 0L, 0L, ntail - 1L)
    s1 <- s2 <- character(0)
    if (is.null(x_rownames)) {
        if (nhead > 0) {
            s1 <- paste0(as.character(p1:nhead))
        }
        if (ntail > 0) {
            s2 <- paste0(as.character((nrow - p2):nrow))
        }
    } else {
        if (nhead > 0) {
            s1 <- paste0(head(x_rownames, nhead))
        }
        if (ntail > 0) {
            s2 <- paste0(tail(x_rownames, ntail))
        }
    }
    c(s1, "...", s2)
}

.make_class_info_for_DataFrame_display <- function(x) {
    matrix(
        unlist(
            lapply(x, function(col) {
                paste0("<", .classNameForDisplay(col), ">")
            }),
            use.names = FALSE
        ),
        nrow = 1L,
        dimnames = list("", colnames(x))
    )
}

.classNameForDisplay <- function(x) class(x)[1L]

# borrowed from S4Vectors (not exported)
.new_DataFrame <- function(listData = list(), nrows = NA, what = "arguments") {
    stopifnot(is.list(listData))
    stopifnot(S4Vectors::isSingleNumberOrNA(nrows))
    if (!is.integer(nrows)) nrows <- as.integer(nrows)
    listData_nrow <- nrow(listData)
    if (is.null(listData_nrow)) {
        if (length(listData) == 0L) {
            if (is.na(nrows)) nrows <- 0L
            names(listData) <- character(0)
        } else {
            if (is.na(nrows)) {
                elt_nrows <- S4Vectors::elementNROWS(listData)
                nrows <- elt_nrows[[1L]]
                if (!all(elt_nrows == nrows)) {
                    stop(
                        S4Vectors::wmsg(
                            what,
                            " imply differing number of rows"
                        )
                    )
                }
            }
            if (is.null(names(listData))) {
                names(listData) <- paste0("V", seq_along(listData))
            }
        }
    } else {
        if (is.na(nrows)) {
            nrows <- listData_nrow
        } else if (nrows != listData_nrow) {
            stop(wmsg(
                "the supplied 'nrows' does not match ",
                "the nb of rows in 'listData'"
            ))
        }
        listData <- as.list(listData)
    }
    new2("DFrame", nrows = nrows, listData = listData, check = FALSE)
}


setAs("grouped_df", "DFrame", function(from) {
    grps <- dplyr::group_vars(from)
    from <- dplyr::ungroup(from)
    rn <- attributes(from)[["row.names"]]
    if (is.integer(rn)) rn <- NULL
    nr <- nrow(from)
    attr(from, "row.names") <- NULL
    class(from) <- NULL
    ans <- .new_DataFrame(from, nrows = nr)
    ans@rownames <- rn
    group_by(ans, !!!grps)
})
