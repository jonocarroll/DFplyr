
#' @import S4Vectors
.show_DF <- function(object) {
  x <- x_src_DF <- object
  nhead <- get_showHeadLines()
  ntail <- get_showTailLines()
  x_nrow <- nrow(x)
  x_ncol <- ncol(x)
  cat("dplyr-compatible DataFrame", " with ",
      x_nrow, " row", ifelse(x_nrow == 1L, "", "s"),
      " and ",
      x_ncol, " column", ifelse(x_ncol == 1L, "", "s"),
      "\n", sep="")
  if (!is.null(group_data(x_src_DF)) & nrow(group_data(x_src_DF)) > 1L) {
    cat("Groups: ", toString(group_vars(x_src_DF)), "\n")
  }
  if (x_nrow != 0L && x_ncol != 0L) {
    x_rownames <- rownames(x)
    if (x_nrow <= nhead + ntail + 1L) {
      m <- as.matrix(x)
      if (!is.null(x_rownames))
        rownames(m) <- x_rownames
    } else {
      m <- rbind(as.matrix(head(x, nhead)),
                 rbind(rep.int("...", x_ncol)),
                 as.matrix(tail(x, ntail)))
      rownames(m) <- .make_rownames_for_display(
        x_rownames, x_nrow,
        nhead, ntail)
    }
    m <- rbind(.make_class_info_for_display(x), m)
    print(m, quote=FALSE, right=TRUE)
  }
  invisible(NULL)
}

setMethod("show", "DataFrame", .show_DF)

.make_rownames_for_display <- function(x_rownames, nrow, nhead, ntail) {
  p1 <- ifelse (nhead == 0L, 0L, 1L)
  p2 <- ifelse (ntail == 0L, 0L, ntail - 1L)
  s1 <- s2 <- character(0)
  if (is.null(x_rownames)) {
    if (nhead > 0)
      s1 <- paste0(as.character(p1:nhead))
    if (ntail > 0)
      s2 <- paste0(as.character((nrow-p2):nrow))
  } else {
    if (nhead > 0)
      s1 <- paste0(head(x_rownames, nhead))
    if (ntail > 0)
      s2 <- paste0(tail(x_rownames, ntail))
  }
  c(s1, "...", s2)
}

.make_class_info_for_display <- function(x) {
  matrix(
    unlist(
      lapply(x, function(col) paste0("<", .classNameForDisplay(col), ">")),
      use.names=FALSE
    ),
    nrow=1L,
    dimnames=list("", colnames(x))
  )
}

.classNameForDisplay <- function(x) class(x)[1L]
