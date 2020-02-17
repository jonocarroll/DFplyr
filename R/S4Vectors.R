
.show_DF <- function(x) {
  nhead <- S4Vectors::get_showHeadLines()
  ntail <- S4Vectors::get_showTailLines()
  x_nrow <- nrow(x)
  x_ncol <- ncol(x)
  cat(S4Vectors::classNameForDisplay(x), " with ",
      x_nrow, " row", ifelse(x_nrow == 1L, "", "s"),
      " and ",
      x_ncol, " column", ifelse(x_ncol == 1L, "", "s"),
      "\n", sep="")
  if (!is.null(group_data(x)) & nrow(group_data(x)) > 1L) {
    cat("Groups: ", toString(group_vars(x)), "\n")
  }
  if (x_nrow != 0L && x_ncol != 0L) {
    x_rownames <- rownames(x)
    if (x_nrow <= nhead + ntail + 1L) {
      m <- S4Vectors:::makeNakedCharacterMatrixForDisplay(x)
      if (!is.null(x_rownames))
        rownames(m) <- x_rownames
    } else {
      m <- rbind(S4Vectors:::makeNakedCharacterMatrixForDisplay(head(x, nhead)),
                 rbind(rep.int("...", x_ncol)),
                 S4Vectors:::makeNakedCharacterMatrixForDisplay(tail(x, ntail)))
      rownames(m) <- S4Vectors:::make_rownames_for_DataTable_display(
        x_rownames, x_nrow,
        nhead, ntail)
    }
    m <- rbind(S4Vectors:::make_class_info_for_DataTable_display(x), m)
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
      lapply(x, function(col) paste0("<", S4Vectors::.classNameForDisplay(col), ">")),
      use.names=FALSE
    ),
    nrow=1L,
    dimnames=list("", colnames(x))
  )
}

.classNameForDisplay <- function(x) class(x)[1L]
