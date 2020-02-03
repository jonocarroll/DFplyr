
.show_DF <- function(object) {
  nhead <- S4Vectors::get_showHeadLines()
  ntail <- S4Vectors::get_showTailLines()
  nr <- nrow(object)
  nc <- ncol(object)
  cat("dplyr-compatible DataFrame", " with ",
      nr, " row", ifelse(nr == 1L, "", "s"),
      " and ",
      nc, " column", ifelse(nc == 1L, "", "s"),
      "\n", sep="")
  if (!is.null(group_data(object)) & nrow(group_data(object)) > 1L) {
    cat("Groups: ", toString(group_vars(object)), "\n")
  }
  if (nr > 0 && nc > 0) {
    nms <- rownames(object)
    if (nr <= (nhead + ntail + 1L)) {
      out <- as.matrix(format(as.data.frame(lapply(object,
                                                   showAsCell), optional = TRUE)))
      if (!is.null(nms))
        rownames(out) <- nms
    }
    else {
      out <- rbind(as.matrix(format(as.data.frame(lapply(object,
                                                         function(x) S4Vectors::showAsCell(head(x, nhead))), optional = TRUE))),
                   rbind(rep.int("...", nc)), as.matrix(format(as.data.frame(lapply(object,
                                                                                    function(x) S4Vectors::showAsCell(tail(x, ntail))), optional = TRUE))))
      rownames(out) <- S4Vectors:::.rownames(nms, nr, nhead, ntail)
    }
    classinfo <- matrix(unlist(lapply(object, function(x) {
      paste0("<", S4Vectors::classNameForDisplay(x)[1], ">")
    }), use.names = FALSE), nrow = 1, dimnames = list("",
                                                      colnames(out)))
    out <- rbind(classinfo, out)
    print(out, quote = FALSE, right = TRUE)
  }
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
