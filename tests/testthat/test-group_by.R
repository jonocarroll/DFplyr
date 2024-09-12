test_that("group_by works with regular columns", {
    d <- S4Vectors::DataFrame(mtcars)

    g <- group_by(d, am, cyl)
    expect_s4_class(g, "DataFrame")
    expect_identical(nrow(g), 32L)
    expect_identical(names(g), names(mtcars))
    expect_equal(g, d, ignore_attr = TRUE)

    expect_identical(group_vars(g), c("am", "cyl"))
    expect_identical(group_data(g)$am, c(rep(0, 3), rep(1, 3)))
    expect_identical(group_data(g)$cyl, rep(c(4, 6, 8), 2))
    expect_identical(lengths(group_data(g)$.rows), c(3L, 4L, 12L, 8L, 3L, 2L))
})

# D <- S4Vectors::DataFrame(mtcars)
# D$gr <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))
# D$gr[17:32] <- GenomicRanges::GRanges("chrY", IRanges::IRanges(1:16, width = 10))
#
# test_that("group_by works with S4 columns", {
#   g <- group_by(D, am, gr)
#   expect_s4_class(f, "DataFrame")
#   expect_identical(nrow(f), 10L)
#   expect_identical(names(f), names(D))
#   expect_true(all(BiocGenerics::end(f$gr) < 20))
#
#   expect_s4_class(f$gr, "GRanges")
#   expect_identical(f$gr, D$gr[1:10])
# })

test_that("filter respects groups", {
    d <- S4Vectors::DataFrame(mtcars)

    gf <- filter(group_by(d, am, cyl), hp == max(hp))
    expect_s4_class(gf, "DataFrame")
    expect_identical(nrow(gf), 8L)
    expect_identical(names(gf), names(mtcars))
    gfd <- dplyr::filter(dplyr::group_by(mtcars, am, cyl), hp == max(hp))
    expect_identical(sum(gf$mpg), sum(gfd$mpg))
})

test_that("mutate respects groups", {
    d <- S4Vectors::DataFrame(mtcars)
    gm <- mutate(group_by(d, am, cyl), newvar = mean(hp))
    expect_s4_class(gm, "DataFrame")
    expect_identical(nrow(gm), 32L)
    expect_identical(names(gm), c(names(mtcars), "newvar"))
    gmd <- dplyr::mutate(dplyr::group_by(mtcars, am, cyl), newvar = mean(hp))
    expect_equal(sum(gm$newvar), sum(gmd$newvar))
})

test_that("tally respects groups", {
    d <- S4Vectors::DataFrame(mtcars)

    gt <- tally(group_by(d, am, cyl))
    expect_s4_class(gt, "DataFrame")
    expect_identical(nrow(gt), 6L)
    expect_identical(names(gt), c("am", "cyl", "n"))
    gtd <- dplyr::tally(group_by(mtcars, am, cyl))
    expect_identical(gt$am, gtd$am)
    expect_identical(gt$cyl, gtd$cyl)
    expect_identical(gt$n, gtd$n)

    gt <- tally(group_by(d, am, cyl), hp)
    expect_s4_class(gt, "DataFrame")
    expect_identical(nrow(gt), 6L)
    expect_identical(names(gt), c("am", "cyl", "n"))
    gtd <- dplyr::tally(group_by(mtcars, am, cyl), hp)
    expect_identical(gt$am, gtd$am)
    expect_identical(gt$cyl, gtd$cyl)
    expect_identical(gt$n, gtd$n)
})


test_that("count respects groups", {
    d <- S4Vectors::DataFrame(mtcars)

    gc <- count(group_by(d, am, cyl))
    expect_s4_class(gc, "DataFrame")
    expect_identical(nrow(gc), 6L)
    expect_identical(names(gc), c("am", "cyl", "n"))
    gcd <- dplyr::tally(group_by(mtcars, am, cyl))
    expect_identical(gc$am, gcd$am)
    expect_identical(gc$cyl, gcd$cyl)
    expect_identical(gc$n, gcd$n)

    gc <- tally(group_by(d, am, cyl), hp)
    expect_s4_class(gc, "DataFrame")
    expect_identical(nrow(gc), 6L)
    expect_identical(names(gc), c("am", "cyl", "n"))
    gcd <- dplyr::tally(group_by(mtcars, am, cyl), hp)
    expect_identical(gc$am, gcd$am)
    expect_identical(gc$cyl, gcd$cyl)
    expect_identical(gc$n, gcd$n)
})
