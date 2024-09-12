test_that("filter works with regular columns", {
    d <- S4Vectors::DataFrame(mtcars)

    f <- filter(d, am == 0)
    expect_s4_class(f, "DataFrame")
    expect_identical(nrow(f), 19L)
    expect_identical(names(f), names(mtcars))
    expect_true(all(f$am == 0))

    f <- filter(d, am == 0, gear == 3)
    expect_identical(nrow(f), 15L)
    expect_identical(names(f), names(mtcars))
    expect_true(all(f$am == 0))
    expect_true(all(f$gear == 3))

    f <- filter(d, mpg > 15)
    expect_identical(nrow(f), 26L)
    expect_identical(names(f), names(mtcars))
    expect_true(all(f$mpg > 15))
})


test_that("filter works with S4 columns", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")
    skip_if_not_installed("BiocGenerics")

    D <- S4Vectors::DataFrame(mtcars)
    D$gr <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))

    f <- filter(D, BiocGenerics::end(gr) < 20)
    expect_s4_class(f, "DataFrame")
    expect_identical(nrow(f), 10L)
    expect_identical(names(f), names(D))
    expect_true(all(BiocGenerics::end(f$gr) < 20))

    expect_s4_class(f$gr, "GRanges")
    expect_identical(f$gr, D$gr[1:10])
})
