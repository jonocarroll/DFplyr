test_that("select works with regular columns", {
    d <- S4Vectors::DataFrame(mtcars)

    s <- select(d, cyl, mpg)
    expect_s4_class(s, "DataFrame")
    expect_identical(nrow(s), 32L)
    expect_identical(names(s), c("cyl", "mpg"))

    s <- select(d, disp:wt)
    expect_s4_class(s, "DataFrame")
    expect_identical(nrow(s), 32L)
    expect_identical(names(s), c("disp", "hp", "drat", "wt"))

    s <- select_at(d, vars(starts_with("c")))
    expect_s4_class(s, "DataFrame")
    expect_identical(nrow(s), 32L)
    expect_identical(names(s), c("cyl", "carb"))
})

test_that("mutate works with S4 columns", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    D <- S4Vectors::DataFrame(mtcars)
    D$gr <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))

    s <- select(D, cyl, gr)
    expect_s4_class(s, "DataFrame")
    expect_identical(nrow(s), 32L)
    expect_identical(names(s), c("cyl", "gr"))
    expect_s4_class(s$gr, "GRanges")
})
