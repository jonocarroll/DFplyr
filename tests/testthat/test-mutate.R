test_that("mutate works with regular columns", {
    d <- S4Vectors::DataFrame(mtcars)

    m <- mutate(d, newvar = cyl + mpg)
    expect_s4_class(m, "DataFrame")
    expect_true("newvar" %in% names(m))
    expect_identical(m$newvar, mtcars$cyl + mtcars$mpg)

    m <- mutate(d, newvar2 = cyl^2)
    expect_s4_class(m, "DataFrame")
    expect_true("newvar2" %in% names(m))
    expect_identical(m$newvar2, mtcars$cyl^2)

    m <- mutate_at(d, vars(starts_with("c")), ~ .^2)
    expect_s4_class(m, "DataFrame")
    expect_identical(names(m), names(d))
    expect_identical(m$cyl, mtcars$cyl^2)
    expect_identical(m$carb, mtcars$carb^2)
    expect_identical(m$mpg, mtcars$mpg)
})

test_that("mutate works with S4 columns", {
    D <- S4Vectors::DataFrame(mtcars)
    D$gr <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))

    m <- mutate(D, chr = GenomeInfoDb::seqnames(gr))
    expect_s4_class(m, "DataFrame")
    expect_identical(nrow(m), 32L)
    expect_true("chr" %in% names(m))
    expect_s4_class(m$chr, "Rle")
    expect_identical(m$chr, IRanges::RleList(factor(rep("chrX", 32)))[[1]])
})
