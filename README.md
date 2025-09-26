
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DFplyr

<!-- badges: start -->

<!-- badges: end -->

The goal of DFplyr is to enable `dplyr` and `ggplot2` support for
`S4Vectors::DataFrame` by providing the appropriate extension methods.
As row names are an important feature of many Bioconductor structures,
these are preserved where possible.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jonocarroll/DFplyr")
```

You can install from [Bioconductor](https://bioconductor.org) with:

``` r
if (!require("BiocManager", quietly =TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("DFplyr")
```

## Examples

First create an S4Vectors `DataFrame`, including S4 columns if desired

``` r
suppressMessages(library(S4Vectors))
m <- mtcars[, c("cyl", "hp", "am", "gear", "disp")]
d <- as(m, "DataFrame")
d$grX <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))
d$grY <- GenomicRanges::GRanges("chrY", IRanges::IRanges(1:32, width = 10))
d$nl <- IRanges::NumericList(lapply(d$gear, function(n) round(rnorm(n), 2)))
d
#> DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                 6       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4       160  chrX:2-11
#> Datsun 710                4        93         1         4       108  chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258  chrX:4-13
#> Hornet Sportabout         8       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa              4       113         1         5      95.1 chrX:28-37
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#> Ferrari Dino              6       175         1         5     145.0 chrX:30-39
#> Maserati Bora             8       335         1         5     301.0 chrX:31-40
#> Volvo 142E                4       109         1         4     121.0 chrX:32-41
#>                          grY                    nl
#>                    <GRanges>         <NumericList>
#> Mazda RX4          chrY:1-10 -1.30, 0.23, 2.60,...
#> Mazda RX4 Wag      chrY:2-11 -0.22,-0.17,-0.23,...
#> Datsun 710         chrY:3-12 -0.03,-0.74,-0.93,...
#> Hornet 4 Drive     chrY:4-13     -0.87,-2.53,-0.26
#> Hornet Sportabout  chrY:5-14     -0.48,-2.95,-0.17
#> ...                      ...                   ...
#> Lotus Europa      chrY:28-37 -0.80,-1.10,-0.41,...
#> Ford Pantera L    chrY:29-38 -0.61,-1.41, 1.43,...
#> Ferrari Dino      chrY:30-39  0.24, 0.49,-0.11,...
#> Maserati Bora     chrY:31-40 -0.56,-0.91,-1.02,...
#> Volvo 142E        chrY:32-41  1.86,-1.82,-1.15,...
```

This will appear in RStudio’s environment pane as a

    Formal class DataFrame (dplyr-compatible)

when using `DFplyr`. No interference with the actual object is required,
but this helps identify that `dplyr`-compatibility is available.

`DataFrame`s can then be used in `dplyr` calls the same as `data.frame`
or `tibble` objects. Support for working with S4 columns is enabled
provided they have appropriate functions. Adding multiple columns will
result in the new columns being created in alphabetical order

``` r
suppressMessages(library(DFplyr))

mutate(d, newvar = cyl + hp)
#> DataFrame with 32 rows and 9 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                 6       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4       160  chrX:2-11
#> Datsun 710                4        93         1         4       108  chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258  chrX:4-13
#> Hornet Sportabout         8       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa              4       113         1         5      95.1 chrX:28-37
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#> Ferrari Dino              6       175         1         5     145.0 chrX:30-39
#> Maserati Bora             8       335         1         5     301.0 chrX:31-40
#> Volvo 142E                4       109         1         4     121.0 chrX:32-41
#>                          grY                      nl    newvar
#>                    <GRanges> <CompressedNumericList> <numeric>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...       116
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...       116
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...        97
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26       116
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17       183
#> ...                      ...                     ...       ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...       117
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...       272
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...       181
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...       343
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...       113

mutate(d, nl2 = nl * 2)
#> DataFrame with 32 rows and 9 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                 6       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4       160  chrX:2-11
#> Datsun 710                4        93         1         4       108  chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258  chrX:4-13
#> Hornet Sportabout         8       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa              4       113         1         5      95.1 chrX:28-37
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#> Ferrari Dino              6       175         1         5     145.0 chrX:30-39
#> Maserati Bora             8       335         1         5     301.0 chrX:31-40
#> Volvo 142E                4       109         1         4     121.0 chrX:32-41
#>                          grY                      nl                     nl2
#>                    <GRanges> <CompressedNumericList> <CompressedNumericList>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...   -2.60, 0.46, 5.20,...
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...   -0.44,-0.34,-0.46,...
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...   -0.06,-1.48,-1.86,...
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26       -1.74,-5.06,-0.52
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17       -0.96,-5.90,-0.34
#> ...                      ...                     ...                     ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...   -1.60,-2.20,-0.82,...
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...   -1.22,-2.82, 2.86,...
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...    0.48, 0.98,-0.22,...
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...   -1.12,-1.82,-2.04,...
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...    3.72,-3.64,-2.30,...

mutate(d, length_nl = lengths(nl))
#> DataFrame with 32 rows and 9 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                 6       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4       160  chrX:2-11
#> Datsun 710                4        93         1         4       108  chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258  chrX:4-13
#> Hornet Sportabout         8       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa              4       113         1         5      95.1 chrX:28-37
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#> Ferrari Dino              6       175         1         5     145.0 chrX:30-39
#> Maserati Bora             8       335         1         5     301.0 chrX:31-40
#> Volvo 142E                4       109         1         4     121.0 chrX:32-41
#>                          grY                      nl length_nl
#>                    <GRanges> <CompressedNumericList> <integer>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...         4
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...         4
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...         4
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26         3
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17         3
#> ...                      ...                     ...       ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...         5
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...         5
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...         5
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...         5
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...         4

mutate(d,
    chr = GenomeInfoDb::seqnames(grX),
    strand_X = BiocGenerics::strand(grX),
    end_X = BiocGenerics::end(grX)
)
#> DataFrame with 32 rows and 11 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                 6       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4       160  chrX:2-11
#> Datsun 710                4        93         1         4       108  chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258  chrX:4-13
#> Hornet Sportabout         8       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa              4       113         1         5      95.1 chrX:28-37
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#> Ferrari Dino              6       175         1         5     145.0 chrX:30-39
#> Maserati Bora             8       335         1         5     301.0 chrX:31-40
#> Volvo 142E                4       109         1         4     121.0 chrX:32-41
#>                          grY                      nl   chr     end_X strand_X
#>                    <GRanges> <CompressedNumericList> <Rle> <integer>    <Rle>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...  chrX        10        *
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...  chrX        11        *
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...  chrX        12        *
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26  chrX        13        *
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17  chrX        14        *
#> ...                      ...                     ...   ...       ...      ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...  chrX        37        *
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...  chrX        38        *
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...  chrX        39        *
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...  chrX        40        *
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...  chrX        41        *
```

the object returned remains a standard `DataFrame`, and further calls
can be piped with `%>%`

``` r
mutate(d, newvar = cyl + hp) %>%
    pull(newvar)
#>  [1] 116 116  97 116 183 111 253  66  99 129 129 188 188 188 213 223 238  70  56
#> [20]  69 101 158 158 253 183  70  95 117 272 181 343 113
```

Some of the variants of the `dplyr` verbs also work

``` r
mutate_if(d, is.numeric, ~ .^2)
#> DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                36     12100         1        16     25600  chrX:1-10
#> Mazda RX4 Wag            36     12100         1        16     25600  chrX:2-11
#> Datsun 710               16      8649         1        16     11664  chrX:3-12
#> Hornet 4 Drive           36     12100         0         9     66564  chrX:4-13
#> Hornet Sportabout        64     30625         0         9    129600  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa             16     12769         1        25   9044.01 chrX:28-37
#> Ford Pantera L           64     69696         1        25 123201.00 chrX:29-38
#> Ferrari Dino             36     30625         1        25  21025.00 chrX:30-39
#> Maserati Bora            64    112225         1        25  90601.00 chrX:31-40
#> Volvo 142E               16     11881         1        16  14641.00 chrX:32-41
#>                          grY                      nl
#>                    <GRanges> <CompressedNumericList>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17
#> ...                      ...                     ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...

mutate_if(d, ~ inherits(., "GRanges"), BiocGenerics::start)
#> DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp       grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric> <integer>
#> Mazda RX4                 6       110         1         4       160         1
#> Mazda RX4 Wag             6       110         1         4       160         2
#> Datsun 710                4        93         1         4       108         3
#> Hornet 4 Drive            6       110         0         3       258         4
#> Hornet Sportabout         8       175         0         3       360         5
#> ...                     ...       ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1        28
#> Ford Pantera L            8       264         1         5     351.0        29
#> Ferrari Dino              6       175         1         5     145.0        30
#> Maserati Bora             8       335         1         5     301.0        31
#> Volvo 142E                4       109         1         4     121.0        32
#>                         grY                      nl
#>                   <integer> <CompressedNumericList>
#> Mazda RX4                 1   -1.30, 0.23, 2.60,...
#> Mazda RX4 Wag             2   -0.22,-0.17,-0.23,...
#> Datsun 710                3   -0.03,-0.74,-0.93,...
#> Hornet 4 Drive            4       -0.87,-2.53,-0.26
#> Hornet Sportabout         5       -0.48,-2.95,-0.17
#> ...                     ...                     ...
#> Lotus Europa             28   -0.80,-1.10,-0.41,...
#> Ford Pantera L           29   -0.61,-1.41, 1.43,...
#> Ferrari Dino             30    0.24, 0.49,-0.11,...
#> Maserati Bora            31   -0.56,-0.91,-1.02,...
#> Volvo 142E               32    1.86,-1.82,-1.15,...
```

Use of `tidyselect` helpers is limited to within `dplyr::vars()` calls
and using the `_at` variants

``` r
mutate_at(d, vars(starts_with("c")), ~ .^2)
#> DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                36       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag            36       110         1         4       160  chrX:2-11
#> Datsun 710               16        93         1         4       108  chrX:3-12
#> Hornet 4 Drive           36       110         0         3       258  chrX:4-13
#> Hornet Sportabout        64       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa             16       113         1         5      95.1 chrX:28-37
#> Ford Pantera L           64       264         1         5     351.0 chrX:29-38
#> Ferrari Dino             36       175         1         5     145.0 chrX:30-39
#> Maserati Bora            64       335         1         5     301.0 chrX:31-40
#> Volvo 142E               16       109         1         4     121.0 chrX:32-41
#>                          grY                      nl
#>                    <GRanges> <CompressedNumericList>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17
#> ...                      ...                     ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...

select_at(d, vars(starts_with("gr")))
#> DataFrame with 32 rows and 2 columns
#>                          grX        grY
#>                    <GRanges>  <GRanges>
#> Mazda RX4          chrX:1-10  chrY:1-10
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11
#> Datsun 710         chrX:3-12  chrY:3-12
#> Hornet 4 Drive     chrX:4-13  chrY:4-13
#> Hornet Sportabout  chrX:5-14  chrY:5-14
#> ...                      ...        ...
#> Lotus Europa      chrX:28-37 chrY:28-37
#> Ford Pantera L    chrX:29-38 chrY:29-38
#> Ferrari Dino      chrX:30-39 chrY:30-39
#> Maserati Bora     chrX:31-40 chrY:31-40
#> Volvo 142E        chrX:32-41 chrY:32-41
```

Importantly, grouped operations are supported. `DataFrame` does not
natively support groups (the same way that `data.frame` does not) so
these are implemented specifically for `DFplyr`

``` r
group_by(d, cyl, am)
#> DataFrame with 32 rows and 8 columns
#> Groups:  cyl, am 
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                 6       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4       160  chrX:2-11
#> Datsun 710                4        93         1         4       108  chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258  chrX:4-13
#> Hornet Sportabout         8       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa              4       113         1         5      95.1 chrX:28-37
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#> Ferrari Dino              6       175         1         5     145.0 chrX:30-39
#> Maserati Bora             8       335         1         5     301.0 chrX:31-40
#> Volvo 142E                4       109         1         4     121.0 chrX:32-41
#>                          grY                      nl
#>                    <GRanges> <CompressedNumericList>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17
#> ...                      ...                     ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...
```

Other verbs are similarly implemented, and preserve row names where
possible

``` r
select(d, am, cyl)
#> DataFrame with 32 rows and 2 columns
#>                          am       cyl
#>                   <numeric> <numeric>
#> Mazda RX4                 1         6
#> Mazda RX4 Wag             1         6
#> Datsun 710                1         4
#> Hornet 4 Drive            0         6
#> Hornet Sportabout         0         8
#> ...                     ...       ...
#> Lotus Europa              1         4
#> Ford Pantera L            1         8
#> Ferrari Dino              1         6
#> Maserati Bora             1         8
#> Volvo 142E                1         4

arrange(d, desc(hp))
#> DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Maserati Bora             8       335         1         5       301 chrX:31-40
#> Ford Pantera L            8       264         1         5       351 chrX:29-38
#> Duster 360                8       245         0         3       360  chrX:7-16
#> Camaro Z28                8       245         0         3       350 chrX:24-33
#> Chrysler Imperial         8       230         0         3       440 chrX:17-26
#> ...                     ...       ...       ...       ...       ...        ...
#> Fiat 128                  4        66         1         4      78.7 chrX:18-27
#> Fiat X1-9                 4        66         1         4      79.0 chrX:26-35
#> Toyota Corolla            4        65         1         4      71.1 chrX:20-29
#> Merc 240D                 4        62         0         4     146.7  chrX:8-17
#> Honda Civic               4        52         1         4      75.7 chrX:19-28
#>                          grY                      nl
#>                    <GRanges> <CompressedNumericList>
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...
#> Duster 360         chrY:7-16        0.67,-0.97,-0.76
#> Camaro Z28        chrY:24-33       -0.87, 1.31,-1.45
#> Chrysler Imperial chrY:17-26        0.38, 1.22,-1.25
#> ...                      ...                     ...
#> Fiat 128          chrY:18-27    1.58,-0.43,-0.97,...
#> Fiat X1-9         chrY:26-35   -1.60,-0.76, 0.39,...
#> Toyota Corolla    chrY:20-29    0.84,-1.32, 2.19,...
#> Merc 240D          chrY:8-17   -0.64, 0.53, 0.49,...
#> Honda Civic       chrY:19-28      0.38,0.24,0.01,...

filter(d, am == 0)
#> DataFrame with 19 rows and 8 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Hornet 4 Drive            6       110         0         3     258.0  chrX:4-13
#> Hornet Sportabout         8       175         0         3     360.0  chrX:5-14
#> Valiant                   6       105         0         3     225.0  chrX:6-15
#> Duster 360                8       245         0         3     360.0  chrX:7-16
#> Merc 240D                 4        62         0         4     146.7  chrX:8-17
#> ...                     ...       ...       ...       ...       ...        ...
#> Toyota Corona             4        97         0         3     120.1 chrX:21-30
#> Dodge Challenger          8       150         0         3     318.0 chrX:22-31
#> AMC Javelin               8       150         0         3     304.0 chrX:23-32
#> Camaro Z28                8       245         0         3     350.0 chrX:24-33
#> Pontiac Firebird          8       175         0         3     400.0 chrX:25-34
#>                          grY                      nl
#>                    <GRanges> <CompressedNumericList>
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17
#> Valiant            chrY:6-15       -0.97, 0.10, 0.76
#> Duster 360         chrY:7-16        0.67,-0.97,-0.76
#> Merc 240D          chrY:8-17   -0.64, 0.53, 0.49,...
#> ...                      ...                     ...
#> Toyota Corona     chrY:21-30          1.74,0.63,0.03
#> Dodge Challenger  chrY:22-31       -1.39, 0.51,-0.90
#> AMC Javelin       chrY:23-32        0.37,-0.96, 1.56
#> Camaro Z28        chrY:24-33       -0.87, 1.31,-1.45
#> Pontiac Firebird  chrY:25-34          1.23,1.43,0.95

slice(d, 3:6)
#> DataFrame with 4 rows and 8 columns
#>                         cyl        hp        am      gear      disp       grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric> <GRanges>
#> Datsun 710                4        93         1         4       108 chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258 chrX:4-13
#> Hornet Sportabout         8       175         0         3       360 chrX:5-14
#> Valiant                   6       105         0         3       225 chrX:6-15
#>                         grY                      nl
#>                   <GRanges> <CompressedNumericList>
#> Datsun 710        chrY:3-12   -0.03,-0.74,-0.93,...
#> Hornet 4 Drive    chrY:4-13       -0.87,-2.53,-0.26
#> Hornet Sportabout chrY:5-14       -0.48,-2.95,-0.17
#> Valiant           chrY:6-15       -0.97, 0.10, 0.76

group_by(d, gear) %>%
    slice(1:2)
#> DataFrame with 6 rows and 8 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Hornet Sportabout         8       175         0         3     360.0  chrX:5-14
#> Merc 450SL                8       180         0         3     275.8 chrX:13-22
#> Mazda RX4                 6       110         1         4     160.0  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4     160.0  chrX:2-11
#> Porsche 914-2             4        91         1         5     120.3 chrX:27-36
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#>                          grY                      nl
#>                    <GRanges> <CompressedNumericList>
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17
#> Merc 450SL        chrY:13-22          0.63,1.48,0.55
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...
#> Porsche 914-2     chrY:27-36    0.47,-1.47, 0.25,...
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...
```

`rename` is itself renamed to `rename2` due to conflicts between {dplyr}
and {S4Vectors}, but works in the {dplyr} sense of taking `new = old`
replacements with NSE syntax

``` r
select(d, am, cyl) %>%
    rename2(foo = am)
#> Warning in rename2(., foo = am): DFplyr now properly supports rename with NSE
#> syntax
```

Row names are not preserved when there may be duplicates or they don’t
make sense, otherwise the first label (according to the current
de-duplication method, in the case of `distinct`, this is via
`BiocGenerics::duplicated`). This may have complications for S4 columns.

``` r
distinct(d)
#> DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp        grX
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>  <GRanges>
#> Mazda RX4                 6       110         1         4       160  chrX:1-10
#> Mazda RX4 Wag             6       110         1         4       160  chrX:2-11
#> Datsun 710                4        93         1         4       108  chrX:3-12
#> Hornet 4 Drive            6       110         0         3       258  chrX:4-13
#> Hornet Sportabout         8       175         0         3       360  chrX:5-14
#> ...                     ...       ...       ...       ...       ...        ...
#> Lotus Europa              4       113         1         5      95.1 chrX:28-37
#> Ford Pantera L            8       264         1         5     351.0 chrX:29-38
#> Ferrari Dino              6       175         1         5     145.0 chrX:30-39
#> Maserati Bora             8       335         1         5     301.0 chrX:31-40
#> Volvo 142E                4       109         1         4     121.0 chrX:32-41
#>                          grY                      nl
#>                    <GRanges> <CompressedNumericList>
#> Mazda RX4          chrY:1-10   -1.30, 0.23, 2.60,...
#> Mazda RX4 Wag      chrY:2-11   -0.22,-0.17,-0.23,...
#> Datsun 710         chrY:3-12   -0.03,-0.74,-0.93,...
#> Hornet 4 Drive     chrY:4-13       -0.87,-2.53,-0.26
#> Hornet Sportabout  chrY:5-14       -0.48,-2.95,-0.17
#> ...                      ...                     ...
#> Lotus Europa      chrY:28-37   -0.80,-1.10,-0.41,...
#> Ford Pantera L    chrY:29-38   -0.61,-1.41, 1.43,...
#> Ferrari Dino      chrY:30-39    0.24, 0.49,-0.11,...
#> Maserati Bora     chrY:31-40   -0.56,-0.91,-1.02,...
#> Volvo 142E        chrY:32-41    1.86,-1.82,-1.15,...

group_by(d, cyl, am) %>%
    tally(gear)
#> DataFrame with 6 rows and 3 columns
#>         cyl        am         n
#>   <numeric> <numeric> <numeric>
#> 1         4         0        11
#> 2         4         1        34
#> 3         6         0        14
#> 4         6         1        13
#> 5         8         0        36
#> 6         8         1        10

count(d, gear, am, cyl)
#> DataFrame with 10 rows and 4 columns
#>        gear    am   cyl         n
#>    <factor> <Rle> <Rle> <integer>
#> 1         3     0     4         1
#> 2         3     0     6         2
#> 3         3     0     8        12
#> 4         4     0     4         2
#> 5         4     0     6         2
#> 6         4     1     4         6
#> 7         4     1     6         2
#> 8         5     1     4         2
#> 9         5     1     6         1
#> 10        5     1     8         2
```

## Joins

Joins attempt to preserve rownames and grouping wherever possible

``` r
Da <- as(starwars[, c("name", "eye_color", "height", "mass")], "DataFrame") |> 
  head(10) |> 
  group_by(eye_color)
Da
#> DataFrame with 10 rows and 4 columns
#> Groups:  eye_color 
#>                  name   eye_color    height      mass
#>           <character> <character> <integer> <numeric>
#> 1      Luke Skywalker        blue       172        77
#> 2               C-3PO      yellow       167        75
#> 3               R2-D2         red        96        32
#> 4         Darth Vader      yellow       202       136
#> 5         Leia Organa       brown       150        49
#> 6           Owen Lars        blue       178       120
#> 7  Beru Whitesun Lars        blue       165        75
#> 8               R5-D4         red        97        32
#> 9   Biggs Darklighter       brown       183        84
#> 10     Obi-Wan Kenobi   blue-gray       182        77

Db <- as(starwars[, c("name", "eye_color", "homeworld")], "DataFrame")
Db
#> DataFrame with 87 rows and 3 columns
#>               name   eye_color   homeworld
#>        <character> <character> <character>
#> 1   Luke Skywalker        blue    Tatooine
#> 2            C-3PO      yellow    Tatooine
#> 3            R2-D2         red       Naboo
#> 4      Darth Vader      yellow    Tatooine
#> 5      Leia Organa       brown    Alderaan
#> ...            ...         ...         ...
#> 83            Finn        dark          NA
#> 84             Rey       hazel          NA
#> 85     Poe Dameron       brown          NA
#> 86             BB8       black          NA
#> 87  Captain Phasma     unknown          NA

left_join(Da, Db)
#> Joining with `by = c("name", "eye_color")`
#> DataFrame with 10 rows and 5 columns
#> Groups:  eye_color 
#>                  name   eye_color    height      mass   homeworld
#>           <character> <character> <integer> <numeric> <character>
#> 1      Luke Skywalker        blue       172        77    Tatooine
#> 2               C-3PO      yellow       167        75    Tatooine
#> 3               R2-D2         red        96        32       Naboo
#> 4         Darth Vader      yellow       202       136    Tatooine
#> 5         Leia Organa       brown       150        49    Alderaan
#> 6           Owen Lars        blue       178       120    Tatooine
#> 7  Beru Whitesun Lars        blue       165        75    Tatooine
#> 8               R5-D4         red        97        32    Tatooine
#> 9   Biggs Darklighter       brown       183        84    Tatooine
#> 10     Obi-Wan Kenobi   blue-gray       182        77     Stewjon

right_join(Da, Db)
#> Joining with `by = c("name", "eye_color")`
#> DataFrame with 87 rows and 5 columns
#> Groups:  eye_color 
#>               name     eye_color    height      mass   homeworld
#>        <character>   <character> <integer> <numeric> <character>
#> 1   Luke Skywalker          blue       172        77    Tatooine
#> 2            C-3PO        yellow       167        75    Tatooine
#> 3            R2-D2           red        96        32       Naboo
#> 4      Darth Vader        yellow       202       136    Tatooine
#> 5      Leia Organa         brown       150        49    Alderaan
#> ...            ...           ...       ...       ...         ...
#> 83             BB8         black        NA        NA          NA
#> 84  Captain Phasma       unknown        NA        NA          NA
#> 85        San Hill          gold        NA        NA  Muunilinst
#> 86        Shaak Ti         black        NA        NA       Shili
#> 87        Grievous green, yellow        NA        NA       Kalee

inner_join(Da, Db[1:3, ])
#> Joining with `by = c("name", "eye_color")`
#> DataFrame with 3 rows and 5 columns
#> Groups:  eye_color 
#>             name   eye_color    height      mass   homeworld
#>      <character> <character> <integer> <numeric> <character>
#> 1 Luke Skywalker        blue       172        77    Tatooine
#> 2          C-3PO      yellow       167        75    Tatooine
#> 3          R2-D2         red        96        32       Naboo

full_join(Da, Db[1:3, ])
#> Joining with `by = c("name", "eye_color")`
#> DataFrame with 10 rows and 5 columns
#> Groups:  eye_color 
#>                  name   eye_color    height      mass   homeworld
#>           <character> <character> <integer> <numeric> <character>
#> 1      Luke Skywalker        blue       172        77    Tatooine
#> 2               C-3PO      yellow       167        75    Tatooine
#> 3               R2-D2         red        96        32       Naboo
#> 4         Leia Organa       brown       150        49          NA
#> 5           Owen Lars        blue       178       120          NA
#> 6  Beru Whitesun Lars        blue       165        75          NA
#> 7         Darth Vader      yellow       202       136          NA
#> 8   Biggs Darklighter       brown       183        84          NA
#> 9      Obi-Wan Kenobi   blue-gray       182        77          NA
#> 10              R5-D4         red        97        32          NA
```

## Coverage

Most `dplyr` functions are implemented.

If you find any which are not, please [file an
issue](https://github.com/jonocarroll/DFplyr/issues/new).
