
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

## Examples

Most `dplyr` functions are implemented. If you find any which are not,
please [file an
issue](https://github.com/jonocarroll/DFplyr/issues/new).

``` r
suppressPackageStartupMessages(
  suppressWarnings({
    library(S4Vectors)
    library(dplyr)
    library(DFplyr)
  }))
```

First create an S4Vectors DataFrame, including S4 columns if desired

``` r
library(S4Vectors)
m <- mtcars[, c("cyl", "hp", "am", "gear", "disp")]
d <- as(m, "DataFrame")
d$grX <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width=10))
d$grY <- GenomicRanges::GRanges("chrY", IRanges::IRanges(1:32, width = 10))
d$nl <- IRanges::NumericList(lapply(d$gear, function(n) round(rnorm(n), 2)))
d
#> dplyr-compatible DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                          grX        grY                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36
#> ...                      ...        ...                   ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...
```

This will appear in RStudio’s environment pane as a `Formal class
DataFrame (dplyr-compatible)` when using `DFplyr`. No interference with
the actual object is required, but this helps identify that
`dplyr`-compatibility is available.

`DataFrame`s can then be used in `dplyr` calls the same as `data.frame`
or `tibble` objects. Support for working with S4 columns is enabled
provided they have appropriate functions. Adding multiple columns will
result in the new columns being created in alphabetical order

``` r
mutate(d, newvar = cyl + hp)
#> dplyr-compatible DataFrame with 32 rows and 9 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                          grX        grY                    nl    newvar
#>                    <GRanges>  <GRanges>         <NumericList> <numeric>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...       116
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...       116
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...        97
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05       116
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36       183
#> ...                      ...        ...                   ...       ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...       117
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...       272
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...       181
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...       343
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...       113

mutate(d, nl2 = nl * 2)
#> dplyr-compatible DataFrame with 32 rows and 9 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                          grX        grY                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36
#> ...                      ...        ...                   ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...
#>                                     nl2
#>                           <NumericList>
#> Mazda RX4         -3.28,-4.08,-2.38,...
#> Mazda RX4 Wag        2.28,0.38,0.56,...
#> Datsun 710          0.06,0.62,-1.84,...
#> Hornet 4 Drive           2.58,4.56,-4.1
#> Hornet Sportabout        -2.52,1.9,0.72
#> ...                                 ...
#> Lotus Europa        -1.78,2.22,1.18,...
#> Ford Pantera L      2.68,0.04,-1.32,...
#> Ferrari Dino       -1.26,-0.2,-0.32,...
#> Maserati Bora       4.06,-1.38,1.42,...
#> Volvo 142E           0.26,-0.26,0.1,...

mutate(d, length_nl = lengths(nl))
#> dplyr-compatible DataFrame with 32 rows and 9 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                          grX        grY                    nl length_nl
#>                    <GRanges>  <GRanges>         <NumericList> <integer>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...         4
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...         4
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...         4
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05         3
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36         3
#> ...                      ...        ...                   ...       ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...         5
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...         5
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...         5
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...         5
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...         4

mutate(d, 
       chr = GenomeInfoDb::seqnames(grX), 
       strand_X = BiocGenerics::strand(grX), 
       end_X = BiocGenerics::end(grX))
#> dplyr-compatible DataFrame with 32 rows and 11 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                          grX        grY                    nl   chr
#>                    <GRanges>  <GRanges>         <NumericList> <Rle>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...  chrX
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...  chrX
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...  chrX
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05  chrX
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36  chrX
#> ...                      ...        ...                   ...   ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...  chrX
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...  chrX
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...  chrX
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...  chrX
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...  chrX
#>                       end_X strand_X
#>                   <integer>    <Rle>
#> Mazda RX4                10        *
#> Mazda RX4 Wag            11        *
#> Datsun 710               12        *
#> Hornet 4 Drive           13        *
#> Hornet Sportabout        14        *
#> ...                     ...      ...
#> Lotus Europa             37        *
#> Ford Pantera L           38        *
#> Ferrari Dino             39        *
#> Maserati Bora            40        *
#> Volvo 142E               41        *
```

the object returned remains a standard `DataFrame`, and further calls
can be piped with `%>%`

``` r
mutate(d, newvar = cyl + hp) %>%
  pull(newvar)
#>  [1] 116 116  97 116 183 111 253  66  99 129 129 188 188 188 213 223 238
#> [18]  70  56  69 101 158 158 253 183  70  95 117 272 181 343 113
```

Some of the variants of the `dplyr` verbs also work

``` r
mutate_if(d, is.numeric, ~.^2)
#> dplyr-compatible DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                36     12100         1        16     25600
#> Mazda RX4 Wag            36     12100         1        16     25600
#> Datsun 710               16      8649         1        16     11664
#> Hornet 4 Drive           36     12100         0         9     66564
#> Hornet Sportabout        64     30625         0         9    129600
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa             16     12769         1        25   9044.01
#> Ford Pantera L           64     69696         1        25    123201
#> Ferrari Dino             36     30625         1        25     21025
#> Maserati Bora            64    112225         1        25     90601
#> Volvo 142E               16     11881         1        16     14641
#>                          grX        grY                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36
#> ...                      ...        ...                   ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...

mutate_if(d, ~inherits(., "GRanges"), BiocGenerics::start)
#> dplyr-compatible DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                         grX       grY                    nl
#>                   <integer> <integer>         <NumericList>
#> Mazda RX4                 1         1 -1.64,-2.04,-1.19,...
#> Mazda RX4 Wag             2         2    1.14,0.19,0.28,...
#> Datsun 710                3         3   0.03,0.31,-0.92,...
#> Hornet 4 Drive            4         4       1.29,2.28,-2.05
#> Hornet Sportabout         5         5       -1.26,0.95,0.36
#> ...                     ...       ...                   ...
#> Lotus Europa             28        28   -0.89,1.11,0.59,...
#> Ford Pantera L           29        29   1.34,0.02,-0.66,...
#> Ferrari Dino             30        30  -0.63,-0.1,-0.16,...
#> Maserati Bora            31        31   2.03,-0.69,0.71,...
#> Volvo 142E               32        32   0.13,-0.13,0.05,...
```

Use of `tidyselect` helpers is limited to within `dplyr::vars()` calls
and using the `_at` variants

``` r
mutate_at(d, vars(starts_with("c")), ~.^2)
#> dplyr-compatible DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                36       110         1         4       160
#> Mazda RX4 Wag            36       110         1         4       160
#> Datsun 710               16        93         1         4       108
#> Hornet 4 Drive           36       110         0         3       258
#> Hornet Sportabout        64       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa             16       113         1         5      95.1
#> Ford Pantera L           64       264         1         5       351
#> Ferrari Dino             36       175         1         5       145
#> Maserati Bora            64       335         1         5       301
#> Volvo 142E               16       109         1         4       121
#>                          grX        grY                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36
#> ...                      ...        ...                   ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...

select_at(d, vars(starts_with("gr")))
#> dplyr-compatible DataFrame with 32 rows and 2 columns
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
#> dplyr-compatible DataFrame with 32 rows and 8 columns
#> Groups:  cyl, am 
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                          grX        grY                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36
#> ...                      ...        ...                   ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...

group_by(d, cyl) %>% 
  top_n(1, disp)
#> dplyr-compatible DataFrame with 3 rows and 8 columns
#>                          cyl        hp        am      gear      disp
#>                    <numeric> <numeric> <numeric> <numeric> <numeric>
#> Merc 240D                  4        62         0         4     146.7
#> Hornet 4 Drive             6       110         0         3       258
#> Cadillac Fleetwood         8       205         0         3       472
#>                           grX        grY                   nl
#>                     <GRanges>  <GRanges>        <NumericList>
#> Merc 240D           chrX:8-17  chrY:8-17 -1.55,-0.57,0.19,...
#> Hornet 4 Drive      chrX:4-13  chrY:4-13      1.29,2.28,-2.05
#> Cadillac Fleetwood chrX:15-24 chrY:15-24    -0.81,-1.12,-0.67
```

Other verbs are similiarly implemented, and preserve row names where
possible

``` r
select(d, am, cyl)
#> dplyr-compatible DataFrame with 32 rows and 2 columns
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

select(d, am, cyl) %>%
  rename(foo = am)
#> dplyr-compatible DataFrame with 32 rows and 2 columns
#>                         foo       cyl
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
#> dplyr-compatible DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Maserati Bora             8       335         1         5       301
#> Ford Pantera L            8       264         1         5       351
#> Duster 360                8       245         0         3       360
#> Camaro Z28                8       245         0         3       350
#> Chrysler Imperial         8       230         0         3       440
#> ...                     ...       ...       ...       ...       ...
#> Fiat 128                  4        66         1         4      78.7
#> Fiat X1-9                 4        66         1         4        79
#> Toyota Corolla            4        65         1         4      71.1
#> Merc 240D                 4        62         0         4     146.7
#> Honda Civic               4        52         1         4      75.7
#>                          grX        grY                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Maserati Bora     chrX:31-40 chrY:31-40  2.03,-0.69,0.71,...
#> Ford Pantera L    chrX:29-38 chrY:29-38  1.34,0.02,-0.66,...
#> Duster 360         chrX:7-16  chrY:7-16     -0.26,-0.09,0.15
#> Camaro Z28        chrX:24-33 chrY:24-33      -0.32,1.28,0.14
#> Chrysler Imperial chrX:17-26 chrY:17-26     -0.35,1.85,-0.56
#> ...                      ...        ...                  ...
#> Fiat 128          chrX:18-27 chrY:18-27 -0.09,-0.6,-2.08,...
#> Fiat X1-9         chrX:26-35 chrY:26-35  0.57,-1.87,0.14,...
#> Toyota Corolla    chrX:20-29 chrY:20-29  -1.75,0.48,0.12,...
#> Merc 240D          chrX:8-17  chrY:8-17 -1.55,-0.57,0.19,...
#> Honda Civic       chrX:19-28 chrY:19-28  -0.37,-0.02,1.7,...

filter(d, am == 0) 
#> dplyr-compatible DataFrame with 19 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> Valiant                   6       105         0         3       225
#> Duster 360                8       245         0         3       360
#> Merc 240D                 4        62         0         4     146.7
#> ...                     ...       ...       ...       ...       ...
#> Toyota Corona             4        97         0         3     120.1
#> Dodge Challenger          8       150         0         3       318
#> AMC Javelin               8       150         0         3       304
#> Camaro Z28                8       245         0         3       350
#> Pontiac Firebird          8       175         0         3       400
#>                          grX        grY                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Hornet 4 Drive     chrX:4-13  chrY:4-13      1.29,2.28,-2.05
#> Hornet Sportabout  chrX:5-14  chrY:5-14      -1.26,0.95,0.36
#> Valiant            chrX:6-15  chrY:6-15       0.89,-0.48,0.7
#> Duster 360         chrX:7-16  chrY:7-16     -0.26,-0.09,0.15
#> Merc 240D          chrX:8-17  chrY:8-17 -1.55,-0.57,0.19,...
#> ...                      ...        ...                  ...
#> Toyota Corona     chrX:21-30 chrY:21-30     -0.57,-1.7,-1.58
#> Dodge Challenger  chrX:22-31 chrY:22-31      1.07,2.18,-0.97
#> AMC Javelin       chrX:23-32 chrY:23-32       -1.26,-0.6,0.5
#> Camaro Z28        chrX:24-33 chrY:24-33      -0.32,1.28,0.14
#> Pontiac Firebird  chrX:25-34 chrY:25-34      1.27,-0.04,1.09

slice(d, 3:6)
#> dplyr-compatible DataFrame with 4 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> Valiant                   6       105         0         3       225
#>                         grX       grY                  nl
#>                   <GRanges> <GRanges>       <NumericList>
#> Datsun 710        chrX:3-12 chrY:3-12 0.03,0.31,-0.92,...
#> Hornet 4 Drive    chrX:4-13 chrY:4-13     1.29,2.28,-2.05
#> Hornet Sportabout chrX:5-14 chrY:5-14     -1.26,0.95,0.36
#> Valiant           chrX:6-15 chrY:6-15      0.89,-0.48,0.7

# dd <- rbind(data.frame(d[1, ], row.names = "MyCar"), d)
# dd
```

Row names are not preserved when there may be duplicates or they don’t
make sense, otherwise the first label (according to the current
de-duplication method, in the case of `distinct`, this is via
`BiocGenerics::duplicated`). This may have complications for S4 columns.

``` r
distinct(d)
#> dplyr-compatible DataFrame with 32 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1         4       160
#> Mazda RX4 Wag             6       110         1         4       160
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4       113         1         5      95.1
#> Ford Pantera L            8       264         1         5       351
#> Ferrari Dino              6       175         1         5       145
#> Maserati Bora             8       335         1         5       301
#> Volvo 142E                4       109         1         4       121
#>                          grX        grY                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrX:1-10  chrY:1-10 -1.64,-2.04,-1.19,...
#> Mazda RX4 Wag      chrX:2-11  chrY:2-11    1.14,0.19,0.28,...
#> Datsun 710         chrX:3-12  chrY:3-12   0.03,0.31,-0.92,...
#> Hornet 4 Drive     chrX:4-13  chrY:4-13       1.29,2.28,-2.05
#> Hornet Sportabout  chrX:5-14  chrY:5-14       -1.26,0.95,0.36
#> ...                      ...        ...                   ...
#> Lotus Europa      chrX:28-37 chrY:28-37   -0.89,1.11,0.59,...
#> Ford Pantera L    chrX:29-38 chrY:29-38   1.34,0.02,-0.66,...
#> Ferrari Dino      chrX:30-39 chrY:30-39  -0.63,-0.1,-0.16,...
#> Maserati Bora     chrX:31-40 chrY:31-40   2.03,-0.69,0.71,...
#> Volvo 142E        chrX:32-41 chrY:32-41   0.13,-0.13,0.05,...

group_by(d, cyl, am) %>%
  tally(gear)
#> dplyr-compatible DataFrame with 6 rows and 3 columns
#>         cyl        am         n
#>   <numeric> <numeric> <numeric>
#> 1         4         0        11
#> 2         4         1        34
#> 3         6         0        14
#> 4         6         1        13
#> 5         8         0        36
#> 6         8         1        10

count(d, gear, am, cyl)
#> dplyr-compatible DataFrame with 18 rows and 4 columns
#>         gear    am   cyl         n
#>     <factor> <Rle> <Rle> <integer>
#> 1          3     0     4         1
#> 2          4     0     4         2
#> 3          5     0     4         0
#> 4          3     1     4         0
#> 5          4     1     4         6
#> ...      ...   ...   ...       ...
#> 14         4     0     8         0
#> 15         5     0     8         0
#> 16         3     1     8         0
#> 17         4     1     8         0
#> 18         5     1     8         2
```

`ggplot2` support is also enabled

``` r
library(ggplot2)
ggplot(d, aes(disp, cyl)) + geom_point()
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

## Implementation

Most of the `dplyr` verbs for `DataFrame`s are implmented by first
converting to `tibble`, performing the verb operation, then converting
back to `DataFrame`. Care has been taken to retain groups and row names
through these operations, but this may introduce some complications. If
you spot any, please [file an
issue](https://github.com/jonocarroll/DFplyr/issues/new).
