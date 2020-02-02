
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

First create an S4Vectors DataFrame

``` r
# m <- mtcars[, c("cyl", "hp", "am", "gear", "disp")]
# d <- as(m, "DataFrame")
# d
library(S4Vectors)
m <- mtcars[, c("cyl", "hp", "am", "gear", "disp")]
d <- as(m, "DataFrame")
d$gr <- GenomicRanges::GRanges("chrY", IRanges::IRanges(1:32, width=10))
d$gr2 <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))
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
#>                           gr        gr2                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Mazda RX4          chrY:1-10  chrX:1-10  0.43,-1.85,0.22,...
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11 0.52,-1.62,-1.08,...
#> Datsun 710         chrY:3-12  chrX:3-12  0.46,-0.41,-2.1,...
#> Hornet 4 Drive     chrY:4-13  chrX:4-13     0.85,-1.15,-1.35
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.4,0.26,-1.55
#> ...                      ...        ...                  ...
#> Lotus Europa      chrY:28-37 chrX:28-37   1.38,1.92,1.67,...
#> Ford Pantera L    chrY:29-38 chrX:29-38  -1.43,0.1,-0.02,...
#> Ferrari Dino      chrY:30-39 chrX:30-39  -0.76,1.74,1.32,...
#> Maserati Bora     chrY:31-40 chrX:31-40  -0.39,0.87,2.21,...
#> Volvo 142E        chrY:32-41 chrX:32-41   1.51,0.42,1.48,...
```

This will appear in RStudio’s environment pane as a `Formal class
DataFrame (dplyr-compatible)` when using `DFplyr`. No interference with
the actual object is required, but this helps identify that
`dplyr`-compatibility is available.

`DataFrame`s can then be used in `dplyr` calls the same as `data.frame`
or `tibble` objects

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
#>                           gr        gr2                   nl    newvar
#>                    <GRanges>  <GRanges>        <NumericList> <numeric>
#> Mazda RX4          chrY:1-10  chrX:1-10  0.43,-1.85,0.22,...       116
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11 0.52,-1.62,-1.08,...       116
#> Datsun 710         chrY:3-12  chrX:3-12  0.46,-0.41,-2.1,...        97
#> Hornet 4 Drive     chrY:4-13  chrX:4-13     0.85,-1.15,-1.35       116
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.4,0.26,-1.55       183
#> ...                      ...        ...                  ...       ...
#> Lotus Europa      chrY:28-37 chrX:28-37   1.38,1.92,1.67,...       117
#> Ford Pantera L    chrY:29-38 chrX:29-38  -1.43,0.1,-0.02,...       272
#> Ferrari Dino      chrY:30-39 chrX:30-39  -0.76,1.74,1.32,...       181
#> Maserati Bora     chrY:31-40 chrX:31-40  -0.39,0.87,2.21,...       343
#> Volvo 142E        chrY:32-41 chrX:32-41   1.51,0.42,1.48,...       113
```

This includes piped calls with `%>%`

``` r
mutate(d, newvar = cyl + hp) %>%
  pull(newvar)
#>  [1] 116 116  97 116 183 111 253  66  99 129 129 188 188 188 213 223 238
#> [18]  70  56  69 101 158 158 253 183  70  95 117 272 181 343 113
```

and `tidyselect` helpers.

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
#>                           gr        gr2                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Mazda RX4          chrY:1-10  chrX:1-10  0.43,-1.85,0.22,...
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11 0.52,-1.62,-1.08,...
#> Datsun 710         chrY:3-12  chrX:3-12  0.46,-0.41,-2.1,...
#> Hornet 4 Drive     chrY:4-13  chrX:4-13     0.85,-1.15,-1.35
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.4,0.26,-1.55
#> ...                      ...        ...                  ...
#> Lotus Europa      chrY:28-37 chrX:28-37   1.38,1.92,1.67,...
#> Ford Pantera L    chrY:29-38 chrX:29-38  -1.43,0.1,-0.02,...
#> Ferrari Dino      chrY:30-39 chrX:30-39  -0.76,1.74,1.32,...
#> Maserati Bora     chrY:31-40 chrX:31-40  -0.39,0.87,2.21,...
#> Volvo 142E        chrY:32-41 chrX:32-41   1.51,0.42,1.48,...
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
#>                           gr        gr2                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Mazda RX4          chrY:1-10  chrX:1-10  0.43,-1.85,0.22,...
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11 0.52,-1.62,-1.08,...
#> Datsun 710         chrY:3-12  chrX:3-12  0.46,-0.41,-2.1,...
#> Hornet 4 Drive     chrY:4-13  chrX:4-13     0.85,-1.15,-1.35
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.4,0.26,-1.55
#> ...                      ...        ...                  ...
#> Lotus Europa      chrY:28-37 chrX:28-37   1.38,1.92,1.67,...
#> Ford Pantera L    chrY:29-38 chrX:29-38  -1.43,0.1,-0.02,...
#> Ferrari Dino      chrY:30-39 chrX:30-39  -0.76,1.74,1.32,...
#> Maserati Bora     chrY:31-40 chrX:31-40  -0.39,0.87,2.21,...
#> Volvo 142E        chrY:32-41 chrX:32-41   1.51,0.42,1.48,...

group_by(d, cyl) %>% 
  top_n(1, disp)
#> dplyr-compatible DataFrame with 3 rows and 8 columns
#>                          cyl        hp        am      gear      disp
#>                    <numeric> <numeric> <numeric> <numeric> <numeric>
#> Merc 240D                  4        62         0         4     146.7
#> Hornet 4 Drive             6       110         0         3       258
#> Cadillac Fleetwood         8       205         0         3       472
#>                            gr        gr2                 nl
#>                     <GRanges>  <GRanges>      <NumericList>
#> Merc 240D           chrY:8-17  chrX:8-17 -0.3,1.06,1.98,...
#> Hornet 4 Drive      chrY:4-13  chrX:4-13   0.85,-1.15,-1.35
#> Cadillac Fleetwood chrY:15-24 chrX:15-24   0.46,-1.23,-1.25
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
#>                           gr        gr2                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Maserati Bora     chrY:31-40 chrX:31-40  -0.39,0.87,2.21,...
#> Ford Pantera L    chrY:29-38 chrX:29-38  -1.43,0.1,-0.02,...
#> Duster 360         chrY:7-16  chrX:7-16       0.99,0.07,0.59
#> Camaro Z28        chrY:24-33 chrX:24-33       1.23,0.08,0.73
#> Chrysler Imperial chrY:17-26 chrX:17-26      -0.91,0.55,0.48
#> ...                      ...        ...                  ...
#> Fiat 128          chrY:18-27 chrX:18-27  0.12,-1.61,0.74,...
#> Fiat X1-9         chrY:26-35 chrX:26-35  -1.03,0.6,-0.97,...
#> Toyota Corolla    chrY:20-29 chrX:20-29   0.4,-0.53,0.26,...
#> Merc 240D          chrY:8-17  chrX:8-17   -0.3,1.06,1.98,...
#> Honda Civic       chrY:19-28 chrX:19-28 -0.48,0.64,-1.94,...

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
#>                           gr        gr2                 nl
#>                    <GRanges>  <GRanges>      <NumericList>
#> Hornet 4 Drive     chrY:4-13  chrX:4-13   0.85,-1.15,-1.35
#> Hornet Sportabout  chrY:5-14  chrX:5-14     0.4,0.26,-1.55
#> Valiant            chrY:6-15  chrX:6-15   -1.03,-1.07,2.49
#> Duster 360         chrY:7-16  chrX:7-16     0.99,0.07,0.59
#> Merc 240D          chrY:8-17  chrX:8-17 -0.3,1.06,1.98,...
#> ...                      ...        ...                ...
#> Toyota Corona     chrY:21-30 chrX:21-30   -0.68,-0.7,-1.83
#> Dodge Challenger  chrY:22-31 chrX:22-31   -0.29,-1.55,0.21
#> AMC Javelin       chrY:23-32 chrX:23-32   -0.23,-0.23,0.12
#> Camaro Z28        chrY:24-33 chrX:24-33     1.23,0.08,0.73
#> Pontiac Firebird  chrY:25-34 chrX:25-34    -0.78,0.46,0.16

slice(d, 3:6)
#> dplyr-compatible DataFrame with 4 rows and 8 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> Valiant                   6       105         0         3       225
#>                          gr       gr2                  nl
#>                   <GRanges> <GRanges>       <NumericList>
#> Datsun 710        chrY:3-12 chrX:3-12 0.46,-0.41,-2.1,...
#> Hornet 4 Drive    chrY:4-13 chrX:4-13    0.85,-1.15,-1.35
#> Hornet Sportabout chrY:5-14 chrX:5-14      0.4,0.26,-1.55
#> Valiant           chrY:6-15 chrX:6-15    -1.03,-1.07,2.49

# dd <- rbind(data.frame(m[1, ], row.names = "MyCar"), d)
# dd
```

Row names are not preserved when there may be duplicates

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
#>                           gr        gr2                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Mazda RX4          chrY:1-10  chrX:1-10  0.43,-1.85,0.22,...
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11 0.52,-1.62,-1.08,...
#> Datsun 710         chrY:3-12  chrX:3-12  0.46,-0.41,-2.1,...
#> Hornet 4 Drive     chrY:4-13  chrX:4-13     0.85,-1.15,-1.35
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.4,0.26,-1.55
#> ...                      ...        ...                  ...
#> Lotus Europa      chrY:28-37 chrX:28-37   1.38,1.92,1.67,...
#> Ford Pantera L    chrY:29-38 chrX:29-38  -1.43,0.1,-0.02,...
#> Ferrari Dino      chrY:30-39 chrX:30-39  -0.76,1.74,1.32,...
#> Maserati Bora     chrY:31-40 chrX:31-40  -0.39,0.87,2.21,...
#> Volvo 142E        chrY:32-41 chrX:32-41   1.51,0.42,1.48,...

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

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

## S4 Columns

Support for S4 columns is tentative. Printing works, but not yet working
with them.

``` r
d2 <- d
library(GenomicRanges)
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following objects are masked from 'package:DFplyr':
#> 
#>     desc, slice
#> The following objects are masked from 'package:dplyr':
#> 
#>     collapse, desc, slice
#> Loading required package: GenomeInfoDb
d2$gr <- GRanges("chrY", IRanges(1:32, width=10))
d2
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
#>                           gr        gr2                   nl
#>                    <GRanges>  <GRanges>        <NumericList>
#> Mazda RX4          chrY:1-10  chrX:1-10  0.43,-1.85,0.22,...
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11 0.52,-1.62,-1.08,...
#> Datsun 710         chrY:3-12  chrX:3-12  0.46,-0.41,-2.1,...
#> Hornet 4 Drive     chrY:4-13  chrX:4-13     0.85,-1.15,-1.35
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.4,0.26,-1.55
#> ...                      ...        ...                  ...
#> Lotus Europa      chrY:28-37 chrX:28-37   1.38,1.92,1.67,...
#> Ford Pantera L    chrY:29-38 chrX:29-38  -1.43,0.1,-0.02,...
#> Ferrari Dino      chrY:30-39 chrX:30-39  -0.76,1.74,1.32,...
#> Maserati Bora     chrY:31-40 chrX:31-40  -0.39,0.87,2.21,...
#> Volvo 142E        chrY:32-41 chrX:32-41   1.51,0.42,1.48,...
```

## Implementation

Most of the `dplyr` verbs for `DataFrame`s are implmented by first
converting to `tibble`, performing the verb operation, then converting
back to `DataFrame`. Care has been taken to retain groups and row names
through these operations, but this may introduce some complications. If
you spot any, please [file an
issue](https://github.com/jonocarroll/DFplyr/issues/new).
