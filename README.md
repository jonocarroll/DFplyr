
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
m <- mtcars[, c("cyl", "hp", "am", "gear", "disp")]
d <- as(m, "DataFrame")
d
#> dplyr-compatible DataFrame with 32 rows and 5 columns
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
```

This will appear in RStudio’s environment pane as a `Formal class
DataFrame (dplyr-compatible)` when using `DFplyr`. No interference with
the actual object is required, but this helps identify that
`dplyr`-compatibility is available.

`DataFrame`s can then be used in `dplyr` calls the same as `data.frame`
or `tibble` objects

``` r
mutate(d, newvar = cyl + hp)
#> dplyr-compatible DataFrame with 32 rows and 6 columns
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
#>                      newvar
#>                   <numeric>
#> Mazda RX4               116
#> Mazda RX4 Wag           116
#> Datsun 710               97
#> Hornet 4 Drive          116
#> Hornet Sportabout       183
#> ...                     ...
#> Lotus Europa            117
#> Ford Pantera L          272
#> Ferrari Dino            181
#> Maserati Bora           343
#> Volvo 142E              113
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
#> dplyr-compatible DataFrame with 32 rows and 5 columns
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
```

Importantly, grouped operations are supported. `DataFrame` does not
natively support groups (the same way that `data.frame` does not) so
these are implemented specifically for `DFplyr`

``` r
group_by(d, cyl, am)
#> dplyr-compatible DataFrame with 32 rows and 5 columns
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

group_by(d, cyl) %>% 
  top_n(1, disp)
#> dplyr-compatible DataFrame with 3 rows and 5 columns
#>                          cyl        hp        am      gear      disp
#>                    <numeric> <numeric> <numeric> <numeric> <numeric>
#> Hornet 4 Drive             6       110         0         3       258
#> Merc 240D                  4        62         0         4     146.7
#> Cadillac Fleetwood         8       205         0         3       472
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
#> dplyr-compatible DataFrame with 32 rows and 5 columns
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

filter(d, am == 0) 
#> dplyr-compatible DataFrame with 19 rows and 5 columns
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

slice(d, 3:6)
#> dplyr-compatible DataFrame with 4 rows and 5 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> Valiant                   6       105         0         3       225

dd <- rbind(DataFrame(d[1, ], row.names = "MyCar"), d)
dd
#> dplyr-compatible DataFrame with 33 rows and 5 columns
#>                      cyl        hp        am      gear      disp
#>                <numeric> <numeric> <numeric> <numeric> <numeric>
#> MyCar                  6       110         1         4       160
#> Mazda RX4              6       110         1         4       160
#> Mazda RX4 Wag          6       110         1         4       160
#> Datsun 710             4        93         1         4       108
#> Hornet 4 Drive         6       110         0         3       258
#> ...                  ...       ...       ...       ...       ...
#> Lotus Europa           4       113         1         5      95.1
#> Ford Pantera L         8       264         1         5       351
#> Ferrari Dino           6       175         1         5       145
#> Maserati Bora          8       335         1         5       301
#> Volvo 142E             4       109         1         4       121
```

Row names are not preserved when there may be duplicates

``` r
distinct(dd)
#> dplyr-compatible DataFrame with 28 rows and 5 columns
#>           cyl        hp        am      gear      disp
#>     <numeric> <numeric> <numeric> <numeric> <numeric>
#> 1           6       110         1         4       160
#> 2           4        93         1         4       108
#> 3           6       110         0         3       258
#> 4           8       175         0         3       360
#> 5           6       105         0         3       225
#> ...       ...       ...       ...       ...       ...
#> 24          4       113         1         5      95.1
#> 25          8       264         1         5       351
#> 26          6       175         1         5       145
#> 27          8       335         1         5       301
#> 28          4       109         1         4       121

group_by(d, cyl, am) %>%
  tally(gear)
#> # A tibble: 6 x 3
#> # Groups:   cyl [3]
#>     cyl    am     n
#>   <dbl> <dbl> <dbl>
#> 1     4     0    11
#> 2     4     1    34
#> 3     6     0    14
#> 4     6     1    13
#> 5     8     0    36
#> 6     8     1    10

count(d, gear, am, cyl)
#> # A tibble: 10 x 4
#>     gear    am   cyl     n
#>    <dbl> <dbl> <dbl> <int>
#>  1     3     0     4     1
#>  2     3     0     6     2
#>  3     3     0     8    12
#>  4     4     0     4     2
#>  5     4     0     6     2
#>  6     4     1     4     6
#>  7     4     1     6     2
#>  8     5     1     4     2
#>  9     5     1     6     1
#> 10     5     1     8     2
```

## S4 Columns

Support for S4 columns is preliminary. Printing works, and a `DataFrame`
containing them may be used in verb operations, but the verbs cannot
(currently) act directly on them (via `...`)

``` r
d2 <- d
d2$gr <- GenomicRanges::GRanges("chrY", IRanges::IRanges(1:32, width=10))
d2$gr2 <- GenomicRanges::GRanges("chrX", IRanges::IRanges(1:32, width = 10))
d2$nl <- IRanges::NumericList(lapply(d2$gear, function(n) round(rnorm(n), 2)))
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
#>                           gr        gr2                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrY:1-10  chrX:1-10   -0.09,-1.1,1.53,...
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11  -0.17,-0.74,2.48,...
#> Datsun 710         chrY:3-12  chrX:3-12    1.15,-0.27,0.4,...
#> Hornet 4 Drive     chrY:4-13  chrX:4-13       -0.34,0.05,-0.6
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.83,-0.24,1.59
#> ...                      ...        ...                   ...
#> Lotus Europa      chrY:28-37 chrX:28-37 -2.84,-0.58,-1.32,...
#> Ford Pantera L    chrY:29-38 chrX:29-38    -0.14,0.3,0.19,...
#> Ferrari Dino      chrY:30-39 chrX:30-39    -0.7,1.93,1.85,...
#> Maserati Bora     chrY:31-40 chrX:31-40    -0.9,0.36,-1.3,...
#> Volvo 142E        chrY:32-41 chrX:32-41  -0.26,-0.32,0.73,...

filter(d2, cyl == 6)
#> dplyr-compatible DataFrame with 7 rows and 8 columns
#>                      cyl        hp        am      gear      disp
#>                <numeric> <numeric> <numeric> <numeric> <numeric>
#> Mazda RX4              6       110         1         4       160
#> Mazda RX4 Wag          6       110         1         4       160
#> Hornet 4 Drive         6       110         0         3       258
#> Valiant                6       105         0         3       225
#> Merc 280               6       123         0         4     167.6
#> Merc 280C              6       123         0         4     167.6
#> Ferrari Dino           6       175         1         5       145
#>                        gr        gr2                   nl
#>                 <GRanges>  <GRanges>        <NumericList>
#> Mazda RX4       chrY:1-10  chrX:1-10  -0.09,-1.1,1.53,...
#> Mazda RX4 Wag   chrY:2-11  chrX:2-11 -0.17,-0.74,2.48,...
#> Hornet 4 Drive  chrY:4-13  chrX:4-13      -0.34,0.05,-0.6
#> Valiant         chrY:6-15  chrX:6-15      1.39,1.11,-0.33
#> Merc 280       chrY:10-19 chrX:10-19 -2.4,-0.55,-1.17,...
#> Merc 280C      chrY:11-20 chrX:11-20  -1.1,-0.94,1.45,...
#> Ferrari Dino   chrY:30-39 chrX:30-39   -0.7,1.93,1.85,...

select(d2, cyl:am)
#> dplyr-compatible DataFrame with 32 rows and 3 columns
#>                         cyl        hp        am
#>                   <numeric> <numeric> <numeric>
#> Mazda RX4                 6       110         1
#> Mazda RX4 Wag             6       110         1
#> Datsun 710                4        93         1
#> Hornet 4 Drive            6       110         0
#> Hornet Sportabout         8       175         0
#> ...                     ...       ...       ...
#> Lotus Europa              4       113         1
#> Ford Pantera L            8       264         1
#> Ferrari Dino              6       175         1
#> Maserati Bora             8       335         1
#> Volvo 142E                4       109         1

mutate(d2, newcol = rnorm(32))
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
#>                           gr        gr2                    nl
#>                    <GRanges>  <GRanges>         <NumericList>
#> Mazda RX4          chrY:1-10  chrX:1-10   -0.09,-1.1,1.53,...
#> Mazda RX4 Wag      chrY:2-11  chrX:2-11  -0.17,-0.74,2.48,...
#> Datsun 710         chrY:3-12  chrX:3-12    1.15,-0.27,0.4,...
#> Hornet 4 Drive     chrY:4-13  chrX:4-13       -0.34,0.05,-0.6
#> Hornet Sportabout  chrY:5-14  chrX:5-14       0.83,-0.24,1.59
#> ...                      ...        ...                   ...
#> Lotus Europa      chrY:28-37 chrX:28-37 -2.84,-0.58,-1.32,...
#> Ford Pantera L    chrY:29-38 chrX:29-38    -0.14,0.3,0.19,...
#> Ferrari Dino      chrY:30-39 chrX:30-39    -0.7,1.93,1.85,...
#> Maserati Bora     chrY:31-40 chrX:31-40    -0.9,0.36,-1.3,...
#> Volvo 142E        chrY:32-41 chrX:32-41  -0.26,-0.32,0.73,...
#>                               newcol
#>                            <numeric>
#> Mazda RX4         -0.646170522988381
#> Mazda RX4 Wag     -0.363572590356091
#> Datsun 710        0.0760921143103913
#> Hornet 4 Drive      1.09176897631427
#> Hornet Sportabout  0.205572114336237
#> ...                              ...
#> Lotus Europa      0.0855617764314784
#> Ford Pantera L      0.32212307936588
#> Ferrari Dino       0.181316280960277
#> Maserati Bora      -1.12336658458553
#> Volvo 142E        -0.321001162771463
```

`DataFrame` objects are still useful as the `.data` argument

``` r
mtcars_desc <- tibble::tribble(
  ~col,    ~desc,
  "mpg",   "Miles/(US) gallon",
  "cyl",   "Number of cylinders",
  "disp",  "Displacement (cu.in.)",
  "hp",    "Gross horsepower",
  "drat",  "Rear axle ratio",
  "wt",    "Weight (1000 lbs)",
  "qsec",  "1/4 mile time",
  "vs",    "Engine (0 = V-shaped, 1 = straight)",
  "am",    "Transmission (0 = automatic, 1 = manual)",
  "gear",  "Number of forward gears",
  "carb",  "Number of carburetors"
)
mcols(d) <- mtcars_desc[match(names(d), mtcars_desc$col), "desc"]
mcols(d)
#> dplyr-compatible DataFrame with 5 rows and 1 column
#>                                          desc
#>                                   <character>
#> cyl                       Number of cylinders
#> hp                           Gross horsepower
#> am   Transmission (0 = automatic, 1 = manual)
#> gear                  Number of forward gears
#> disp                    Displacement (cu.in.)

mutate(mcols(d), varname = tools::toTitleCase(desc))
#> dplyr-compatible DataFrame with 5 rows and 2 columns
#>                                          desc
#>                                   <character>
#> cyl                       Number of cylinders
#> hp                           Gross horsepower
#> am   Transmission (0 = automatic, 1 = manual)
#> gear                  Number of forward gears
#> disp                    Displacement (cu.in.)
#>                                       varname
#>                                   <character>
#> cyl                       Number of Cylinders
#> hp                           Gross Horsepower
#> am   Transmission (0 = Automatic, 1 = Manual)
#> gear                  Number of Forward Gears
#> disp                    Displacement (Cu.in.)
```

Some exceptions apply - the `pull` verb *does* work with an S4 column

``` r
pull(d2, gr)
#> GRanges object with 32 ranges and 0 metadata columns:
#>        seqnames    ranges strand
#>           <Rle> <IRanges>  <Rle>
#>    [1]     chrY      1-10      *
#>    [2]     chrY      2-11      *
#>    [3]     chrY      3-12      *
#>    [4]     chrY      4-13      *
#>    [5]     chrY      5-14      *
#>    ...      ...       ...    ...
#>   [28]     chrY     28-37      *
#>   [29]     chrY     29-38      *
#>   [30]     chrY     30-39      *
#>   [31]     chrY     31-40      *
#>   [32]     chrY     32-41      *
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

pull(d2, nl)
#> NumericList of length 32
#> [[1]] -0.09 -1.1 1.53 0.98
#> [[2]] -0.17 -0.74 2.48 0.19
#> [[3]] 1.15 -0.27 0.4 -0.76
#> [[4]] -0.34 0.05 -0.6
#> [[5]] 0.83 -0.24 1.59
#> [[6]] 1.39 1.11 -0.33
#> [[7]] -1.3 0.32 1.05
#> [[8]] -0.95 1 0.59 -0.55
#> [[9]] 0.47 -1.35 -0.17 1.8
#> [[10]] -2.4 -0.55 -1.17 0.62
#> ...
#> <22 more elements>
```

## `ggplot2`

As an added bonus, `ggplot2` support is also enabled

``` r
library(ggplot2)
ggplot(d, aes(disp, cyl)) + geom_point()
```

<img src="man/figures/README-unnamed-chunk-13-1.png" width="100%" />

## Implementation

Most of the `dplyr` verbs for `DataFrame`s are implmented by first
converting to `tibble`, performing the verb operation, then converting
back to `DataFrame`. Care has been taken to retain groups and row names
and preserve S4 columns through these operations, but this is fragile
and may introduce some complications. If you spot any, please [file an
issue](https://github.com/jonocarroll/DFplyr/issues/new).
