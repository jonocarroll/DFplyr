
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DFplyr

<!-- badges: start -->

<!-- badges: end -->

The goal of DFplyr is to enable `dplyr` support for
`S4Vectors::DataFrame`.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jonocarroll/DFplyr")
```

## Examples

``` r
suppressPackageStartupMessages({
  library(S4Vectors)
  library(dplyr)
  library(DFplyr)
})
#> Warning: replacing previous import 'S4Vectors::union' by 'dplyr::union'
#> when loading 'DFplyr'
#> Warning: replacing previous import 'S4Vectors::intersect' by
#> 'dplyr::intersect' when loading 'DFplyr'
#> Warning: replacing previous import 'S4Vectors::setdiff' by 'dplyr::setdiff'
#> when loading 'DFplyr'
#> Warning: replacing previous import 'S4Vectors::first' by 'dplyr::first'
#> when loading 'DFplyr'
#> Warning: replacing previous import 'S4Vectors::setequal' by
#> 'dplyr::setequal' when loading 'DFplyr'
#> Warning: replacing previous import 'S4Vectors::rename' by 'dplyr::rename'
#> when loading 'DFplyr'
#> Warning: multiple methods tables found for 'rename'

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
#>  1     3     0     4     0
#>  2     3     0     6     0
#>  3     3     0     8     0
#>  4     4     0     4     0
#>  5     4     0     6     0
#>  6     4     1     4     0
#>  7     4     1     6     0
#>  8     5     1     4     0
#>  9     5     1     6     0
#> 10     5     1     8     0

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
#> Mazda RX4                 8       335         1         5       301
#> Mazda RX4 Wag             8       264         1         5       351
#> Datsun 710                8       245         0         3       360
#> Hornet 4 Drive            8       245         0         3       350
#> Hornet Sportabout         8       230         0         3       440
#> ...                     ...       ...       ...       ...       ...
#> Lotus Europa              4        66         1         4      78.7
#> Ford Pantera L            4        66         1         4        79
#> Ferrari Dino              4        65         1         4      71.1
#> Maserati Bora             4        62         0         4     146.7
#> Volvo 142E                4        52         1         4      75.7

# ## row names are not preserved as there may be duplicates 
# rbind(data.frame(m[1, ], row.names = "MyCar"), d) %>%
#   distinct()

filter(d, am == 0) 
#> dplyr-compatible DataFrame with 19 rows and 5 columns
#>           cyl        hp        am      gear      disp
#>     <numeric> <numeric> <numeric> <numeric> <numeric>
#> 1           6       110         0         3       258
#> 2           8       175         0         3       360
#> 3           6       105         0         3       225
#> 4           8       245         0         3       360
#> 5           4        62         0         4     146.7
#> ...       ...       ...       ...       ...       ...
#> 15          4        97         0         3     120.1
#> 16          8       150         0         3       318
#> 17          8       150         0         3       304
#> 18          8       245         0         3       350
#> 19          8       175         0         3       400

slice(d, 3:6)
#> dplyr-compatible DataFrame with 4 rows and 5 columns
#>                         cyl        hp        am      gear      disp
#>                   <numeric> <numeric> <numeric> <numeric> <numeric>
#> Datsun 710                4        93         1         4       108
#> Hornet 4 Drive            6       110         0         3       258
#> Hornet Sportabout         8       175         0         3       360
#> Valiant                   6       105         0         3       225
```

`ggplot2` support is also enabled

``` r
library(ggplot2)
#> 
#> Attaching package: 'ggplot2'
#> The following object is masked from 'package:DFplyr':
#> 
#>     fortify
ggplot(d, aes(disp, cyl)) + geom_point()
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="100%" />
