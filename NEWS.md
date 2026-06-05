# DFplyr 1.7.0

NEW FEATURES

* summarise, count, mutate can handle non-standard column names (#26, @smped)
* mutate can refer to previous columns within the same call to mutate (#26, @smped)
* GroupedDataFrame is valid for 1 group with 1 value, and above (#26, @smped)
* wt and sort enabled for count (#26, @smped)
* Explicit dependency on R >= 4.1 (due to use of `|>`)

# DFplyr 1.6.0

NEW FEATURES

* Fixed column-wise slicing
* Aligned to S4Vectors updates

SIGNIFICANT USER-VISIBLE CHANGES

* Grouped DataFrame now appears as a new class GroupedDataFrame with its own dispatched methods 

# DFplyr 1.4.0

NEW FEATURES

* `[` and `rbind` methods now respect groups
* Added `left_join()`, `right_join()`, `inner_join()`, and `full_join()`

SIGNIFICANT USER-VISIBLE CHANGES

* `rename2` deprecated; `rename` is now a `DataFrame` method on `S4Vectors` generic

# DFplyr 0.0.3 / 1.2.0

NEW FEATURES

* Preparation for Bioconductor release

SIGNIFICANT USER-VISIBLE CHANGES

BUG FIXES

* `summarise()` now prints the return value

# DFplyr 0.0.2

NEW FEATURES

* Added `tally()`
* Better importing of generics

SIGNIFICANT USER-VISIBLE CHANGES

* `rename()` is now `rename2()` to avoid conflicts with `S4Vectors::rename`

BUG FIXES

* `slice()` now respects group information

# DFplyr 0.0.1

NEW FEATURES

* Initial working version
* Some group information is not respected


