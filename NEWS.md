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


