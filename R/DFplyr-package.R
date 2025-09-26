#' Treat a `S4Vectors::DataFrame` as a `dplyr` data source
#'
#' Add \pkg{dplyr} compatibility to `S4Vectors::DataFrame`
#' for use with a selection of \pkg{dplyr} verbs.
#'
#' @param x A `S4Vectors::DataFrame` object
#'
#' @importFrom dplyr mutate filter select arrange tbl_vars
#' @importFrom dplyr count group_by_drop_default summarise summarize
#' @importFrom dplyr group_data group_vars group_by ungroup
#' @importFrom dplyr distinct pull slice tally
#' @importFrom dplyr left_join right_join inner_join full_join
#' @import S4Vectors
#'
#' @md
#' @examples
#' library(S4Vectors)
#' library(dplyr)
#'
#' d <- as(mtcars, "DataFrame")
#'
#' mutate(d, newvar = cyl + hp)
#'
#' mutate_at(d, vars(starts_with("c")), ~ .^2)
#'
#' group_by(d, cyl, am) %>%
#'     tally(gear)
#'
#' count(d, gear, am, cyl)
#'
#' select(d, am, cyl)
#'
#' select(d, am, cyl) %>%
#'     rename2(foo = am)
#'
#' arrange(d, desc(hp))
#'
#' rbind(DataFrame(mtcars[1, ], row.names = "MyCar"), d) %>%
#'     distinct()
#'
#' filter(d, am == 0)
#'
#' slice(d, 3:6)
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
