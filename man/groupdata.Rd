% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr.R
\name{group_data}
\alias{group_data}
\alias{set_group_data}
\alias{get_group_data}
\title{Set and Get Group Data on a DataFrame}
\usage{
set_group_data(x, g, .drop = group_by_drop_default(x))

get_group_data(x)
}
\arguments{
\item{x}{A \link[S4Vectors:DataFrame-class]{S4Vectors::DataFrame} on which to set group data.}

\item{g}{Group data (a \code{data.frame}).}

\item{.drop}{Drop groups formed by factor levels that don't appear in
the data?}
}
\value{
For \code{set_group_data}, the input \code{x} with group data set as metadata.
For \code{get_group_data}, the group data that is set on \code{x}.
}
\description{
The location of group data is an internal implemnetation
detail, so these get and set methods enable interfacing
with that data.
}
\keyword{internal}
