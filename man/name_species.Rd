% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{name_species}
\alias{name_species}
\title{Add a column with species names}
\usage{
name_species(data, sp.names)
}
\arguments{
\item{data}{An object of class nelder-biculture or goelz.}

\item{sp.names}{A vector of species names. Must be of length 2 when data is of class nelder-biculture and length 3 when data is of class goelz.}
}
\value{
An object of the same class as data.
}
\description{
Adds a column with species names.
}
\author{
Kevin J Wolz, \email{kevin@savannainstitute.org}
}
