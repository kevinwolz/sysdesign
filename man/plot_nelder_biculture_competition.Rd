% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_nelder_biculture_competition}
\alias{plot_nelder_biculture_competition}
\title{Plot a diagnositc display of compeitition in a Nelder Fan biculture experimental design}
\usage{
plot_nelder_biculture_competition(data, plot = TRUE)
}
\arguments{
\item{data}{An object of class 'nelder-biculture' created via \code{\link{goelz}} with \code{split} = FALSE.}

\item{plot}{If \code{TRUE}, the default, a ggplot object is returned.
If \code{FALSE}, the data that would create the plot is returned.}
}
\value{
If \code{plot = TRUE}, returns a ggplot object, otherwise the data that would create the plot is returned.
}
\description{
Plots a diagnositc display of compeitition in a Nelder Fan biculture experimental design.
The optimal competition scenario is to have counts consistent across arcs. See \code{\link{nelder_optim}} for more information.
}
\examples{
dat <- nelder()
dat.bi <- nelder_biculture(data = dat)
plot_nelder_biculture_competition(data = dat.bi)
}
\seealso{
Other plot functions: 
\code{\link{plot_fitness_trajectory}()},
\code{\link{plot_goelz}()},
\code{\link{plot_nelder}()}
}
\author{
Kevin J Wolz, \email{kevin@savannainstitute.org}
}
\concept{plot functions}
