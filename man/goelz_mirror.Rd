% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/goelz.R
\name{goelz_mirror}
\alias{goelz_mirror}
\title{Mirror a Goelz Triangle experimental design}
\usage{
goelz_mirror(data, joining.borders = max(data$border.num, na.rm = TRUE))
}
\arguments{
\item{data}{An object of class goelz.}

\item{joining.borders}{The number of border rows to maintain between the two triangles.}
}
\value{
A data frame (tibble) of class sysd, goelz, and goelz-mirror.
}
\description{
Mirrors a Goelz Triangle experimental design and combines it with the original design to create a
parallelogram design.
}
\details{
When establishing two or more replicate Goelz Triangle plots, the plots can be stacked next to each other,
with species matched, to further reduce edge effects. This function takes a single Goelz Triangle, mirrors it to
create a second triangle, and then combines the two together into a parallelogram design. When merging the two
designs, the border rows between the two triangles, if there are any, "double up". This may or may not be desired.
Border rows between the two triangles is not really necessary for reducing edge effects, since this edge is now on
the interior of the parallelogram. However, at least some borders between triangles may be desired to "separate"
the two replicates and "maintain" independence. Any approach can be created by specifiying the exact number of
border rows to maintain betwen the two triangles via the \code{joining.borders} argument. The default is to maintain
only the number of border rows that each triangle originally had (i.e. half the number of borders initially
between the two triangles when merging them.)
}
\examples{
dat <- goelz()
dat.border <- goelz_add_border(data = dat, n = 3)
dat.border.mirror <- goelz_mirror(data = dat.border)
}
\references{
\itemize{
 \item Goelz JCG (2001) Systematic experimental designs for mixed species plantings. Native Plants Journal 2:90–96.
 \url{http://npj.uwpress.org/content/2/2/90.short}
}
}
\seealso{
Other definition functions: 
\code{\link{goelz_add_border}()},
\code{\link{goelz_corners}()},
\code{\link{goelz_guides}()},
\code{\link{goelz_optim}()},
\code{\link{goelz_starts}()},
\code{\link{goelz}()},
\code{\link{nelder_biculture_competition}()},
\code{\link{nelder_biculture_optim}()},
\code{\link{nelder_biculture}()},
\code{\link{nelder_decision}()},
\code{\link{nelder_interspoke_distance}()},
\code{\link{nelder}()},
\code{\link{select_optimal_goelz}()},
\code{\link{select_optimal_nelder_biculture}()}
}
\author{
Kevin J Wolz, \email{kevin@savannainstitute.org}
}
\concept{definition functions}
