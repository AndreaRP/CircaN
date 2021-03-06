% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{jtk_wrapper}
\alias{jtk_wrapper}
\title{jtk_wrapper}
\usage{
jtk_wrapper(data, s2c, min_per = 20, max_per = 28)
}
\arguments{
\item{data}{The data matrix with the omics data.}

\item{s2c}{Data frame with at least the columns sample (which must coincide with the samples in data) and time (recording
time corresponding to each sample.}
}
\description{
Wraps JTK method. The period is set between 20 and 28h.
}
\examples{
jtk_results <- jtk_wrapper(my_data, s2c)
}
\keyword{jtk}
