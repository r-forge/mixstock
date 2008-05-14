\name{mysum}
\alias{mysum}
\alias{mysumvec}
\title{summaries of MCMC output}
\description{
  calculates summary statistics of one or more columns of
  MCMC output
}
\usage{
mysum(x, names = NULL)
mysumvec(x,names=NULL)
}
\arguments{
  \item{x}{a numeric matrix representing MCMC output}
  \item{names}{vector of column names}
}
\details{
  a slight variant on the standard numeric summary: calculates
  mean, median, standard deviation, and lower and upper 90 and 95 percent
  quantiles. \code{mysum} gives a table (original data as rows,
  various stats as columns), while \code{mysumvec} collapses the
  answer to a vector
}
\value{
  a numeric vector of results, appropriately named
}
\examples{
  x = matrix(runif(1000),ncol=2,dimnames=list(NULL,c("A","B")))
  mysum(x)
  mysumvec(x)
}
\keyword{misc}
