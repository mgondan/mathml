% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mathml.R
\name{prod_over}
\alias{prod_over}
\title{product over a range. On the R side, this function just returns the product
of the first argument, but allows for decorations.}
\usage{
prod_over(x, from, to)
}
\arguments{
\item{x}{the object to be multiplied}

\item{from}{decoration for prod_from^to x_i}

\item{to}{decoration for prod_from^to x_i}
}
\value{
The function returns prod(x)
}
\description{
product over a range. On the R side, this function just returns the product
of the first argument, but allows for decorations.
}
\examples{
mathjax(quote(prod_over(x[i], i=1L, N)))

}
\seealso{
\code{\link[=prod]{prod()}}, \code{\link[=sum_over]{sum_over()}}
}
