% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mathml.R
\name{fontstyles}
\alias{fontstyles}
\alias{plain}
\alias{italic}
\alias{bold}
\title{Identity functions for different font styles}
\usage{
plain(x)

italic(x)

bold(x)
}
\arguments{
\item{x}{the expression to render}
}
\value{
x
}
\description{
Identity functions for different font styles
}
\examples{
plain(1) + bold(2) + italic(3)

mathml(term=quote(plain(abc) + bold(def) + italic(ghi)))
}
