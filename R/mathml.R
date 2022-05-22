.onAttach <- function(libname, pkgname)
{
  if(!requireNamespace("rolog", quietly=TRUE))
    stop("Could not attach library rolog.")

  rolog::consult(system.file("pl/mathml.pl", package=pkgname))
}

#' MathML output
#'
#' @param term
#' an R call or symbol/number. This function translates _term_ into a
#' MathML string.
#'
#' @return
#' A string with the MathML representation or _term_.
#'
#' @md
#'
#' @seealso [mathjax()]
#'
#' @examples
#' mathml(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L))
#'
mathml <- function(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L))
{
  t = rolog::once(call("r2mathml", term, expression(X)))
  cat(paste(t$X, collapse=""))
}

#' Mathjax output
#'
#' @param term
#' an R call or symbol/number. This function translates _term_ into a
#' LaTeX/MathJax string.
#'
#' @return
#' A string with the MathJax representation or _term_.
#'
#' @md
#'
#' @seealso [mathml()]
#'
#' @examples
#' mathjax(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L))
#'
mathjax <- function(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L))
{
  t = rolog::once(call("r2mathjax", term, expression(X)))
  cat(paste(t$X, collapse=""))
}

#' Calligraphic font
#'
#' @param x
#' an R symbol. This function is used to render the content in calligraphic font
#' in MathJax. In MathML, script font is used.
#'
#' @return
#' The function cal is a wrapper for the identity function.
#'
#' @md
#'
#' @seealso [identity()]
#'
#' @examples
#' mathjax(quote(K %in% cal(K)))
#'
cal <- identity

#' Sum over a range, e.g., sum x_i for i=1 to N
#'
#' @param index
#' an R symbol, e.g., i.
#'
#' @param from
#' an R symbol, e.g., 1.
#'
#' @param to
#' an R symbol, e.g., N.
#'
#' @param fun
#' an R call. This is the return value of the function.
#'
#' @return
#' The function over is a wrapper for the identity function, returning _fun_
#'
#' @md
#'
#' @seealso [identity()]
#'
#' @examples
#' mathjax(quote(over(index=i, from=1, to=N, fun=sum(1:10))))
#'
over <- function(index=quote(i), from=1, to=quote(N), fun=sum(1:10))
{
  return(fun)
}

#' Canonicalize an R call: Reorder and named the function arguments
#'
#' @param term
#' an R call.
#'
#' @return
#' The R function, with arguments rearranged
#'
#' @md
#'
#' @examples
#' canonical(term=quote(`%in%`(table=Table, x=X)))
#'
canonical <- function(term=quote(`%in%`(table=Table, x=X)))
{
  if(is.call(term))
  {
    f <- match.fun(term[[1]])
    if(!is.primitive(f))
      term <- match.call(f, term)
    term[-1] <- lapply(term[-1], canonical)
  }
  return(term)
}
