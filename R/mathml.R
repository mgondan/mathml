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
#' @param flags (default NULL)
#' list of flags that control the translation
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
mathml <- function(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L), flags=NULL)
{
  flags <- c(attributes(term), flags, list(cat=FALSE))
  t <- rolog::once(call("r2mathml", flags, term, expression(X)),
    options=list(preproc=mathml_preproc))
  r <- paste(t$X, collapse="")
  if(flags$cat)
    return(cat(r))

  return(r)
}

# Prolog representation of not equal etc. (left: R, right: Prolog)
mathml_operators = c(
  "!=" = "\\=",
  "<=" = "=<",
  "%.%" = "cdot",
  "%/%" = "div",
  "%+-%" = "pm",
  "%*%" = "times",
  "%~~%" = "approx",
  "%==%" = "equiv",
  "%=~%" = "cong",
  "%prop%" = "propto",
  "%<->%" = "leftrightarrow",
  "%->%" = "rightarrow",
  "%<-%" = "leftarrow",
  "%up%" = "uparrow",
  "%down%" = "downarrow",
  "%<=>%" = "iff",
  "%=>%" = "rArr",
  "%<=%" = "lArr",
  "%dblup%" = "uArr",
  "%dbldown%" = "dArr"
)

mathml_preproc <- function(query=quote(2 != 2))
{
  if(is.call(query))
  {
    args <- as.list(query)

    index <- which(args[[1]] == names(mathml_operators))
    if(length(index) == 1)
      args[[1]] <- mathml_operators[index]

    args[-1] <- lapply(args[-1], FUN=mathml_preproc)
    return(as.call(args))
  }

  if(is.list(query))
    return(lapply(query, FUN=mathml_preproc))

  return(query)
}

#' Mathjax output
#'
#' @param term
#' an R call or symbol/number. This function translates _term_ into a
#' LaTeX/MathJax string.
#'
#' @param flags (default NULL)
#' list of flags that control the translation
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
mathjax <- function(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L), flags=NULL)
{
  flags <- c(attributes(term), flags, list(cat=FALSE))
  t <- rolog::once(call("r2mathjax", flags, term, expression(X)),
    options=list(preproc=mathml_preproc))
  r <- paste(t$X, collapse="")
  if(flags$cat)
    return(cat(r))

  return(r)
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

#' Subscript. On the R side, this function is a wrapper of identity, but allows
#' for decorations.
#'
#' @param sub
#' an R symbol or call, e.g., i
#'
#' @param fun
#' an R call or symbol, e.g. sum(x). This is the return value of the function.
#'
#' @return
#' The function over is a wrapper for the identity function, returning _fun_
#'
#' @md
#'
#' @seealso [identity()]
#'
#' @examples
#' mathjax(quote(subscript(sub=i, fun=x)))
#'
subscript <- function(fun=quote(x), sub=quote(i))
{
  return(fun)
}

#' Superscript. This is a wrapper for the identity function, but decorates the
#' result with a superscript.
#'
#' @param sup
#' an R symbol, e.g., "*"
#'
#' @param fun
#' an R call or symbol, e.g. x. This is the return value of the function.
#'
#' @return
#' The function over is a wrapper for the identity function, returning _fun_
#'
#' @md
#'
#' @seealso [identity()]
#'
#' @examples
#' mathjax(quote(superscript(fun=A, sup="*")))
#'
superscript <- function(fun=quote(A), sup="*")
{
  return(fun)
}

#' Subsupscript. This is a wrapper for the identity function, but decorates the
#' result with a sub- and a superscript.
#'
#' @md
#'
#' @param sub
#' an R symbol, e.g., `i=1`
#'
#' @param sup
#' an R symbol, e.g., `N`
#'
#' @param fun
#' an R call or symbol, e.g. `sum(x[i])`. This is the return value.
#'
#' @return
#' The function over is a wrapper for the identity function, returning _fun_
#'
#' @seealso [identity()]
#'
#' @examples
#' N <- 10
#' i <- 1:N
#' x <- rnorm(N)
#' mathjax(call("subsupscript", fun=sum(x[i]), sub=quote(`=`(i, 1L)), sup=quote(N)))
#'
subsupscript <- function(fun=quote(sum(x[i])), sub=quote(`=`(i, 1)), sup=quote(N))
{
  return(fun)
}

#' Canonicalize an R call: Reorder the function arguments
#'
#' @param term
#' an R call.
#'
#' @param drop
#' whether to drop the argument names or not
#'
#' @return
#' The R function, with arguments rearranged
#'
#' @md
#'
#' @examples
#' canonical(term=quote(`%in%`(table=Table, x=X)))
#'
canonical <- function(term=quote(`%in%`(table=Table, x=X)), drop=TRUE)
{
  attr <- attributes(term)

  if(is.call(term))
  {
    f <- match.fun(term[[1]])
    if(!is.primitive(f))
      term <- match.call(f, term)
    term[-1] <- lapply(term[-1], canonical, drop=drop)
  }

  if(drop)
    term <- unname(term)

  attributes(term) <- attr
  return(term)
}

#' Hook for custom symbols
#'
#' @param term
#' an R call or symbol/number. This is the expression to replace.
#'
#' @param display
#' an R call or symbol/number. This is shown instead of _term_.
#'
#' @return
#' TRUE on success
#'
#' @md
#'
#' @examples
#' hook(term=quote(t0), display=quote(subscript(t, 0)))
#' mathml(quote(t0))
#'
hook <- function(term=quote(t0), display=quote(subscript(t, 0)))
{
  r <- rolog::once(call("assert", call("math_hook", term, display)))
  if(isFALSE(r))
    return(FALSE)

  invisible(r)
}

#' Multiplication
#'
#' @name dot
#'
#' @param e1
#' numerator
#'
#' @param e2
#' denominator
#'
#' @return
#' e1 * e2
#'
dot <- function(e1, e2)
  e1 * e2

#' @rdname dot
#' @export
nodot <- dot

#' @rdname dot
#' @export
times <- dot

#' Division displayed as fraction
#'
#' @param e1
#' numerator
#'
#' @param e2
#' denominator
#'
#' @return
#' e1 / e2
#'
frac <- function(e1, e2)
  e1 / e2

#' Division displayed as fraction
#'
#' @param e1
#' numerator
#'
#' @param e2
#' denominator
#'
#' @return
#' e1 / e2
#'
over <- frac

#' Division displayed as large fraction
#'
#' @param e1
#' numerator
#'
#' @param e2
#' denominator
#'
#' @return
#' e1 / e2
#'
#' @md
#' @seealso [frac()], [over()]
dfrac <- frac

#' Return function body
#'
#' @param fname
#' not clear
#'
#' @param body
#' not clear
#'
#' @return
#' body
#'
fname <- function(fname, body)
{
  return(body)
}

#' Plus Minus
#'
#' @param x
#'
#' @param y
#'
#' @return c(x - y, x + y)
#'
'%+-%' <- function(x, y)
{
  return(c(x - y, x + y))
}


#'x dot y
#'
#' @param x
#'
#' @param y
#'
#' @return e1*e2
#'
'%.%' <- function(e1, e2)
e1*e2

#'x approx y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%~~%' <- function(x, y)
  isTRUE(all.equal(x, y))

#'x equiv y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%==%' <- function(x, y)
  x=y


#'x cong y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%=~%' <- function(x, y)
  x=y

#'x propto y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%prop%' <- function(x, y)
  x=y

#'x leftrightarrow y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%<->%' <- function(x, y)
  x=y

#'x rightarrow y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%->%' <- function(x, y)
  x=y

#'x leftarrow y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%<-%' <- function(x, y)
  x=y

#'x uparrow y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%up%' <- function(x, y)
  x=y

#'x downarrow y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%down%' <- function(x, y)
  x=y

#'x iff y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%<=>%' <- function(x, y)
  x=y

#'x rArr y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%<=%' <- function(x, y)
  x=y

#'x lArr y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%=>%' <- function(x, y)
  x=y

#'x uArr y
#'
#' @param x
#'
#' @param y
#'
#' @return e1=e2
#'
'%dblup%' <- function(x, y)
  x=y

#'x dArr y
#'
#' @param x
#'
#' @param y
#'
#' @return x=y
#'
'%dbldown%' <- function(x, y)
  x=y

#' denote
#' This is a function that allows the user to insert abbreviations in the formula,
#' explain them and make the needed computations
#'
#' @param abbr
#' Abbreviation used in the text to refer to the calculation, for example 's_p' for the pooled
#' variance.
#'
#' @param expr
#' Expression: calculations to be made in order to obtain the value to which the abbreviation
#' refers to.
#'
#' @param info
#' Information: Explanation of the formula used to provide the value of the abbreviation.
#' e.g. 'the pooled variance'
#'
#' @return expr
#'
denote <- function(abbr, expr, info)
  return(expr)

#' omit_left
#' This is a function that allows the user to highlight the mistakes,
#' in particular the omissions in the left-hand side of the expression
#'
#' @param expr
#'
#' @return substitute(expr)[[3]], e.g., b from a + b
#'
omit_left <- function(expr)
{
  # use third element of [-, A, B]
  eval(substitute(expr)[[3]])
}

#' omit_right
#' This is a function that allows the user to highlight the mistakes,
#' in particular the omissions in the right-hand side of the expression
#'
#' @param expr
#'
#' @return substitute(expr)[[2]], e.g., a from a + b
#'
omit_right <- function(expr)
{
  eval(substitute(expr)[[2]])
}

#' omit
#'
#' This is a function that allows the user to highlight the mistakes,
#' in particular the omission of an element from a list.
#'
#' @param expr
#'
#' @return NULL
#'
omit <- function(expr)
{
  NULL
}


#' add_left
#'
#' This is a function that allows the user to highlight the mistakes,
#' in particular the redundancies in the left-hand side of the expression.
#'
#' @param expr
#'
#' @return expr
#'
add_left <- function(expr)
{
  return(expr)
}

#' add_right
#'
#' This is a function that allows the user to highlight the mistakes,
#' in particular the redundancies in the right-hand side of the expression.
#'
#' @param expr
#'
#' @return expr
#'
add_right <- function(expr)
{
  return(expr)
}

#' add
#'
#' This is a function that allows the user to highlight the mistakes,
#' in particular an extra element in a list
#'
#' @param expr
#'
#' @return expr
#'
add <- function(expr)
{
  return(expr)
}

#' instead
#'
#' This is a function that allows the user to highlight the mistakes,
#' in particular adds a curly bracket under the wrong term and it provides the
#' correct solutions.
#'
#' @param inst
#'
#' @param of
#'
#' @return inst
#'
instead <- function(inst, of)
{
  return(inst)
}
