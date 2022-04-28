.onAttach <- function(libname, pkgname)
{
  library(rolog, quietly=TRUE)
  consult(system.file("pl/mathml.pl", package=pkgname))
}

mathml = function(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L))
{
  t = once(call("r2mathml", term, expression(X)))
  cat(paste(t$X, collapse=""))
}
