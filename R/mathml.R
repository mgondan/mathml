.onAttach <- function(libname, pkgname)
{
  library(rolog)
  consult(system.file("pl/mathml.pl", package=pkgname))
}

mathml = function(term)
{
  t = once(call("r2mathml", term, expression(X)))
  cat(paste(t$X, collapse=""))
}
