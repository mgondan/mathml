tml <- function(term=quote((a + b)^2L == a^2L + 2L*a*b + b^2L))
{
  t = rolog::once(call("r2mathml", term, expression(X)))
  paste(t$X, collapse="")
}

test_that("strings",
{
  q <- tml("text")
  expect_equal(q, "<math><mtext>text</mtext></math>")
})

test_that("atoms",
{
  q <- tml(as.name(atom))
  expect_equal(q, "<math><mi>atom</mi></math>")
})
