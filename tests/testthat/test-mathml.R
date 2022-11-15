tml <- function(term)
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
  q <- tml(as.symbol("atom"))
  expect_equal(q, "<math><mi>atom</mi></math>")
})

test_that("greek",
{
  q <- tml(as.symbol("pi"))
  expect_equal(q, "<math><mi>&pi;</mi></math>")
})
