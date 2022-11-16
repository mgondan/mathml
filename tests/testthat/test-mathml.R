library(testthat)

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

test_that("calligraphy",
{
  q <- tml(quote(K %in% cal(K)))
  expect_equal(q, "<math><mrow><mi>K</mi><mo>&isin;</mo><mi mathvariant=\"script\">K</mi></mrow></math>")
})

test_that("calligraphy",
{
  q <- tml(call("%in%", as.symbol("K"), call("cal", as.symbol("K"))))
  expect_equal(q, "<math><mrow><mi>K</mi><mo>&isin;</mo><mi mathvariant=\"script\">K</mi></mrow></math>")
})

test_that("absx",
{
  q <- tml(quote(abs(x)))
  expect_equal(q, "<math><mrow><mo>&vert;</mo><mi>x</mi><mo>&vert;</mo></mrow></math>")
})

test_that("signx",
{
  q <- tml(quote(sign(x)))
  expect_equal(q, "<math><mrow><mi>sgn</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
  })


test_that("sqrtalpha",
{
  q <- tml(quote(sqrt(alpha)))
  expect_equal(q, "<math><msqrt><mi>&alpha;</mi></msqrt></math>")
})


test_that("sin(pi/2)",
{
  q <- tml(quote(sin(pi/2L)))
  expect_equal(q, "<math><mrow><mi>sin</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&pi;</mi><mo>/</mo><mn>2</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("cos(pi/2)",
{
  q <- tml(quote(cos(pi/2L)))
  expect_equal(q, "<math><mrow><mi>cos</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&pi;</mi><mo>/</mo><mn>2</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("tan(pi/2)",
{
  q <- tml(quote(tan(pi/2L)))
  expect_equal(q, "<math><mrow><mi>tan</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&pi;</mi><mo>/</mo><mn>2</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("asin(1/2*sqrt2",
{
  q <- tml(quote(asin(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>sin</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("acos(1/2*sqrt2",
{
  q <- tml(quote(acos(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>cos</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("atan(1/2*sqrt2",
{
  q <- tml(quote(atan(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>tan</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("allx",
{
  q <- tml(quote(all(x)))
  expect_equal(q, "<math><mrow><mo>&ForAll;</mo><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("anyx",
{
  q <- tml(quote(any(x)))
  expect_equal(q, "<math><mrow><mo>&Exists;</mo><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("besselInu(x)",
{
  q <- tml(quote(besselI(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>I</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselKnu(x)",
{
  q <- tml(quote(besselK(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>K</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselKnu(x)",
{
  q <- tml(quote(besselK(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>K</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselJnu(x)",
{
  q <- tml(quote(besselJ(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>J</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselYnu(x)",
{
  q <- tml(quote(besselY(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>Y</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("digammax",
{
  q <- tml(quote(digamma(x)))
  expect_equal(q, "<math><mrow><mfrac><mi>d</mi><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></mrow></mrow></math>")
})

test_that("trigammax",
{
  q <- tml(quote(trigamma(x)))
  expect_equal(q, "<math><mrow><mfrac><msup><mi>d</mi><mn>2</mn></msup><msup><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow><mn>2</mn></msup></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></mrow></mrow></math>")
})

test_that("bincoef",
{
  q <- tml(quote(choose(n, k)))
  expect_equal(q, "<math><mrow><mo>(</mo><mfrac linethickness=\"0\"><mi>n</mi><mi>k</mi></mfrac><mo>)</mo></mrow></math>")
})

test_that("BinCoef",
{
  q <- tml(quote(choose(k=K, n=N)))
  expect_equal(q, "<math><mrow><mo>(</mo><mfrac linethickness=\"0\"><mrow><mi>k</mi><mo>=</mo><mi>K</mi></mrow><mrow><mi>n</mi><mo>=</mo><mi>N</mi></mrow></mfrac><mo>)</mo></mrow></math>")
 })

test_that("logbincoef",
{
  q <- tml(quote(lchoose(n, k)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>(</mo><mfrac linethickness=\"0\"><mi>n</mi><mi>k</mi></mfrac><mo>)</mo></mrow></mrow></math>")
})

test_that("factorial",
{
  q <- tml(quote(factorial(x)))
  expect_equal(q, "<math><mrow><mi>x</mi><mo>!</mo></mrow></math>")
})

test_that("lfactorial",
{
  q <- tml(quote(lfactorial(x)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mi>x</mi><mo>!</mo></mrow></mrow></math>")
})

test_that("AND",
{
  q <- tml(quote(A & B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&and;</mo><mi>B</mi></mrow></math>")
})

test_that("OR",
{
  q <- tml(quote(A | B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&or;</mo><mi>B</mi></mrow></math>")
})

f <- as.function(alist(a=, b=2, a+b))
test_that("function",
{
  q <- tml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("ifelse",
{
  q <- tml(quote(ifelse(a > b, a, b)))
  expect_equal(q,"<math><mrow><mo>{</mo><mtable columnalign=\"left\"><mtr><mi>a</mi><mrow><mtext>if</mtext><mspace width=\"thinmathspace\"></mspace><mrow><mi>a</mi><mo>&gt;</mo><mi>b</mi></mrow></mrow></mtr><mtr><mi>b</mi><mtext>otherwise</mtext></mtr></mtable></mrow></math>")
})






