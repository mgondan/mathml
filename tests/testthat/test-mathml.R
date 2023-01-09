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

test_that("absgreek",
{
  q <- tml(quote(abs(alpha)))
  expect_equal(q, "<math><mrow><mo>&vert;</mo><mi>&alpha;</mi><mo>&vert;</mo></mrow></math>")
})

test_that("signx",
{
  q <- tml(quote(sign(x)))
  expect_equal(q, "<math><mrow><mi>sgn</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("signgreek",
{
  q <- tml(quote(sign(alpha)))
  expect_equal(q, "<math><mrow><mi>sgn</mi><mo>&af;</mo><mi>&alpha;</mi></mrow></math>")
})

test_that("sqrtnum",
{
  q <- tml(quote(sqrt(2L)))
  expect_equal(q, "<math><msqrt><mn>2</mn></msqrt></math>")
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

test_that("asin(1/2*sqrt2)",
{
  q <- tml(quote(asin(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>sin</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("acos(1/2*sqrt2)",
{
  q <- tml(quote(acos(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>cos</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("atan(1/2*sqrt2)",
{
  q <- tml(quote(atan(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>tan</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("atan(1/2*sqrt3)",
{
  q <- tml(quote(atan(1L/2L * sqrt(3L))))
  expect_equal(q, "<math><mrow><msup><mi>tan</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>3</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("singreek(num)",
{
  q <- tml(quote(sinpi(2L)))
  expect_equal(q, "<math><mrow><mi>sin</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("cosgreek(num)",
{
  q <- tml(quote(cospi(2L)))
  expect_equal(q, "<math><mrow><mi>cos</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("tangreek(num)",
{
  q <- tml(quote(tanpi(2L)))
  expect_equal(q, "<math><mrow><mi>tan</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("sinhx",
{
  q <- tml(quote(sinh(x)))
  expect_equal(q, "<math><mrow><mi>sinh</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("coshx",
{
  q <- tml(quote(cosh(x)))
  expect_equal(q, "<math><mrow><mi>cosh</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("tanhx",
{
  q <- tml(quote(tanh(x)))
  expect_equal(q, "<math><mrow><mi>tanh</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("asinhx",
{
  q <- tml(quote(asinh(x)))
  expect_equal(q, "<math><mrow><msup><mi>sinh</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("acoshx",
{
  q <- tml(quote(acosh(x)))
  expect_equal(q, "<math><mrow><msup><mi>cosh</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("atanhx",
{
  q <- tml(quote(atanh(x)))
  expect_equal(q, "<math><mrow><msup><mi>tanh</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
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

test_that("beta(a,b)",
{
  q <- tml(quote(beta(a, b)))
  expect_equal(q, "<math><mrow><mi>B</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>a</mi><mo>,</mo><mi>b</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("lbeta(a,b)",
{
  q <- tml(quote(lbeta(a, b)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>B</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>a</mi><mo>,</mo><mi>b</mi></mrow><mo>)</mo></mrow></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("gammax",
{
  q <- tml(quote(gamma(x)))
  expect_equal(q, "<math><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("lgammax",
{
  q <- tml(quote(lgamma(x)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></math>")
})

test_that("digammax",
{
  q <- tml(quote(digamma(x)))
  expect_equal(q, "<math><mrow><mfrac><mi>d</mi><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
})

test_that("psigammax",
{
  q <- tml(quote(psigamma(x, deriv=psi)))
  expect_equal(q, "<math><mrow><mfrac><msup><mi>d</mi><mrow><mrow><mo>(</mo><mrow><mi>deriv</mi><mo>=</mo><mi>&psi;</mi></mrow><mo>)</mo></mrow><mo>+</mo><mn>2</mn></mrow></msup><msup><mrow><mo>(</mo><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow><mo>)</mo></mrow><mrow><mrow><mo>(</mo><mrow><mi>deriv</mi><mo>=</mo><mi>&psi;</mi></mrow><mo>)</mo></mrow><mo>+</mo><mn>2</mn></mrow></msup></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
})

test_that("alogb",
{
  q <- tml(quote(a*log(b)))
  expect_equal(q, "<math><mrow><mi>a</mi><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mi>b</mi></mrow></mrow></math>")
})

test_that("bloga",
{
  q <- tml(quote(log(b)*a))
  expect_equal(q, "<math><mrow><mrow><mi>log</mi><mo>&af;</mo><mi>b</mi></mrow><mo>&sdot;</mo><mi>a</mi></mrow></math>")
})

test_that("sumab",
{
  q <- tml(quote(a*sum(b)))
  expect_equal(q, "<math><mrow><mi>a</mi><mo>&sdot;</mo><mrow><mo>&sum;</mo><mo>&af;</mo><mi>b</mi></mrow></mrow></math>")
})

test_that("sumba",
{
  q <- tml(quote(sum(b)*a))
  expect_equal(q, "<math><mrow><mrow><mo>&sum;</mo><mo>&af;</mo><mi>b</mi></mrow><mo>&sdot;</mo><mi>a</mi></mrow></math>")
})

test_that("digammax",
{
  q <- tml(quote(digamma(x)))
  expect_equal(q, "<math><mrow><mfrac><mi>d</mi><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
})

test_that("trigammax",
{
  q <- tml(quote(trigamma(x)))
  expect_equal(q, "<math><mrow><mfrac><msup><mi>d</mi><mn>2</mn></msup><msup><mrow><mo>(</mo><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow><mo>)</mo></mrow><mn>2</mn></msup></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
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
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>!</mo></mrow><mo>)</mo></mrow></mrow></math>")
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

test_that("!A",
{
  q <- tml(quote(!A))
  expect_equal(q, "<math><mrow><mo>&not;</mo><mi>A</mi></mrow></math>")
})

test_that("xor",
{
  q <- tml(quote(xor(A, B)))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&veebar;</mo><mi>B</mi></mrow></math>")
})

test_that("expnumgreeki",
{
  q <- tml(quote(exp(2L*pi*i)))
  expect_equal(q, "<math><mrow><mi>exp</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>&sdot;</mo><mi>i</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("expm1numgreeki",
{
  q <- tml(quote(expm1(2L*pi*i)))
  expect_equal(q, "<math><mrow><mrow><mi>exp</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>&sdot;</mo><mi>i</mi></mrow><mo>)</mo></mrow></mrow><mo>-</mo><mn>1</mn></mrow></math>")
})

test_that("logx",
{
  q <- tml(quote(log(x)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("log10x",
{
  q <- tml(quote(log10(x)))
  expect_equal(q, "<math><mrow><msub><mi>log</mi><mn>10</mn></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("log2x",
{
  q <- tml(quote(log2(x)))
  expect_equal(q, "<math><mrow><msub><mi>log</mi><mn>2</mn></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("logex",
{
  q <- tml(quote(logb(x, e)))
  expect_equal(q, "<math><mrow><msub><mi>log</mi><mi>e</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("log1px",
{
  q <- tml(quote(log1p(x)))
  expect_equal(q, "<math><mrow><mn>1</mn><mo>+</mo><mrow><mi>log</mi><mo>&af;</mo><mi>x</mi></mrow></mrow></math>")
})

test_that("ceilingpi",
{
  q <- tml(quote(ceiling(pi)))
  expect_equal(q, "<math><mrow><mo>&lceil;</mo><mi>&pi;</mi><mo>&rceil;</mo></mrow></math>")
})

test_that("floorpi",
{
  q <- tml(quote(floor(pi)))
  expect_equal(q, "<math><mrow><mo>&lfloor;</mo><mi>&pi;</mi><mo>&rfloor;</mo></mrow></math>")
})

f <- as.function(alist(a=, b=2, a+b))
test_that("function",
{
  q <- tml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

f <- function(x) sin(x)
test_that("functionsin",
{
  q <- tml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

f <- function(x) { sin(x) ; tan(x)}
test_that("functionsintan",
{
  q <- tml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

f <- function(x) {sin(x) ; cos(x)}
test_that("functionsincos",
{
  q <- tml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("identical",
{
  q <- tml(quote(identical(1L, 2L)))
  expect_equal(q,"<math><mrow><mn>1</mn><mo>=</mo><mn>2</mn></mrow></math>")
 })

test_that("ifelse",
{
  q <- tml(quote(ifelse(a > b, a, b)))
  expect_equal(q,"<math><mrow><mo>{</mo><mtable columnalign=\"left\"><mtr><mi>a</mi><mrow><mtext>if</mtext><mspace width=\"thinmathspace\"></mspace><mrow><mi>a</mi><mo>&gt;</mo><mi>b</mi></mrow></mrow></mtr><mtr><mi>b</mi><mtext>otherwise</mtext></mtr></mtable></mrow></math>")
})

test_that("in",
{
  q <- tml(quote(a %in% A))
  expect_equal(q,"<math><mrow><mi>a</mi><mo>&isin;</mo><mi>A</mi></mrow></math>")
})

test_that("intersect",
{
  q <- tml(quote(intersect(A, B)))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>&cap;</mo><mi>B</mi></mrow></math>")
})

test_that("union",
{
  q <- tml(quote(union(A, B)))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>&cup;</mo><mi>B</mi></mrow></math>")
})

test_that("null",
{
  q <- tml(quote(is.null(intersect(A, B))))
  expect_equal(q,"<math><mrow><mrow><mi>A</mi><mo>&cap;</mo><mi>B</mi></mrow><mo>=</mo><mi>&empty;</mi></mrow></math>")
})

test_that("column",
{
q <- tml(quote(a:b))
            expect_equal(q,"<math><mrow><mi>a</mi><mo>&#58;</mo><mi>b</mi></mrow></math>")
})

test_that("diff",
{
  q <- tml(quote(setdiff(A, B)))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>-</mo><mi>B</mi></mrow></math>")
})

test_that("x",
{
  q <- tml(quote(X %x% Y))
  expect_equal(q,"<math><mrow><mi>X</mi><mo>&CircleTimes;</mo><mi>Y</mi></mrow></math>")
})

test_that("expr",
{
  q <- tml(quote(a*b + c/d - e %*% f))
  expect_equal(q,"<math><mrow><mrow><mrow><mi>a</mi><mo>&#x2062;</mo><mi>b</mi></mrow><mo>+</mo><mrow><mi>c</mi><mo>/</mo><mi>d</mi></mrow></mrow><mo>-</mo><mrow><mi>e</mi><mo>&times;</mo><mi>f</mi></mrow></mrow></math>")
})

test_that("crossprod",
{
  q <- tml(quote(crossprod(A, B)))
  expect_equal(q,"<math><mrow><msup><mi>A</mi><mtext>T</mtext></msup><mo>&times;</mo><mi>B</mi></mrow></math>")
})

test_that("tcrossprod",
{
  q <- tml(quote(tcrossprod(A, B)))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&times;</mo><msup><mi>B</mi><mtext>T</mtext></msup></mrow></math>")
})

test_that("tilde",
{
  q <- tml(quote(A ~ B))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>&Tilde;</mo><mi>B</mi></mrow></math>")
})

test_that("expx",
{
  q <- tml(quote(exp(x)))
  expect_equal(q,"<math><mrow><mi>exp</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("expm1x",
{
  q <- tml(quote(expm1(x)))
  expect_equal(q,"<math><mrow><mrow><mi>exp</mi><mo>&af;</mo><mi>x</mi></mrow><mo>-</mo><mn>1</mn></mrow></math>")
})

test_that("dbinom",
{
  q <- tml(quote(dbinom(k, N, pi)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Bi</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>=</mo><mi>k</mi></mrow><mo>;</mo><mrow><mi>N</mi><mo>,</mo><mi>&pi;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("pbinom",
{
  q <- tml(quote(pbinom(k, N, pi)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Bi</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mrow><mi>N</mi><mo>,</mo><mi>&pi;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qbinom",
{
  q <- tml(quote(qbinom(k, N, pi)))
  expect_equal(q,"<math><mrow><msub><mtext>arg min</mtext><mi>k</mi></msub><mo>&af;</mo><mrow><mo>[</mo><mrow><mrow><msub><mi>P</mi><mtext>Bi</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mrow><mi>N</mi><mo>,</mo><mi>&pi;</mi></mrow></mrow><mo>)</mo></mrow></mrow><mo>&gt;</mo><mi>k</mi></mrow><mo>]</mo></mrow></mrow></math>")  
})

test_that("dpois",
{
  q <- tml(quote(dpois(k, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Po</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>=</mo><mi>k</mi></mrow><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("ppois",
{
  q <- tml(quote(ppois(k, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Po</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qpois",
{
  q <- tml(quote(qpois(k, lambda)))
  expect_equal(q,"<math><mrow><msub><mtext>arg max</mtext><mi>k</mi></msub><mo>&af;</mo><mrow><mo>[</mo><mrow><mrow><msub><mi>P</mi><mtext>Po</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow><mo>&gt;</mo><mi>k</mi></mrow><mo>]</mo></mrow></mrow></math>")  
})

test_that("dexp",
{
  q <- tml(quote(dexp(x, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>f</mi><mtext>Exp</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("pexp",
{
  q <- tml(quote(pexp(x, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>F</mi><mtext>Exp</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qexp",
{
  q <- tml(quote(qexp(x, lambda)))
  expect_equal(q,"<math><mrow><msubsup><mi>F</mi><mtext>Exp</mtext><mrow><mo>-</mo><mn>1</mn></mrow></msubsup><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("dnorm",
{
  q <- tml(quote(dnorm(x, mu, sigma)))
  expect_equal(q,"<math><mrow><mi>&phi;</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mrow><mi>&mu;</mi><mo>,</mo><mi>&sigma;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("pnorm",
{
  q <- tml(quote(pnorm(x, mu, sigma)))
  expect_equal(q,"<math><mrow><mi>&Phi;</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mrow><mi>&mu;</mi><mo>,</mo><mi>&sigma;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qnorm",
{
  q <- tml(quote(qnorm(alpha/2)))
  expect_equal(q,"<math><mrow><msup><mi>&Phi;</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&alpha;</mi><mo>/</mo><mn>2.00</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("1mpchisq",
{
  q <- tml(quote(1 - pchisq(x, 1)))
  expect_equal(q,"<math><mrow><mn>1.00</mn><mo>-</mo><mrow><msub><mi>F</mi><mrow><msup><mi>&chi;</mi><mn>2</mn></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>1.00</mn><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow><mo>)</mo></mrow></mrow></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></mrow></math>")
})

test_that("pchisq1malpha",
{
  q <- tml(quote(qchisq(1 - alpha, 1)))
  expect_equal(q,"<math><mrow><msubsup><mi>F</mi><mrow><msup><mi>&chi;</mi><mn>2</mn></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>1.00</mn><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow><mo>)</mo></mrow></mrow><mrow><mo>-</mo><mn>1</mn></mrow></msubsup><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>1.00</mn><mo>-</mo><mi>&alpha;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("ptnm1",
{
  q <- tml(quote(pt(t, N - 1)))
  expect_equal(q, "<math><mrow><mi>P</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>T</mi><mo>&le;</mo><mi>t</mi></mrow><mo>;</mo><mrow><mrow><mi>N</mi><mo>-</mo><mn>1.00</mn></mrow><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("qtalpha2nm1",
{
  q <- tml(quote(qt(alpha/2, N - 1)))
  expect_equal(q, "<math><mrow><msub><mi>T</mi><mrow><mi>&alpha;</mi><mo>/</mo><mn>2.00</mn></mrow></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>N</mi><mo>-</mo><mn>1.00</mn></mrow><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow><mo>)</mo></mrow></mrow></math>")
})



