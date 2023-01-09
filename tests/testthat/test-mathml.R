library(testthat)

test_that("strings",
{
  q <- mathml("text")
  expect_equal(q, "<math><mtext>text</mtext></math>")
})

test_that("atoms",
{
  q <- mathml(as.symbol("atom"))
  expect_equal(q, "<math><mi>atom</mi></math>")
})

test_that("greek",
{
  q <- mathml(as.symbol("pi"))
  expect_equal(q, "<math><mi>&pi;</mi></math>")
})

test_that("calligraphy",
{
  q <- mathml(quote(K %in% cal(K)))
  expect_equal(q, "<math><mrow><mi>K</mi><mo>&isin;</mo><mi mathvariant=\"script\">K</mi></mrow></math>")
})

test_that("calligraphy",
{
  q <- mathml(call("%in%", as.symbol("K"), call("cal", as.symbol("K"))))
  expect_equal(q, "<math><mrow><mi>K</mi><mo>&isin;</mo><mi mathvariant=\"script\">K</mi></mrow></math>")
})

test_that("absx",
{
  q <- mathml(quote(abs(x)))
  expect_equal(q, "<math><mrow><mo>&vert;</mo><mi>x</mi><mo>&vert;</mo></mrow></math>")
})

test_that("absgreek",
{
  q <- mathml(quote(abs(alpha)))
  expect_equal(q, "<math><mrow><mo>&vert;</mo><mi>&alpha;</mi><mo>&vert;</mo></mrow></math>")
})

test_that("signx",
{
  q <- mathml(quote(sign(x)))
  expect_equal(q, "<math><mrow><mi>sgn</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("signgreek",
{
  q <- mathml(quote(sign(alpha)))
  expect_equal(q, "<math><mrow><mi>sgn</mi><mo>&af;</mo><mi>&alpha;</mi></mrow></math>")
})

test_that("sqrtnum",
{
  q <- mathml(quote(sqrt(2L)))
  expect_equal(q, "<math><msqrt><mn>2</mn></msqrt></math>")
})

test_that("sqrtalpha",
{
  q <- mathml(quote(sqrt(alpha)))
  expect_equal(q, "<math><msqrt><mi>&alpha;</mi></msqrt></math>")
})

test_that("sin(pi/2)",
{
  q <- mathml(quote(sin(pi/2L)))
  expect_equal(q, "<math><mrow><mi>sin</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&pi;</mi><mo>/</mo><mn>2</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("cos(pi/2)",
{
  q <- mathml(quote(cos(pi/2L)))
  expect_equal(q, "<math><mrow><mi>cos</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&pi;</mi><mo>/</mo><mn>2</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("tan(pi/2)",
{
  q <- mathml(quote(tan(pi/2L)))
  expect_equal(q, "<math><mrow><mi>tan</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&pi;</mi><mo>/</mo><mn>2</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("asin(1/2*sqrt2)",
{
  q <- mathml(quote(asin(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>sin</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("acos(1/2*sqrt2)",
{
  q <- mathml(quote(acos(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>cos</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("atan(1/2*sqrt2)",
{
  q <- mathml(quote(atan(1L/2L * sqrt(2L))))
  expect_equal(q, "<math><mrow><msup><mi>tan</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>2</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("atan(1/2*sqrt3)",
{
  q <- mathml(quote(atan(1L/2L * sqrt(3L))))
  expect_equal(q, "<math><mrow><msup><mi>tan</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>1</mn><mo>/</mo><mn>2</mn></mrow><mo>&sdot;</mo><msqrt><mn>3</mn></msqrt></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("singreek(num)",
{
  q <- mathml(quote(sinpi(2L)))
  expect_equal(q, "<math><mrow><mi>sin</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("cosgreek(num)",
{
  q <- mathml(quote(cospi(2L)))
  expect_equal(q, "<math><mrow><mi>cos</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("tangreek(num)",
{
  q <- mathml(quote(tanpi(2L)))
  expect_equal(q, "<math><mrow><mi>tan</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("sinhx",
{
  q <- mathml(quote(sinh(x)))
  expect_equal(q, "<math><mrow><mi>sinh</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("coshx",
{
  q <- mathml(quote(cosh(x)))
  expect_equal(q, "<math><mrow><mi>cosh</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("tanhx",
{
  q <- mathml(quote(tanh(x)))
  expect_equal(q, "<math><mrow><mi>tanh</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("asinhx",
{
  q <- mathml(quote(asinh(x)))
  expect_equal(q, "<math><mrow><msup><mi>sinh</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("acoshx",
{
  q <- mathml(quote(acosh(x)))
  expect_equal(q, "<math><mrow><msup><mi>cosh</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("atanhx",
{
  q <- mathml(quote(atanh(x)))
  expect_equal(q, "<math><mrow><msup><mi>tanh</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("allx",
{
  q <- mathml(quote(all(x)))
  expect_equal(q, "<math><mrow><mo>&ForAll;</mo><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("anyx",
{
  q <- mathml(quote(any(x)))
  expect_equal(q, "<math><mrow><mo>&Exists;</mo><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("besselInu(x)",
{
  q <- mathml(quote(besselI(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>I</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselKnu(x)",
{
  q <- mathml(quote(besselK(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>K</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselKnu(x)",
{
  q <- mathml(quote(besselK(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>K</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselJnu(x)",
{
  q <- mathml(quote(besselJ(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>J</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("besselYnu(x)",
{
  q <- mathml(quote(besselY(x, nu)))
  expect_equal(q, "<math><mrow><msub><mi>Y</mi><mi>&nu;</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("beta(a,b)",
{
  q <- mathml(quote(beta(a, b)))
  expect_equal(q, "<math><mrow><mi>B</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>a</mi><mo>,</mo><mi>b</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("lbeta(a,b)",
{
  q <- mathml(quote(lbeta(a, b)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>B</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>a</mi><mo>,</mo><mi>b</mi></mrow><mo>)</mo></mrow></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("gammax",
{
  q <- mathml(quote(gamma(x)))
  expect_equal(q, "<math><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("lgammax",
{
  q <- mathml(quote(lgamma(x)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></math>")
})

test_that("digammax",
{
  q <- mathml(quote(digamma(x)))
  expect_equal(q, "<math><mrow><mfrac><mi>d</mi><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
})

test_that("psigammax",
{
  q <- mathml(quote(psigamma(x, deriv=psi)))
  expect_equal(q, "<math><mrow><mfrac><msup><mi>d</mi><mrow><mrow><mo>(</mo><mrow><mi>deriv</mi><mo>=</mo><mi>&psi;</mi></mrow><mo>)</mo></mrow><mo>+</mo><mn>2</mn></mrow></msup><msup><mrow><mo>(</mo><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow><mo>)</mo></mrow><mrow><mrow><mo>(</mo><mrow><mi>deriv</mi><mo>=</mo><mi>&psi;</mi></mrow><mo>)</mo></mrow><mo>+</mo><mn>2</mn></mrow></msup></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
})

test_that("alogb",
{
  q <- mathml(quote(a*log(b)))
  expect_equal(q, "<math><mrow><mi>a</mi><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mi>b</mi></mrow></mrow></math>")
})

test_that("bloga",
{
  q <- mathml(quote(log(b)*a))
  expect_equal(q, "<math><mrow><mrow><mi>log</mi><mo>&af;</mo><mi>b</mi></mrow><mo>&sdot;</mo><mi>a</mi></mrow></math>")
})

test_that("sumab",
{
  q <- mathml(quote(a*sum(b)))
  expect_equal(q, "<math><mrow><mi>a</mi><mo>&sdot;</mo><mrow><mo>&sum;</mo><mo>&af;</mo><mi>b</mi></mrow></mrow></math>")
})

test_that("sumba",
{
  q <- mathml(quote(sum(b)*a))
  expect_equal(q, "<math><mrow><mrow><mo>&sum;</mo><mo>&af;</mo><mi>b</mi></mrow><mo>&sdot;</mo><mi>a</mi></mrow></math>")
})

test_that("digammax",
{
  q <- mathml(quote(digamma(x)))
  expect_equal(q, "<math><mrow><mfrac><mi>d</mi><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
})

test_that("trigammax",
{
  q <- mathml(quote(trigamma(x)))
  expect_equal(q, "<math><mrow><mfrac><msup><mi>d</mi><mn>2</mn></msup><msup><mrow><mo>(</mo><mrow><mi>d</mi><mo>&#x2062;</mo><mi>x</mi></mrow><mo>)</mo></mrow><mn>2</mn></msup></mfrac><mo>&sdot;</mo><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>[</mo><mrow><mi>&Gamma;</mi><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow><mo>]</mo></mrow></mrow></mrow></math>")
})

test_that("bincoef",
{
  q <- mathml(quote(choose(n, k)))
  expect_equal(q, "<math><mrow><mo>(</mo><mfrac linethickness=\"0\"><mi>n</mi><mi>k</mi></mfrac><mo>)</mo></mrow></math>")
})

test_that("BinCoef",
{
  q <- mathml(quote(choose(k=K, n=N)))
  expect_equal(q, "<math><mrow><mo>(</mo><mfrac linethickness=\"0\"><mrow><mi>k</mi><mo>=</mo><mi>K</mi></mrow><mrow><mi>n</mi><mo>=</mo><mi>N</mi></mrow></mfrac><mo>)</mo></mrow></math>")
 })

test_that("logbincoef",
{
  q <- mathml(quote(lchoose(n, k)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>(</mo><mfrac linethickness=\"0\"><mi>n</mi><mi>k</mi></mfrac><mo>)</mo></mrow></mrow></math>")
})

test_that("factorial",
{
  q <- mathml(quote(factorial(x)))
  expect_equal(q, "<math><mrow><mi>x</mi><mo>!</mo></mrow></math>")
})

test_that("lfactorial",
{
  q <- mathml(quote(lfactorial(x)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>!</mo></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("AND",
{
  q <- mathml(quote(A & B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&and;</mo><mi>B</mi></mrow></math>")
})

test_that("OR",
{
  q <- mathml(quote(A | B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&or;</mo><mi>B</mi></mrow></math>")
})

test_that("!A",
{
  q <- mathml(quote(!A))
  expect_equal(q, "<math><mrow><mo>&not;</mo><mi>A</mi></mrow></math>")
})

test_that("xor",
{
  q <- mathml(quote(xor(A, B)))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&veebar;</mo><mi>B</mi></mrow></math>")
})

test_that("expnumgreeki",
{
  q <- mathml(quote(exp(2L*pi*i)))
  expect_equal(q, "<math><mrow><mi>exp</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>&sdot;</mo><mi>i</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("expm1numgreeki",
{
  q <- mathml(quote(expm1(2L*pi*i)))
  expect_equal(q, "<math><mrow><mrow><mi>exp</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mn>2</mn><mo>&#x2062;</mo><mi>&pi;</mi></mrow><mo>&sdot;</mo><mi>i</mi></mrow><mo>)</mo></mrow></mrow><mo>-</mo><mn>1</mn></mrow></math>")
})

test_that("logx",
{
  q <- mathml(quote(log(x)))
  expect_equal(q, "<math><mrow><mi>log</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("log10x",
{
  q <- mathml(quote(log10(x)))
  expect_equal(q, "<math><mrow><msub><mi>log</mi><mn>10</mn></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("log2x",
{
  q <- mathml(quote(log2(x)))
  expect_equal(q, "<math><mrow><msub><mi>log</mi><mn>2</mn></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("logex",
{
  q <- mathml(quote(logb(x, e)))
  expect_equal(q, "<math><mrow><msub><mi>log</mi><mi>e</mi></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("log1px",
{
  q <- mathml(quote(log1p(x)))
  expect_equal(q, "<math><mrow><mn>1</mn><mo>+</mo><mrow><mi>log</mi><mo>&af;</mo><mi>x</mi></mrow></mrow></math>")
})

test_that("ceilingpi",
{
  q <- mathml(quote(ceiling(pi)))
  expect_equal(q, "<math><mrow><mo>&lceil;</mo><mi>&pi;</mi><mo>&rceil;</mo></mrow></math>")
})

test_that("floorpi",
{
  q <- mathml(quote(floor(pi)))
  expect_equal(q, "<math><mrow><mo>&lfloor;</mo><mi>&pi;</mi><mo>&rfloor;</mo></mrow></math>")
})

f <- as.function(alist(a=, b=2, a+b))
test_that("function",
{
  q <- mathml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

f <- function(x) sin(x)
test_that("functionsin",
{
  q <- mathml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

f <- function(x) { sin(x) ; tan(x)}
test_that("functionsintan",
{
  q <- mathml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

f <- function(x) {sin(x) ; cos(x)}
test_that("functionsincos",
{
  q <- mathml(quote(canonical(f)))
  expect_equal(q,"<math><mrow><mi>canonical</mi><mo>&af;</mo><mrow><mo>(</mo><mi>f</mi><mo>)</mo></mrow></mrow></math>")
})

test_that("identical",
{
  q <- mathml(quote(identical(1L, 2L)))
  expect_equal(q,"<math><mrow><mn>1</mn><mo>=</mo><mn>2</mn></mrow></math>")
 })

test_that("ifelse",
{
  q <- mathml(quote(ifelse(a > b, a, b)))
  expect_equal(q,"<math><mrow><mo>{</mo><mtable columnalign=\"left\"><mtr><mi>a</mi><mrow><mtext>if</mtext><mspace width=\"thinmathspace\"></mspace><mrow><mi>a</mi><mo>&gt;</mo><mi>b</mi></mrow></mrow></mtr><mtr><mi>b</mi><mtext>otherwise</mtext></mtr></mtable></mrow></math>")
})

test_that("in",
{
  q <- mathml(quote(a %in% A))
  expect_equal(q,"<math><mrow><mi>a</mi><mo>&isin;</mo><mi>A</mi></mrow></math>")
})

test_that("intersect",
{
  q <- mathml(quote(intersect(A, B)))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>&cap;</mo><mi>B</mi></mrow></math>")
})

test_that("union",
{
  q <- mathml(quote(union(A, B)))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>&cup;</mo><mi>B</mi></mrow></math>")
})

test_that("null",
{
  q <- mathml(quote(is.null(intersect(A, B))))
  expect_equal(q,"<math><mrow><mrow><mi>A</mi><mo>&cap;</mo><mi>B</mi></mrow><mo>=</mo><mi>&empty;</mi></mrow></math>")
})

test_that("column",
{
q <- mathml(quote(a:b))
            expect_equal(q,"<math><mrow><mi>a</mi><mo>&#58;</mo><mi>b</mi></mrow></math>")
})

test_that("diff",
{
  q <- mathml(quote(setdiff(A, B)))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>-</mo><mi>B</mi></mrow></math>")
})

test_that("x",
{
  q <- mathml(quote(X %x% Y))
  expect_equal(q,"<math><mrow><mi>X</mi><mo>&CircleTimes;</mo><mi>Y</mi></mrow></math>")
})

test_that("expr",
{
  q <- mathml(quote(a*b + c/d - e %*% f))
  expect_equal(q,"<math><mrow><mrow><mrow><mi>a</mi><mo>&#x2062;</mo><mi>b</mi></mrow><mo>+</mo><mrow><mi>c</mi><mo>/</mo><mi>d</mi></mrow></mrow><mo>-</mo><mrow><mi>e</mi><mo>&times;</mo><mi>f</mi></mrow></mrow></math>")
})

test_that("crossprod",
{
  q <- mathml(quote(crossprod(A, B)))
  expect_equal(q,"<math><mrow><msup><mi>A</mi><mtext>T</mtext></msup><mo>&times;</mo><mi>B</mi></mrow></math>")
})

test_that("tcrossprod",
{
  q <- mathml(quote(tcrossprod(A, B)))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&times;</mo><msup><mi>B</mi><mtext>T</mtext></msup></mrow></math>")
})

test_that("tilde",
{
  q <- mathml(quote(A ~ B))
  expect_equal(q,"<math><mrow><mi>A</mi><mo>&Tilde;</mo><mi>B</mi></mrow></math>")
})

test_that("expx",
{
  q <- mathml(quote(exp(x)))
  expect_equal(q,"<math><mrow><mi>exp</mi><mo>&af;</mo><mi>x</mi></mrow></math>")
})

test_that("expm1x",
{
  q <- mathml(quote(expm1(x)))
  expect_equal(q,"<math><mrow><mrow><mi>exp</mi><mo>&af;</mo><mi>x</mi></mrow><mo>-</mo><mn>1</mn></mrow></math>")
})

test_that("dbinom",
{
  q <- mathml(quote(dbinom(k, N, pi)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Bi</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>=</mo><mi>k</mi></mrow><mo>;</mo><mrow><mi>N</mi><mo>,</mo><mi>&pi;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("pbinom",
{
  q <- mathml(quote(pbinom(k, N, pi)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Bi</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mrow><mi>N</mi><mo>,</mo><mi>&pi;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qbinom",
{
  q <- mathml(quote(qbinom(k, N, pi)))
  expect_equal(q,"<math><mrow><msub><mtext>arg min</mtext><mi>k</mi></msub><mo>&af;</mo><mrow><mo>[</mo><mrow><mrow><msub><mi>P</mi><mtext>Bi</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mrow><mi>N</mi><mo>,</mo><mi>&pi;</mi></mrow></mrow><mo>)</mo></mrow></mrow><mo>&gt;</mo><mi>k</mi></mrow><mo>]</mo></mrow></mrow></math>")  
})

test_that("dpois",
{
  q <- mathml(quote(dpois(k, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Po</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>=</mo><mi>k</mi></mrow><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("ppois",
{
  q <- mathml(quote(ppois(k, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>P</mi><mtext>Po</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qpois",
{
  q <- mathml(quote(qpois(k, lambda)))
  expect_equal(q,"<math><mrow><msub><mtext>arg max</mtext><mi>k</mi></msub><mo>&af;</mo><mrow><mo>[</mo><mrow><mrow><msub><mi>P</mi><mtext>Po</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>X</mi><mo>&le;</mo><mi>k</mi></mrow><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow><mo>&gt;</mo><mi>k</mi></mrow><mo>]</mo></mrow></mrow></math>")  
})

test_that("dexp",
{
  q <- mathml(quote(dexp(x, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>f</mi><mtext>Exp</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("pexp",
{
  q <- mathml(quote(pexp(x, lambda)))
  expect_equal(q,"<math><mrow><msub><mi>F</mi><mtext>Exp</mtext></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qexp",
{
  q <- mathml(quote(qexp(x, lambda)))
  expect_equal(q,"<math><mrow><msubsup><mi>F</mi><mtext>Exp</mtext><mrow><mo>-</mo><mn>1</mn></mrow></msubsup><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mi>&lambda;</mi></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("dnorm",
{
  q <- mathml(quote(dnorm(x, mu, sigma)))
  expect_equal(q,"<math><mrow><mi>&phi;</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mrow><mi>&mu;</mi><mo>,</mo><mi>&sigma;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("pnorm",
{
  q <- mathml(quote(pnorm(x, mu, sigma)))
  expect_equal(q,"<math><mrow><mi>&Phi;</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>x</mi><mo>;</mo><mrow><mi>&mu;</mi><mo>,</mo><mi>&sigma;</mi></mrow></mrow><mo>)</mo></mrow></mrow></math>")  
})

test_that("qnorm",
{
  q <- mathml(quote(qnorm(alpha/2)))
  expect_equal(q,"<math><mrow><msup><mi>&Phi;</mi><mrow><mo>-</mo><mn>1</mn></mrow></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mi>&alpha;</mi><mo>/</mo><mn>2.00</mn></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("1mpchisq",
{
  q <- mathml(quote(1 - pchisq(x, 1)))
  expect_equal(q,"<math><mrow><mn>1.00</mn><mo>-</mo><mrow><msub><mi>F</mi><mrow><msup><mi>&chi;</mi><mn>2</mn></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>1.00</mn><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow><mo>)</mo></mrow></mrow></msub><mo>&af;</mo><mrow><mo>(</mo><mi>x</mi><mo>)</mo></mrow></mrow></mrow></math>")
})

test_that("pchisq1malpha",
{
  q <- mathml(quote(qchisq(1 - alpha, 1)))
  expect_equal(q,"<math><mrow><msubsup><mi>F</mi><mrow><msup><mi>&chi;</mi><mn>2</mn></msup><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>1.00</mn><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow><mo>)</mo></mrow></mrow><mrow><mo>-</mo><mn>1</mn></mrow></msubsup><mo>&af;</mo><mrow><mo>(</mo><mrow><mn>1.00</mn><mo>-</mo><mi>&alpha;</mi></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("ptnm1",
{
  q <- mathml(quote(pt(t, N - 1)))
  expect_equal(q, "<math><mrow><mi>P</mi><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>T</mi><mo>&le;</mo><mi>t</mi></mrow><mo>;</mo><mrow><mrow><mi>N</mi><mo>-</mo><mn>1.00</mn></mrow><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("qtalpha2nm1",
{
  q <- mathml(quote(qt(alpha/2, N - 1)))
  expect_equal(q, "<math><mrow><msub><mi>T</mi><mrow><mi>&alpha;</mi><mo>/</mo><mn>2.00</mn></mrow></msub><mo>&af;</mo><mrow><mo>(</mo><mrow><mrow><mi>N</mi><mo>-</mo><mn>1.00</mn></mrow><mspace width=\"thinmathspace\"></mspace><mtext>df</mtext></mrow><mo>)</mo></mrow></mrow></math>")
})

test_that("AXB",
{
  q <- mathml(quote(A %*% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&times;</mo><mi>B</mi></mrow></math>")
})

test_that("AdotB",
{
  q <- mathml(quote(A %.% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&sdot;</mo><mi>B</mi></mrow></math>")
})

test_that("AxB",
{
  q <- mathml(quote(A %x% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&CircleTimes;</mo><mi>B</mi></mrow></math>")
})

test_that("A/B",
{
  q <- mathml(quote(A %/% B))
  expect_equal(q, "<math><mrow><mo>&lfloor;</mo><mrow><mi>A</mi><mo>/</mo><mi>B</mi></mrow><mo>&rfloor;</mo></mrow></math>")
})

test_that("A mod B",
{
  q <- mathml(quote(A %% B))
  expect_equal(q, "<math><mrow><mo>&lceil;</mo><mrow><mi>A</mi><mo>/</mo><mi>B</mi></mrow><mo>&rceil;</mo></mrow></math>")
})

test_that("A&B",
{
  q <- mathml(quote(A & B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&and;</mo><mi>B</mi></mrow></math>")
})

test_that("AvelB",
{
  q <- mathml(quote(A | B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&or;</mo><mi>B</mi></mrow></math>")
})

test_that("AeqeqB",
{
  q <- mathml(quote(A == B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>=</mo><mi>B</mi></mrow></math>")
})

test_that("AarroweqB",
{
  q <- mathml(quote(A <- B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>=</mo><mi>B</mi></mrow></math>")
})

test_that("AdiffB",
{
  q <- mathml(quote(A != B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&ne;</mo><mi>B</mi></mrow></math>")
})

test_that("AtildetildeB",
{
  q <- mathml(quote(A %~~% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&approx;</mo><mi>B</mi></mrow></math>")
})

test_that("AequivB",
{
  q <- mathml(quote(A %==% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&equiv;</mo><mi>B</mi></mrow></math>")
})

test_that("AcongB",
{
  q <- mathml(quote(A %=~% B))
  expect_equal(q,  "<math><mrow><mi>A</mi><mo>&cong;</mo><mi>B</mi></mrow></math>")
})

test_that("ApropB",
{
  q <- mathml(quote(A %prop% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&prop;</mo><mi>B</mi></mrow></math>")
})

test_that("AisinB",
{
  q <- mathml(quote(A %in% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&isin;</mo><mi>B</mi></mrow></math>")
})

test_that("AlrarrowB",
{
  q <- mathml(quote(A %<->% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&leftrightarrow;</mo><mi>B</mi></mrow></math>")
})

test_that("ArarrowB",
{
  q <- mathml(quote(A %->% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&rightarrow;</mo><mi>B</mi></mrow></math>")
})

test_that("AlarrowB",
{
  q <- mathml(quote(A %<-% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&leftarrow;</mo><mi>B</mi></mrow></math>")
})

test_that("AuparrowB",
{
  q <- mathml(quote(A %up% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&leftarrow;</mo><mi>B</mi></mrow></math>")
})

test_that("AdownarrowB",
{
  q <- mathml(quote(A %down% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&downarrow;</mo><mi>B</mi></mrow></math>")
})

test_that("AlrdoublearrowB",
{
  q <- mathml(quote(A %<=>% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&iff;</mo><mi>B</mi></mrow></math>")
})

test_that("AldoublearrowB",
{
  q <- mathml(quote(A %<=% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&lArr;</mo><mi>B</mi></mrow></math>")
})

test_that("ArdoublearrowB",
{
  q <- mathml(quote(A %=>% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&rArr;</mo><mi>B</mi></mrow></math>")
})

test_that("AupdoublearrowB",
{
  q <- mathml(quote(A %dblup%  B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&uArr;</mo><mi>B</mi></mrow></math>")
})

test_that("AdowndoublearrowB",
{
  q <- mathml(quote(A %dbldown% B))
  expect_equal(q, "<math><mrow><mi>A</mi><mo>&dArr;</mo><mi>B</mi></mrow></math>")
})


