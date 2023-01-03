:- discontiguous math/2, math/3, math/4, current/3, paren/3, prec/3, type/3, denoting/3, ml/3, jax/3, mathml/0.
:- use_module(library(http/html_write)).

% Here you can define your own macros
%
% Example: assert(math_hook(t0, subscript(t, 0)))
%
:- dynamic math_hook/2.

%
% R interface: Translate R expression to MathML string
%
r2mathml(R, S)
=> r2mathml([], R, S).

r2mathml(Flags, R, S)
 => mathml(Flags, R, M),
    html(M, H, []),
    maplist(atom_string, H, S).

% Same for MathJax/LaTeX
r2mathjax(R, S)
 => r2mathjax([], R, S).

r2mathjax(Flags, R, S)
 => mathjax(Flags, R, S).

%
% Translate R expression to HTML term
%
mathml(Flags, R, M)
 => ml(Flags, R, M0),
    denoting(Flags, R, Denoting),
    ml(Flags, with(Denoting), With),
    !, M = [math(M0), With].

mathjax(Flags, R, M)
 => jax(Flags, R, M0),
    denoting(Flags, R, Denoting),
    jax(Flags, with(Denoting), With),
    !, format(string(M), "$~w$~w", [M0, With]).

%
% Macros
%
% macro(Flags, R, New, M): translates the R expression to another R
% expression M, checking for Flags and eventually changing Flags to New
%
% Calls math/2,3,4 macros
%
macro(Flags, R, New, M) :-
    math_hook(R, M0),
    !, New = Flags,
    M = M0.

macro(Flags, R, New, M) :-
    math(Flags, R, New, M),
    dif(Flags-R, New-M).

macro(Flags, R, Flags, M) :-
    math(Flags, R, M),
    dif(R, M).

macro(Flags, R, Flags, M) :-
    math(R, M),
    dif(R, M).

%
% Main MathML translation
%
% Flags: to control some aspects of the output
% R: R expression
% M: HTML term
%
% This predicate only checks if a macro can be applied. Add ml/3 predicates for
% R expressions with their translation below.
%
ml(R, M, Flags),
    macro(R, New, R0, Flags)
 => ml(New, R0, M).

% Same for MathJax/LaTeX
jax(Flags, R, M),
    macro(Flags, R, New, R0)
 => jax(New, R0, M).

% Same for precedence checks, parentheses and types
prec(Flags, R, Prec),
    macro(Flags, R, New, R0)
 => prec(New, R0, Prec).

paren(Flags, R, Paren),
    macro(Flags, R, New, R0)
 => paren(New, R0, Paren).

type(Flags, R, Paren),
    macro(Flags, R, New, R0)
 => type(New, R0, Paren).

%
% Check if the current element is to be replaced
%
math(Flags, R, New, M),
    member(replace(R, _), Flags)
 => select(replace(R, M), Flags, New).

%
% Examples
%
mathml(R) :-
    r2mathml([], R, M),
    atomic_list_concat(M, S),
    writeln(R-S).

mathjax(R) :-
    r2mathjax([], R, M),
    atomic_list_concat(M, S),
    writeln(R-S).

%
% Content starts here
%

% Summation sign from to (subscript and superscript)
math(Flags, Sum, New, M),
    compound(Sum),
    compound_name_arguments(Sum, sum, [Arg]),
    select(subscript(From), Flags, New0),
    select(superscript(To), New0, New1)
 => New = New1,
    M = fn(subsupscript(sum, From, To), [Arg]).

% Summation sign from (only subscript)
math(Flags, Sum, New, M),
    compound(Sum),
    compound_name_arguments(Sum, sum, [Arg]),
    select(subscript(Idx), Flags, New0)
 => New = New0,
    M = fn(subscript(sum, Idx), [Arg]).

mathml :-
    mathml(sum('['(x, i))).

% Same for product sign
math(Flags, Prod, New, M),
    compound(Prod),
    compound_name_arguments(Prod, prod, Arg),
    select(subscript(Idx), Flags, New0),
    select(superscript(Pwr), New0, New1)
 => New = New1,
    M = fn(subsupscript(prod, Idx, Pwr), Arg).

math(Flags, Prod, New, M),
    compound(Prod),
    compound_name_arguments(Prod, prod, Arg),
    select(subscript(Idx), Flags, New0)
 => New = New0,
    M = fn(subscript(prod, Idx), Arg).

mathml :-
    mathml(prod('['(x, i))).

%
% Subscript and superscript
%
% Moves sub- and superscript into the flags, so that the next available
% function can handle them. This is needed, e.g., if colors are to be
% skipped.
%
math(Flags, subsupscript(R, Idx, Pwr), New, M)
 => New = [subscript(Idx), superscript(Pwr) | Flags],
    M = R.

% Render
ml(R, M, Flags),
    select(subscript(Idx), Flags, New0),
    select(superscript(Pwr), New0, New)
 => ml(New, R, X),
    ml(New, Idx, Y),
    ml(New, Pwr, Z),
    M = msubsup([X, Y, Z]).

jax(Flags, R, M),
    select(subscript(Idx), Flags, New0),
    select(superscript(Pwr), New0, New)
 => jax(New, R, X),
    jax(New, Idx, Y),
    jax(New, Pwr, Z),
    format(string(M), "{~w}_{~w}^{~w}", [X, Y, Z]).

type(Flags, R, Type),
    member(subscript(Idx), Flags),
    member(superscript(Pwr), Flags)
 => Type = subsupscript(R, Idx, Pwr).

mathml :-
    mathml(subsupscript(x, i, 2)).

%
% Subscript like s_D
%
math(Sub, M),
    compound(Sub),
    compound_name_arguments(Sub, '[', [R | Idx])
 => M = subscript(R, list("", Idx)).

math(Flags, subscript(R, Idx), New, M)
 => New = [subscript(Idx) | Flags],
    M = R.

% Render
ml(R, M, Flags),
    select(subscript(Idx), Flags, New)
 => ml(New, R, X),
    ml(New, Idx, Y),
    M = msub([X, Y]).

jax(Flags, R, M),
    select(subscript(Idx), Flags, New)
 => jax(New, R, X),
    jax(New, Idx, Y),
    format(string(M), "{~w}_{~w}", [X, Y]).

type(Flags, R, Type),
    member(subscript(Idx), Flags)
 => Type = subscript(R, Idx).

mathml :-
    mathml(subscript(x, i)).

%
% Superscripts like s^D
%
math(Flags, superscript(R, Idx), New, M)
 => New = [superscript(Idx) | Flags],
    M = R.

% Render
ml(R, M, Flags),
    select(superscript(Pwr), Flags, New),
    prec(Flags, R, Prec),
    current_op(P, xfy, ^),
    Prec >= P
 => ml(New, paren(R), X),
    ml(New, Pwr, Y),
    M = msup([X, Y]).

ml(R, M, Flags),
    select(superscript(Pwr), Flags, New)
 => ml(New, R, X),
    ml(New, Pwr, Y),
    M = msup([X, Y]).

jax(Flags, R, M),
    prec(Flags, R, Prec),
    current_op(P, xfy, ^),
    Prec >= P,
    select(superscript(Pwr), Flags, New)
 => jax(New, paren(R), X),
    jax(New, Pwr, Y),
    format(string(M), "{~w}^{~w}", [X, Y]).

jax(Flags, R, M),
    select(superscript(Pwr), Flags, New)
 => jax(New, R, X),
    jax(New, Pwr, Y),
    format(string(M), "{~w}^{~w}", [X, Y]).

type(Flags, R, Type),
    member(superscript(Pwr), Flags)
 => Type = superscript(R, Pwr).

mathml :-
    mathml(superscript(x, 2)).

%
% Suppress the names of function arguments from R
%
math(name(_)=R, M)
 => M = R.

%
% Upright text
%
math(R, M),
    string(R)
 => M = text(R).

ml(text(R), M, _Flags)
 => M = mtext(R).

jax(_Flags, text(R), M)
 => format(string(M), "\\mathrm{~w}", [R]).

type(_Flags, text(_), T)
 => T = atomic.

mathml :-
    mathml("text").

mathjax :-
    mathjax("text").

%
% Greek letters
%
math(R, M),
    atom(R),
    memberchk(R, [alpha, beta, gamma, delta, epsilon, varepsilon, zeta, eta,
        theta, vartheta, iota, kappa, lambda, mu, nu, xi, pi, rho, sigma,
        varsigma, tau, upsilon, phi, varphi, chi, psi, omega, 'Gamma', 'Delta',
        'Theta', 'Lambda', 'Xi', 'Pi', 'Sigma', 'Upsilon', 'Phi', 'Psi',
        'Omega'])
 => M = greek(R).

ml(greek(R), M, _Flags)
 => M = mi(&(R)).

jax(_Flags, greek(R), M)
 => format(string(M), "\\~w", [R]).

type(_Flags, greek(_), T)
 => T = atomic.

mathml :-
    mathml(alpha).

%
% Booleans
%
math('TRUE', M)
 => M = boolean('T').

math('FALSE', M)
 => M = boolean('F').

ml(boolean(R), M, _Flags)
 => M = mi(R).

jax(_Flags, boolean(R), M)
 => format(string(M), "~w", [R]).

type(_Flags, boolean(_), T)
 => T = atomic.

mathml :-
    mathml('TRUE'),
    mathml('FALSE').

%
% Set
%
math('is.null'(R), M)
 => M = (R == 'NULL').

math('NULL', M)
 => M = set(empty).

ml(set(empty), M, _Flags)
 => M = mi(&(empty)).

jax(_Flags, set(empty), M)
 => M = "\\emptyset".

type(_Flags, set(empty), T)
 => T = atomic.

%
% sin^2(x) etc.
%
math(Flags, sin(A), New, M),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(sin^Pwr, [A]).

math(Flags, sinpi(A), New, M),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(sinpi^Pwr, [A]).

math(Flags, cos(A), New, M),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(cos^Pwr, [A]).

math(Flags, cospi(A), New, M),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(cospi^Pwr, [A]).

math(Flags, tan(A), New, M),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(tan^Pwr, [A]).

math(Flags, tanpi(A), New, M),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(tanpi^Pwr, [A]).

%
% Special functions
%
special(R) :-
    atom(R),
    memberchk(R, [sgn, sin, cos, tan, asin, acos, atan, arctan2, sinh, cosh, tanh,
               arsinh, arcosh, artanh, log, exp, sum, prod, min, max]).

math(R, M),
    special(R)
 => M = special(R).

ml(special(sum), M, _Flags)
 => M = mo(&(sum)).

prec(_Flags, special(sum), Prec)
 => current(P, yfx, *),
    Prec is P + 1.

ml(special(prod), M, _Flags)
 => M = mo(&(prod)).

prec(_Flags, special(prod), Prec)
 => current(P, yfx, *),
    Prec is P + 1.

ml(special(R), M, _Flags)
 => M = mi(R).

jax(_Flags, special(sum), M)
 => M = "{\\sum}".

jax(_Flags, special(prod), M)
 => M = "{\\prod}".

jax(_Flags, special(sgn), M)
 => M = "{\\mathrm{sgn}\\,}".

jax(_Flags, special(R), M)
 => format(string(M), "\\~w", [R]).

type(_Flags, special(_), T)
 => T = special.

prec(_Flags, special(_), Prec)
 => Prec = 0. % current(Prec, yfx, *).

pmathml :-
    mathml(exp(x)),
    mathml(exp(x + y)).

%
% Space
%
math(space, M)
 => M = space(thinmathspace).

ml(space(W), M, _Flags)
 => M = mspace(width(W), []).

jax(_Flags, space(thinmathspace), M)
 => M = "\\,".

jax(_Flags, space(_Width), M)
 => M = "\\ ".

%
% Symbols/Identifiers
%
math(R, M),
    atom(R)
 => M = ident(R).

ml(ident(R), M, Flags),
    member(mathvariant(calligraphy), Flags)
 => M = mi(mathvariant(script), R).

ml(ident(R), M, _Flags)
 => M = mi(R).

jax(Flags, ident(R), M),
    member(mathvariant(calligraphy), Flags)
 => format(string(M), "\\mathcal{~w}", [R]).

jax(_Flags, ident(R), M)
 => format(string(M), "~w", [R]).

type(_Flags, ident(_), T)
 => T = atomic.

%
% Linear model
%
math(lm(F, _Data), M)
 => M = F.

%
% R package base
%
math(return(X), M)
 => M = X.

math(length(R), M)
 => M = abs(R).

ml(abs(R), M, Flags)
 => ml(R, X, Flags),
    M = mrow([mo(&(vert)), X, mo(&(vert))]).

jax(Flags, abs(R), M)
 => jax(Flags, R, X),
    format(string(M), "\\left|~w\\right|", [X]).

paren(_Flags, abs(_), P)
 => P = 0.

prec(Flags, abs(R), P)
 => prec(Flags, paren(R), P).

math(sign(R), M)
 => M = fn(sgn, [R]).

ml(sqrt(R), M, Flags)
 => ml(R, X, Flags),
    M = msqrt(X).

jax(Flags, sqrt(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\sqrt{~w}", [X]).

paren(_Flags, sqrt(_), P)
 => P = 0.

prec(_Flags, sqrt(_), P)
 => current_op(P0, xfy, ^),
    P is P0 + 1.

math(sin(A), M)
 => M = fn(sin, [A]).

math(cos(A), M)
 => M = fn(cos, [A]).

math(tan(A), M)
 => M = fn(tan, [A]).

math(asin(A), M)
 => M = fn(superscript(sin, -1), [A]).

math(acos(A), M)
 => M = fn(superscript(cos, -1), [A]).

math(atan(A), M)
 => M = fn(superscript(tan, -1), [A]).

math(atan2(A, B), M)
 => M = fn(superscript(tan, -1), [A, B]).

math(sinpi(A), M)
 => M = fn(sin, [A*pi]).

math(cospi(A), M)
 => M = fn(cos, [A*pi]).

math(tanpi(A), M)
 => M = fn(tan, [A*pi]).

math(sinh(A), M)
 => M = fn(sinh, [A]).

math(cosh(A), M)
 => M = fn(cosh, [A]).

math(tanh(A), M)
 => M = fn(tanh, [A]).

math(asinh(A), M)
 => M = fn(superscript(sinh, -1), [A]).

math(acosh(A), M)
 => M = fn(superscript(cosh, -1), [A]).

math(atanh(A), M)
 => M = fn(superscript(tanh, -1), [A]).

math(all(A), M)
 => M = forall(A).

ml(forall(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo(&('ForAll')), mo(&(af)), X]).

jax(Flags, forall(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\forall{~w}", [X]).

paren(Flags, forall(A), P)
 => paren(Flags, A, P).

prec(_Flags, forall(_), P)
 => current(P, yfx, *).

math(any(A), M)
 => M = exists(A).

ml(exists(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo(&('Exists')), mo(&(af)), X]).

jax(Flags, exists(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\exists{~w}", [X]).

paren(Flags, exists(A), P)
 => paren(Flags, A, P).

prec(_Flags, exists(_), P)
 => current(P, yfx, *).

% todo: expon.scaled
math(besselI(x=X, nu=Nu), M)
 => M = besselI(X, Nu).

math(besselI(X, Nu), M)
 => M = fn(subscript('I', Nu), [paren(X)]).

% todo: expon.scaled
math(besselK(x=X, nu=Nu), M)
 => M = besselK(X, Nu).

math(besselK(X, Nu), M)
 => M = fn(subscript('K', Nu), [paren(X)]).

math(besselJ(x=X, nu=Nu), M)
 => M = besselJ(X, Nu).

math(besselJ(X, Nu), M)
 => M = fn(subscript('J', Nu), [paren(X)]).

math(besselY(x=X, nu=Nu), M)
 => M = besselY(X, Nu).

math(besselY(X, Nu), M)
 => M = fn(subscript('Y', Nu), [paren(X)]).

math(beta(a=A, b=B), M)
 => M = beta(A, B).

math(beta(A, B), M)
 => M = fn('B', [A, B]).

math(lbeta(a=A, b=B), M)
 => M = lbeta(A, B).

math(lbeta(A, B), M)
 => M = log(beta(A, B)).

math(gamma(A), M)
 => M = fn('Gamma', [paren(A)]).

math(lgamma(A), M)
 => M = log(gamma(A)).

math(digamma(A), M)
 => M = frac(d, d*A) * log(gamma(A)).

math(trigamma(A), M)
 => M = frac(d^2, (d*A)^2) * log(gamma(A)).

math(psigamma(x=A, deriv=Deriv), M)
 => M = psigamma(A, Deriv).

math(psigamma(A, Deriv), M)
 => M = frac(d^(Deriv+2), (d*A)^(Deriv+2)) * log(gamma(A)).

math(choose(n=N, k=K), M)
 => M = choose(N, K).

ml(choose(N, K), M, Flags)
 => ml(N, X, Flags),
    ml(K, Y, Flags),
    M = mrow([mo('('), mfrac([linethickness(0)], [X, Y]), mo(')')]).

jax(Flags, choose(N, K), M)
 => jax(Flags, N, X),
    jax(Flags, K, Y),
    format(string(M), "\\binom{~w}{~w}", [X, Y]).

paren(_Flags, choose(_, _), P)
 => P = 1.

prec(_Flags, choose(_, _), P)
 => P = 0.

type(_Flags, choose(_, _), T)
 => T = paren.

math(lchoose(N, K), M)
 => M = log(choose(N, K)).

math(factorial(N), M)
 => current(Prec, xfy, ^),
    M = yf(Prec, !, N).

math(lfactorial(N), M)
 => M = log(factorial(N)).

math(and(A, B), M)
 => current(Prec, xfy, ','),
    M = xfy(Prec, and, A, B).

math(or(A, B), M)
 => current(Prec, xfy, ';'),
    M = xfy(Prec, or, A, B).

math(!(A), M)
 => current(Prec, xfy, ^),
    M = fy(Prec, not, A).

math(xor(x=A, y=B), M)
 => M = xor(A, B).

math(xor(A, B), M)
 => current(Prec, xfy, ';'),
    M = xfy(Prec, veebar, A, B).

math(exp(A), M)
 => M = fn(exp, [A]).

math(expm1(A), M)
 => M = exp(A) - 1.

math(log(X), M)
 => M = fn(log, [X]).

math(log10(X), M)
 => M = logb(X, 10).

math(log2(X), M)
 => M = logb(X, 2).

% unclear why the cases need to be handled separately
math(logb(X, B), M)
 => M = fn(subscript(log, B), [X]).

math(log1p(A), M)
 => M = 1 + log(A).

ml(ceiling(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo(&(lceil)), X, mo(&(rceil))]).

jax(Flags, ceiling(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\lceil{~w}\\rceil", [X]).

paren(_Flags, ceiling(_), P)
 => P is 0.

ml(floor(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo(&(lfloor)), X, mo(&(rfloor))]).

jax(Flags, floor(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\lfloor{~w}\\rfloor", [X]).

paren(_Flags, floor(_), P)
 => P is 0.

math(_Flags, (_F :- Body), M)
 => M = Body.

math(_Flags, Function, M),
    compound(Function),
    compound_name_arguments(Function, function, Args)
 => M = lambda(Args).

math(Flags, lambda(Args), M),
    member(name-N, Flags)
 => M = fn(N, Args).

math(Flags, lambda(Args), M)
 => option(name(N), Flags),
    M = fn(N, Args).

math(_Flags, lambda(Args), M)
 => M = fn(lambda, Args).

type(_Flags, lambda(_), T)
 => T = special.

math('<-'(R, S), M)
 => M = (R == S).

math(function(na, Body, _), M),
    compound(Body),
    compound_name_arguments(Body, '{', Args)
 => M = body(Args).

math(function(na, Body, _), M)
 => M = body([[Body]]).

math(Curly, M),
    compound(Curly),
    compound_name_arguments(Curly, '{', Args)
 => M = body(Args).

ml(body([R]), M, Flags)
 => ml(R, M, Flags).

ml(body(Body), M, Flags)
 => maplist(ml(Flags), Body, R),
    M = mrow([mo('{'), mtable(columnalign(left), R)]).

jax(Flags, body([R]), M)
 => jax(Flags, R, M).

jax(Flags, body(Body), M)
 => maplist(jax(Flags), Body, Rows),
    format(string(M), "\\left\\{\\begin{array}{l}~w\\end{array}\\right.", [Rows]).

% Vectors
math(Hash, M),
    compound(Hash),
    compound_name_arguments(Hash, Name, Elements),
    member(Name, ['#', '$$', '%', '!'])
 => M = paren(Elements).

% Matrices
ml(Matrix, M, Flags),
    compound(Matrix),
    compound_name_arguments(Matrix, Name, Rows),
    member(Name, ['##', '$$$', '%%', '!!'])
 => maplist(ml_row(Flags), Rows, R),
    M = mrow([mo('('), mtable(columnalign(left), R), mo(')')]).

ml_row(Flags, Row, M),
    compound(Row),
    compound_name_arguments(Row, Name, Cells),
    member(Name, ['#', '$$', '%', '!'])
 => maplist(ml_cell(Flags), Cells, C),
    M = mtr(C).

ml_cell(Cell, M, Flags)
 => ml(Cell, C, Flags),
    M = mtd(C).

jax(Flags, Matrix, M),
    compound(Matrix),
    compound_name_arguments(Matrix, Name, [Row1 | Rows]),
    member(Name, ['##', '$$$', '%%', '!!'])
 => findall(c, arg(_, Row1, _), Ls),
    atomic_list_concat(Ls, LLL),
    maplist(jax_row(Flags), [Row1 | Rows], R),
    atomic_list_concat(R, Lines),
    format(string(M), "\\left(\\begin{array}{~w}~w\\end{array}\\right)", [LLL, Lines]).

jax_row(Flags, Row, M),
    compound(Row),
    compound_name_arguments(Row, Name, Cells),
    member(Name, ['#', '$$', '%', '!'])
 => maplist(jax_cell(Flags), Cells, C),
    atomic_list_concat(C, ' & ', R),
    format(string(M), "~w\\\\\n", [R]).

jax_cell(Flags, C, M)
 => jax(Flags, C, X),
    format(string(M), "~w", [X]).

math(Identical, M),
    compound(Identical),
    compound_name_arguments(Identical, identical, [X, Y])
 => M = (X == Y).

ml(ifelse(T, Y, N), M, Flags)
 => ml(T, Test, Flags),
    ml(Y, Yes, Flags),
    ml(N, No, Flags),
    ml(space, S, Flags),
    M = mrow([mo('{'),
      mtable(columnalign(left),
      [ mtr([Yes, mrow([mtext("if"), S, Test])]),
        mtr([No, mtext("otherwise")])
      ])]).

jax(Flags, ifelse(T, Y, N), M)
 => jax(Flags, T, Test),
    jax(Flags, Y, Yes),
    jax(Flags, N, No),
    format(string(M),
      "\\left\\{\\begin{array}{ll} {~w} & \\mathrm{if}~~{~w}\\\\ {~w} & \\mathrm{otherwise}\\end{array}\\right.",
      [Yes, Test, No]).

paren(_Flags, ifelse(_, _, _), P)
 => P is 0.

ml(if(T, Y), M, Flags)
 => ml(T, Test, Flags),
    ml(Y, Yes, Flags),
    ml(space, S, Flags),
    M = mrow([Yes, mtext(","), S, mtext("if"), S, Test]).

jax(Flags, if(T, Y), M)
 => jax(Flags, T, Test),
    jax(Flags, Y, Yes),
    format(string(M), "{~w},\\ \\mathrm{if}\\ {~w}", [Yes, Test]).

paren(_Flags, if(_, _), P)
 => P is 0.

math('%in%'(X, Y), M)
 => M = isin(X, Y).

math(setdiff(X, Y), M)
 => M = X - Y.

math('%x%'(X, Y), M)
 => M = kronecker(X, Y).

math('&'(A, B), M)
 => M = and(A, B).

math('|'(A, B), M)
 => M = or(A, B).

ml(Prod, M, Flags),
    compound(Prod),
    compound_name_arguments(Prod, prod, Args)
 => maplist(ml(Flags), Args, MX),
    M = mrow([mo(&(prod)), mrow(MX)]).

jax(Flags, prod(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\prod{~w}", [X]).

jax(Flags, Prod, M),
    compound(Prod),
    compound_name_arguments(Prod, prod, Args)
 => maplist(jax(Flags), Args, X),
    format(string(M), "\\prod{~w}", [X]).

paren(Flags, Prod, P),
    compound(Prod),
    compound_name_arguments(Prod, prod, Args)
 => maplist(paren(Flags), Args, PX),
    max_list(PX, P).

prec(_Flags, Prod, P),
    compound(Prod),
    compound_name_arity(Prod, prod, _)
 => current(P, yfx, *).

math(mean(A), M)
 => M = overline(A).

math(Min, M),
    compound(Min),
    compound_name_arguments(Min, min, Args)
 => M = fn(min, Args).

math(Max, M),
    compound(Max),
    compound_name_arguments(Max, max, Args)
 => M = fn(max, Args).

math(t(A), M)
 => M = A^"T".

math(Which, M),
    compound(Which),
    compound_name_arguments(Which, which, Args)
 => M = subscript("I", Args).

math('which.max'(A), M)
 => M = fn("argmax", [A]).

math('which.min'(A), M)
 => M = fn("argmin", [A]).

% calligraphic letters
math(Flags, cal(A), New, M)
 => New = [mathvariant(calligraphy) | Flags],
    M = A.

%
% Extract value from a result (e.g., integrate)
%
math(Flags, $(Fn, "value"), New, M)
 => Flags = New,
    M = Fn.

%
% Integrate over range
%
% Case A: Fn is a function
math(Flags, integrate(Fn, Lower, Upper), New, M),
    Fn = (Head :- _Body),
    compound(Head),
    compound_name_arguments(Head, function, [DX | _]),
    member(name-Name, Flags)
 => Flags = New,
    M = integrate(fn(Name, [DX]), Lower, Upper, DX).

math(Flags, integrate(Fn, Lower, Upper), New, M),
    Fn = (Head :- _Body),
    compound(Head),
    compound_name_arguments(Head, function, [DX | _])
 => Flags = New,
    M = integrate(fn(lambda, [DX]), Lower, Upper, DX).

% Case B: Fn is an atom
math(Flags, integrate(Fn, Lower, Upper), New, M),
    atom(Fn)
 => Flags = New,
    r_eval('['(formalArgs(args(Fn)), 1), Arg1),
    atom_string(DX, Arg1),
    M = integrate(fn(Fn, [DX]), Lower, Upper, DX).

% Internal
ml(integrate(Fn, From, To, DX), M, Flags)
 => ml(Fn, XFn, Flags),
    ml(From, XFrom, Flags),
    ml(To, XTo, Flags),
    ml(DX, XDX, Flags),
    ml(space, Space, Flags),
    M = mrow([munderover([mo(&(int)), XFrom, XTo]), XFn, Space, mi(d), XDX]).

jax(Flags, integrate(Fn, From, To, DX), M)
 => jax(Flags, Fn, XFn),
    jax(Flags, From, XFrom),
    jax(Flags, To, XTo),
    jax(Flags, DX, XDX),
    format(string(M), "\\int_{~w}^{~w}{~w}\\,{d{~w}}", [XFrom, XTo, XFn, XDX]).

paren(Flags, integrate(_, _, _, A), Paren)
 => paren(Flags, A, Paren).

prec(_Flags, integrate(_, _, _, _), Prec)
 => current(Prec, yfx, *).

% hats
ml(hat(A), M, Flags)
 => ml(A, X, Flags),
    M = mover(accent(true), [X, mo(&('Hat'))]).

jax(Flags, hat(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\hat{~w}", [X]).

paren(Flags, hat(A), Paren)
 => paren(Flags, A, Paren).

prec(Flags, hat(A), Prec)
 => prec(Flags, A, Prec).

type(Flags, hat(A), Type)
 => type(Flags, A, Type).

% tilde
ml(tilde(A), M, Flags)
 => ml(A, X, Flags),
	M = mover(accent(true), [X, mo(&(tilde))]).

paren(Flags, tilde(A), Paren)
 => paren(Flags, A, Paren).

prec(Flags, tilde(A), Prec)
 => prec(Flags, A, Prec).

type(Flags, tilde(A), Type)
 => type(Flags, A, Type).

%
% Mathematical signs
%
ml(op(le), M, _Flags)
 => M = mo(&(le)).

jax(_Flags, op(le), M)
 => M = "\\le".

ml(op(ge), M, _Flags)
 => M = mo(&(ge)).

jax(_Flags, op(ge), M)
 => M = "\\ge".

ml(op(ne), M, _Flags)
 => M = mo(&(ne)).

jax(_Flags, op(ne), M)
 => M = "\\ne".

ml(op(cdot), M, _Flags)
 => M = mo(&(sdot)).

jax(_Flags, op(cdot), M)
 => M = "\\cdot".

ml(op(pm), M, _Flags)
 => M = mo(&(pm)).

jax(_Flags, op(pm), M)
 => M = "\\pm".

ml(op(times), M, _Flags)
 => M = mo(&(times)).

jax(_Flags, op(times), M)
 => M = "\\times".

ml(op(sum), M, _Flags)
 => M = mo(&(sum)).

jax(_Flags, op(sum), M)
 => M = "\\sum".

ml(op(prod), M, _Flags)
 => M = mo(&(prod)).

jax(_Flags, op(prod), M)
 => M = "\\prod".

ml(op('#58'), M, _Flags)
 => M = mo(&('#58')).

jax(_Flags, op('#58'), M)
 => M = ":".

ml(op(','), M, _Flags)
 => M = mo(',').

jax(_Flags, op(','), M)
 => M = ",".

ml(op('CircleTimes'), M, _Flags)
 => M = mo(&('CircleTimes')).

jax(_Flags, op('CircleTimes'), M)
 => M = "\\otimes".

ml(op('#x2062'), M, _Flags)
 => M = mo(&('#x2062')).

jax(_Flags, op('#x2062'), M)
 => M = "{}".

ml(op('Tilde'), M, _Flags)
 => M = mo(&('Tilde')).

jax(_Flags, op('Tilde'), M)
 => M = "\\sim".

ml(op(leftrightarrow), M, _Flags)
 => M = mo(&(leftrightarrow)).

jax(_Flags, op(leftrightarrow), M)
 => M = "\\leftrightarrow".

ml(op(iff), M, _Flags)
 => M = mo(&(iff)).

jax(_Flags, op(iff), M)
 => M = "\\iff".

ml(op(rightarrow), M, _Flags)
 => M = mo(&(rightarrow)).

jax(_Flags, op(rightarrow), M)
 => M = "\\rightarrow".

ml(op(rArr), M, _Flags)
 => M = mo(&(rArr)).

jax(_Flags, op(rArr), M)
 => M = "\\Rightarrow".

ml(op(leftarrow), M, _Flags)
 => M = mo(&(leftarrow)).

jax(_Flags, op(leftarrow), M)
 => M = "\\leftarrow".

ml(op(lArr), M, _Flags)
 => M = mo(&(lArr)).

jax(_Flags, op(lArr), M)
 => M = "\\Leftarrow".

ml(op(uparrow), M, _Flags)
 => M = mo(&(uparrow)).

jax(_Flags, op(uparrow), M)
 => M = "\\uparrow".

ml(op(uArr), M, _Flags)
 => M = mo(&(uArr)).

jax(_Flags, op(uArr), M)
 => M = "\\Uparrow".

ml(op(downarrow), M, _Flags)
 => M = mo(&(downarrow)).

jax(_Flags, op(downarrow), M)
 => M = "\\downarrow".

ml(op(dArr), M, _Flags)
 => M = mo(&(dArr)).

jax(_Flags, op(dArr), M)
 => M = "\\Downarrow".

ml(op(approx), M, _Flags)
 => M = mo(&(approx)).

jax(_Flags, op(approx), M)
 => M = "\\approx".

ml(op(equiv), M, _Flags)
 => M = mo(&(equiv)).

jax(_Flags, op(equiv), M)
 => M = "\\equiv".

ml(op(cong), M, _Flags)
 => M = mo(&(cong)).

jax(_Flags, op(cong), M)
 => M = "\\cong".

ml(op(propto), M, _Flags)
 => M = mo(&(prop)).

jax(_Flags, op(propto), M)
 => M = "\\propto".

ml(op(and), M, _Flags)
 => M = mo(&(and)).

jax(_Flags, op(and), M)
 => M = "\\land".

ml(op(or), M, _Flags)
 => M = mo(&(or)).

jax(_Flags, op(or), M)
 => M = "\\lor".

ml(op(not), M, _Flags)
 => M = mo(&(not)).

jax(_Flags, op(not), M)
 => M = "\\lnot".

ml(op(veebar), M, _Flags)
 => M = mo(&(veebar)).

jax(_Flags, op(veebar), M)
 => M = "\\veebar".

ml(op(isin), M, _Flags)
 => M = mo(&(isin)).

jax(_Flags, op(isin), M)
 => M = "\\in".

ml(op(notin), M, _Flags)
 => M = mo(&(notin)).

jax(_Flags, op(notin), M)
 => M = "\\notin".

ml(op(cap), M, _Flags)
 => M = mo(&(cap)).

jax(_Flags, op(cap), M)
 => M = "\\cap".

ml(op(cup), M, _Flags)
 => M = mo(&(cup)).

jax(_Flags, op(cup), M)
 => M = "\\cup".

ml(op(A), M, _Flags)
 => M = mo(A).

jax(_Flags, op(A), M)
 => format(string(M), "~w", [A]).

prec(_Flags, op(A), P),
    current(P0, _Fix, A)
 => P = P0.

current(0, fy, op(sum)).

denoting(_Flags, op(_), D)
 => D = [].


%
% Numbers
%
math(A, M),
    integer(A),
    A >= 0
 => M = posint(A).

math(A, M),
    integer(A)
 => M = integer(A).

math(integer(A), M),
    A >= 0
 => M = posint(A).

math(integer(A), M)
 => Abs is abs(A),
    M = -posint(Abs).

math(A, M),
    number(A),
    A >= 0
 => M = pos(A).

math(A, M),
    number(A)
 => M = number(A).

ml(posint(A), M, _Flags)
 => M = mn(A).

ml(pos(1.0Inf), M, _Flags)
 => M = mi(&('#x221E')).

ml(pos(A), M, Flags)
 => option(round(D), Flags, 2),
    format(atom(Mask), '~~~wf', [D]),
    format(string(X), Mask, [A]),
    M = mn(X).

jax(_Flags, posint(A), M)
 => format(string(M), "~w", [A]).

jax(_Flags, pos(1.0Inf), M)
 => M = "\\infty".

jax(Flags, pos(A), M)
 => option(round(D), Flags, 2),
    format(atom(Mask), '~~~wf', [D]),
    format(string(M), Mask, [A]).

type(_Flags, pos(_), Type)
 => Type = atomic.

type(_Flags, posint(_), Type)
 => Type = atomic.

math(number(A), M),
    A < 0
 => Abs is abs(A),
    M = -pos(Abs).

math(number(A), M)
 => M = pos(A).

%
% Hypothenuse
%
math(hypo(A, B), M)
 => M = sqrt(A^2 + B^2).

%
% Operators
%
math(Flags, isin(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfx(Prec, isin, A, B).

math(Flags, notin(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfx(Prec, notin, A, B).

math(Flags, intersect(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, cap, A, B).

math(Flags, union(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, cup, A, B).

math(Flags, ':'(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, '#58', A, B).

math(Flags, kronecker(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, 'CircleTimes', A, B).

math(Flags, '=='(A, B), New, X)
 => math(Flags, A = B, New, X).

math(Flags, A = B, New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, =, A, B).

math(Flags, A \= B, New, X)
 => New = Flags,
    current_op(Prec, xfx, \=),
    X = xfx(Prec, ne, A, B).

math(Flags, A < B, New, X)
 => New = Flags,
    current_op(Prec, xfx, <),
    X = yfy(Prec, <, A, B).

math(Flags, A =< B, New, X)
 => New = Flags,
    current_op(Prec, xfx, =<),
    X = yfy(Prec, le, A, B).

math(Flags, ~(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, 'Tilde', A, B).

math(Flags, leftrightarrow(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, leftrightarrow, A, B).

math(Flags, iff(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, iff, A, B).

math(Flags, rightarrow(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, rightarrow, A, B).

math(Flags, rArr(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, rArr, A, B).

math(Flags, leftarrow(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, leftarrow, A, B).

math(Flags, lArr(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, lArr, A, B).

math(Flags, uparrow(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, uparrow, A, B).

math(Flags, uArr(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, uArr, A, B).

math(Flags, downarrow(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, downarrow, A, B).

math(Flags, dArr(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, dArr, A, B).

math(Flags, equiv(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, equiv, A, B).

math(Flags, cong(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, cong, A, B).

math(Flags, propto(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, propto, A, B).

math(Flags, A > B, New, X)
 => New = Flags,
    current_op(Prec, xfx, >),
    X = yfy(Prec, >, A, B).

math(Flags, A >= B, New, X)
 => New = Flags,
    current_op(Prec, xfx, >=),
    X = yfy(Prec, ge, A, B).

math(Flags, +A, New, X)
 => New = Flags,
    current_op(Prec, yfx, +),
    X = fy(Prec, +, A).

math(Flags, A + B, New, X)
 => New = Flags,
    current_op(Prec, yfx, +),
    X = yfy(Prec, +, A, B).

math(Flags, -A, New, X)
 => New = Flags,
    current_op(Prec, yfx, -),
    X = fy(Prec, -, A).

math(Flags, A - B, New, X)
 => New = Flags,
    current_op(Prec, yfx, -),
    X = yfy(Prec, -, A, B).

% Use dot or no dot instead of asterisk
math(Flags, A * B, New, X),
    type(Flags, A, TypeA),
    TypeA = atomic,
    type(Flags, B, TypeB),
    TypeB = atomic
 => New = Flags,
    X = nodot(A, B).

math(Flags, A * B, New, M)
 => New = Flags,
    M = cdot(A, B).

math(_Flags, '%*%'(A, B), M)
 => M = times(A, B).

math(_Flags, crossprod(A, B), M)
 => M = '%*%'(t(A), B).

math(_Flags, tcrossprod(A, B), M)
 => M = '%*%'(A, t(B)).

math(Flags, approx(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, approx, A, B).

math(Flags, ~(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, 'Tilde', A, B).

math(Flags, dot(A, B), New, X)
 => New = Flags,
    X = cdot(A, B).

math(Flags, cdot(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfy(Prec, cdot, A, B).

math(Flags, pm(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, +),
    X = yfy(Prec, pm, A, B).

math(Flags, nodot(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfy(Prec, '#x2062', A, B).

math(times(A, B), M)
 => current_op(Prec, yfx, *),
    M = yfy(Prec, times, A, B).

math(Flags, A / B, New, X)
 => New = Flags,
    current_op(Prec, yfx, /),
    X = yfx(Prec, /, A, B).

math(Flags, (A ; B), New, X)
 => New = Flags,
    current_op(Prec, xfy, ;),
    X = xfy(Prec, ;, A, B).

math(Flags, A^B, New, X)
 => New = Flags,
    X = superscript(A, B).

%
% Render
%
ml(fy(Prec, Op, A), M, Flags)
 => ml(op(Op), S, Flags),
    ml(right(Prec, A), X, Flags),
    M = mrow([S, X]).

ml(yf(Prec, Op, A), M, Flags)
 => ml(op(Op), S, Flags),
    ml(left(Prec, A), X, Flags),
    M = mrow([X, S]).

ml(xfx(Prec, Op, A, B), M, Flags)
 => ml(left(Prec-1, A), X, Flags),
    ml(op(Op), S, Flags),
    ml(right(Prec-1, B), Y, Flags),
    M = mrow([X, S, Y]).

ml(yfx(Prec, Op, A, B), M, Flags)
 => ml(left(Prec, A), X, Flags),
    ml(op(Op), S, Flags),
    ml(right(Prec-1, B), Y, Flags),
    M = mrow([X, S, Y]).

ml(xfy(Prec, Op, A, B), M, Flags)
 => ml(left(Prec-1, A), X, Flags),
    ml(op(Op), S, Flags),
    ml(right(Prec, B), Y, Flags),
    M = mrow([X, S, Y]).

ml(yfy(Prec, Op, A, B), M, Flags)
 => ml(left(Prec, A), X, Flags),
    ml(op(Op), S, Flags),
    ml(right(Prec, B), Y, Flags),
    M = mrow([X, S, Y]).

jax(Flags, fy(Prec, Op, A), M)
 => jax(Flags, op(Op), S),
    jax(Flags, right(Prec, A), X),
    format(string(M), "{~w}{~w}", [S, X]).

jax(Flags, yf(Prec, Op, A), M)
 => jax(Flags, op(Op), S),
    jax(Flags, left(Prec, A), X),
    format(string(M), "{~w}{~w}", [X, S]).

jax(Flags, xfx(Prec, Op, A, B), M)
 => jax(Flags, left(Prec-1, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec-1, B), Y),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

jax(Flags, yfx(Prec, Op, A, B), M)
 => jax(Flags, left(Prec, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec-1, B), Y),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

jax(Flags, xfy(Prec, Op, A, B), M)
 => jax(Flags, left(Prec-1, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec, B), Y),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

jax(Flags, yfy(Prec, Op, A, B), M)
 => jax(Flags, left(Prec, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec, B), Y),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

denoting(Flags, fy(_, _, A), D)
 => denoting(Flags, A, D).

denoting(Flags, yf(_, _, A), D)
 => denoting(Flags, A, D).

denoting(Flags, xfx(_, _, A, B), D)
 => denoting(Flags, A, DA),
    denoting(Flags, B, DB),
    append(DA, DB, D).

denoting(Flags, xfy(_, _, A, B), D)
 => denoting(Flags, A, DA),
    denoting(Flags, B, DB),
    append(DA, DB, D).

denoting(Flags, yfx(_, _, A, B), D)
 => denoting(Flags, A, DA),
    denoting(Flags, B, DB),
    append(DA, DB, D).

denoting(Flags, yfy(_, _, A, B), D)
 => denoting(Flags, A, DA),
    denoting(Flags, B, DB),
    append(DA, DB, D).

paren(Flags, fy(_, _, A), P)
 => paren(Flags, A, P).

paren(Flags, yf(_, _, A), P)
 => paren(Flags, A, P).

paren(Flags, xfx(_, _, A, B), P)
 => paren(Flags, A, PA),
    paren(Flags, B, PB),
    P is max(PA, PB).

paren(Flags, yfx(_, _, A, B), P)
 => paren(Flags, A, PA),
    paren(Flags, B, PB),
    P is max(PA, PB).

paren(Flags, xfy(_, _, A, B), P)
 => paren(Flags, A, PA),
    paren(Flags, B, PB),
    P is max(PA, PB).

paren(Flags, yfy(_, _, A, B), P)
 => paren(Flags, A, PA),
    paren(Flags, B, PB),
    P is max(PA, PB).

prec(_Flags, fy(Prec, _, _), P)
 => P = Prec.

prec(_Flags, yf(Prec, _, _), P)
 => P = Prec.

prec(_Flags, xfx(Prec, _, _, _), P)
 => P = Prec.

prec(_Flags, yfx(Prec, _, _, _), P)
 => P = Prec.

prec(_Flags, xfy(Prec, _, _, _), P)
 => P = Prec.

prec(_Flags, yfy(Prec, _, _, _), P)
 => P = Prec.

type(_Flags, fy(_, _, _), Type)
 => Type = op.

type(_Flags, yf(_, _, _), Type)
 => Type = op.

type(_Flags, xfx(_, _, _, _), Type)
 => Type = op.

type(_Flags, yfx(_, _, _, _), Type)
 => Type = op.

type(_Flags, xfy(_, _, _, _), Type)
 => Type = op.

type(_Flags, yfy(_, _, _, _), Type)
 => Type = op.


math(Flags, left(Prec, A), M),
    prec(Flags, A, P),
    P > Prec
 => M = paren(A).

math(_Flags, left(_, A), M)
 => M = A.

math(right(Prec, A), M)
 => P is Prec, % - 1,
    M = left(P, A).

denoting(Flags, left(_, A), D)
 => denoting(Flags, A, D).

denoting(Flags, right(_, A), D)
 => denoting(Flags, A, D).

%
% Named elements (see, e.g., integrate)
%
math(Flags, name(A, Name), New, M)
 => New = [name(Name) | Flags],
    M = A.

%
% Abbreviations
%
% with s^2_pool denoting the pooled variance
%
ml(denote(A, _, _), X, Flags)
 => ml(A, X, Flags).

jax(Flags, denote(A, _, _), X)
 => jax(Flags, A, X).

paren(Flags, denote(A, _, _), Paren)
 => paren(Flags, A, Paren).

prec(Flags, denote(A, _, _), Prec)
 => prec(Flags, A, Prec).

type(Flags, denote(A, _, _), Type)
 => type(Flags, A, Type).

denoting(Flags, denote(A, Expr, Info), Den)
 => denoting(Flags, Expr, T),
    Den = [denoting(A, Expr, Info) | T].

%
% Expand abbreviations
%
ml(denoting(A, Expr, Info), X, Flags)
 => ml(A = Expr, AExpr, Flags),
    X = span([math(AExpr), " denoting ", Info]).

jax(Flags, denoting(A, Expr, Info), X)
 => jax(Flags, A = Expr, AExpr),
    format(string(X), "$~w$ denoting ~w", [AExpr, Info]).

type(Flags, denoting(A, _, _), Type)
 => type(Flags, A, Type).

denoting(_Flags, denoting(_, _, _), Den)
 => Den = [].

%
% Collect abbreviations
%
ml(with(Abbreviations), X, Flags)
 => sort(Abbreviations, Sorted), % remove duplicates
    ml(with_(Sorted), X, Flags).

ml(with_([]), W, _Flags)
 => W = "".

ml(with_([A]), W, Flags)
 => ml(A, X, Flags),
    W = span([", with", &(nbsp), X]).

ml(with_([A, B | T]), W, Flags)
 => ml(A, X, Flags),
    ml(and([B | T]), Y, Flags),
    W = span([", with", &(nbsp), X | Y]).

ml(and([]), W, Flags)
 => W = ".".

ml(and([A | T]), W, Flags)
 => ml(A, X, Flags),
    ml(and(T), Y, Flags),
    W = span([", and", &(nbsp), X | Y]).

jax(Flags, with(Abbreviations), X)
 => sort(Abbreviations, Sorted), % remove duplicates
    jax(Flags, with_(Sorted), X).

jax(_Flags, with_([]), W)
 => W = "".

jax(Flags, with_([A]), W)
 => jax(Flags, A, X),
    format(string(W), ", with ~w", [X]).

jax(Flags, with_([A, B | T]), W)
 => jax(Flags, A, X),
    jax(Flags, and([B | T]), Y),
    format(string(W), ", with ~w~w", [X, Y]).

jax(_Flags, and([]), W)
 => W = ".".

jax(Flags, and([A | T]), W)
 => jax(Flags, A, X),
    jax(Flags, and(T), Y),
    format(string(W), ", and ~w~w", [X, Y]).

%
% No parentheses
%
math({}(A), M)
 => M = A.

%
% Parentheses
%
math('('(A), M)
 => M = paren(A).

ml(paren(A), M, Flags),
    paren(Flags, A, P),
    2 is P mod 3
 => ml(braces(A), M, Flags).

ml(paren(A), M, Flags),
    paren(Flags, A, P),
    1 is P mod 3
 => ml(brackets(A), M, Flags).

ml(paren(A), M, Flags)
 => ml(parentheses(A), M, Flags).

jax(Flags, paren(A), M),
    paren(Flags, A, P),
    2 is P mod 3
 => jax(Flags, braces(A), M).

jax(Flags, paren(A), M),
    paren(Flags, A, P),
    1 is P mod 3
 => jax(Flags, brackets(A), M).

jax(Flags, paren(A), M)
 => jax(Flags, parentheses(A), M).

paren(Flags, paren(A), P)
 => paren(Flags, A, P0),
    succ(P0, P).

type(_Flags, paren(_), T)
 => T = paren.

ml(parentheses(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo('('), X, mo(')')]).

jax(Flags, parentheses(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\left(~w\\right)", [X]).

paren(_Flags, parentheses(_), P)
 => P = 1.

type(_Flags, parentheses(_), T)
 => T = paren.

ml(brackets(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo('['), X, mo(']')]).

jax(Flags, brackets(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\left[~w\\right]", [X]).

paren(_Flags, brackets(_), P)
 => P = 2.

type(_Flags, brackets(_), T)
 => T = paren.

ml(braces(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo('{'), X, mo('}')]).

jax(Flags, braces(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\left\\{~w\\right\\}", [X]).

paren(_Flags, braces(_), P)
 => P = 3.

type(_Flags, braces(_), T)
 => T = paren.

%
% Lists of things
%
math(Flags, [H | T], New, M)
 => Flags = New,
    M = list(space, [H | T]).

ml(list(_, [A]), M, Flags)
 => ml(A, M, Flags).

ml(list(Sep, [A, B | T]), M, Flags)
 => ml(A, X, Flags),
    ml(tail(Sep, [B | T]), Y, Flags),
    M = mrow([X | Y]).

ml(tail(Sep, [A]), M, Flags)
 => ml(Sep, S, Flags),
    ml(A, X, Flags),
    M = [S, X].

ml(tail(Sep, [A, B | T]), M, Flags)
 => ml(Sep, S, Flags),
    ml(A, X, Flags),
    ml(tail(Sep, [B | T]), Y, Flags),
    M = [S, X | Y].

jax(Flags, list(_, [A]), M)
 => jax(Flags, A, M).

jax(Flags, list(Sep, [A, B | T]), M)
 => jax(Flags, A, X),
    jax(Flags, tail(Sep, [B | T]), Y),
    format(string(M), "{~w}{~w}", [X, Y]).

jax(Flags, tail(Sep, [A]), M)
 => jax(Flags, Sep, S),
    jax(Flags, A, X),
    format(string(M), "{~w}{~w}", [S, X]).

jax(Flags, tail(Sep, [A, B | T]), M)
 => jax(Flags, Sep, S),
    jax(Flags, A, X),
    jax(Flags, tail(Sep, [B | T]), Y),
    format(string(M), "{~w}{~w}{~w}", [S, X, Y]).

paren(Flags, list(_, List), P)
 => maplist(paren(Flags), List, P0),
    max_list(P0, P).

prec(Flags, list(_, [A]), P)
 => prec(Flags, A, P).

prec(Flags, list(Sep, [_, _ | _]), P)
 => prec(Flags, Sep, P).

denoting(Flags, list(_, L), D)
 => maplist(denoting(Flags), L, List),
    append(List, D).

%
% Fractions
%
ml(frac(N, D), M, Flags)
 => ml(N, X, Flags),
    ml(D, Y, Flags),
    M = mfrac([X, Y]).

jax(Flags, frac(N, D), M)
 => jax(Flags, N, X),
    jax(Flags, D, Y),
    format(string(M), "\\frac{~w}{~w}", [X, Y]).

paren(_Flags, frac(_, _), P)
 => P = 0.

prec(_Flags, frac(_, _), P)
 => current(P0, yfx, /),
    P is P0. % was - 1

type(_Flags, frac(_, _), Type)
  => Type = fraction.

%
% Large fraction
%
math(dfrac(Num, Den), M)
 => M = display(frac(Num, Den)).

math(over(Num, Den), M)
 => M = frac(Num, Den).

% Integer division
math(div(Num, Den), M)
 => M = floor(Num / Den).

%
% Integrate over range
%
ml(integrate(Fn, From, To, DX), M, Flags)
 => ml(Fn, XFn, Flags),
    ml(From, XFrom, Flags),
    ml(To, XTo, Flags),
    ml(DX, XDX, Flags),
    ml(space, Space, Flags),
    M = mrow([munderover([mo(&(int)), XFrom, XTo]), XFn, Space, mi(d), XDX]).

paren(Flags, integrate(_, _, _, A), Paren)
 => paren(Flags, A, Paren).

prec(_Flags, integrate(_, _, _, _), Prec)
 => current(Prec, yfx, +).

%
% Large font ("displaystyle")
%
ml(display(A), M, Flags)
 => ml(A, X, Flags),
    M = mstyle(displaystyle(true), X).

jax(Flags, display(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\displaystyle{~w}", [X]).

prec(Flags, display(A), P)
 => prec(Flags, A, P).

type(Flags, display(A), T)
 => type(Flags, A, T).

%
% Decorations
%
ml(overline(A), M, Flags)
 => ml(A, X, Flags),
    M = mover(accent(true), [X, mo(&(macr))]).

jax(Flags, overline(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\overline{~w}", [X]).

paren(Flags, overline(A), Paren)
 => paren(Flags, A, Paren).

% Put overline(x)^2 in parentheses
prec(_Flags, overline(_), Prec)
 => current(P, yfx, *),
    Prec = P.

type(Flags, overline(A), Type)
 => type(Flags, A, Type).

%
% Cancel out
%
ml(cancel(A), M, Flags)
 => ml(A, X, Flags),
    M = menclose(notation(updiagonalstrike), X).

jax(Flags, cancel(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\cancel{~w}", [X]).

paren(Flags, cancel(A), Paren)
 => paren(Flags, A, Paren).

prec(Flags, cancel(A), Prec)
 => prec(Flags, A, Prec).

type(Flags, cancel(A), Type)
 => type(Flags, A, Type).

%
% Box
%
ml(box(A), M, Flags)
 => ml(A, X, Flags),
    M = menclose(notation(roundedbox), X).

jax(Flags, box(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\boxed{~w}", [X]).

paren(Flags, box(A), Paren)
 => paren(Flags, A, Paren).

prec(Flags, box(A), Prec)
 => prec(Flags, A, Prec).

type(Flags, box(A), Type)
 => type(Flags, A, Type).

%
% Underbrace
%
ml(underbrace(A, U), M, Flags)
 => ml(A, X, Flags),
    ml(U, Y, Flags),
    M = munder([munder(accentunder(true),
                  [X, mo(stretchy(true), &('UnderBrace'))]), Y]).

paren(Flags, underbrace(A, _), Paren)
 => paren(Flags, A, Paren).

prec(Flags, underbrace(A, _), Prec)
 => prec(Flags, A, Prec).

type(Flags, underbrace(A, _), Type)
 => type(Flags, A, Type).

%
% Mistakes
%
math(Flags, omit_left(Expr), New, M),
    option(error(ignore), Flags, highlight)
 => Flags = New,
    M = Expr.

math(Flags, omit_left(Expr), New, M),
    option(error(asis), Flags, highlight),
    Expr =.. [_Op, _L, R]
 => Flags = New,
    M = R.

math(Flags, omit_left(Expr), New, M),
    option(error(fix), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [box(list(space, [L, op(Op)])), R]).

math(Flags, omit_left(Expr), New, M),
    option(error(highlight), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [cancel(list(space, [L, op(Op)])), R]).

math(Flags, omit_right(Expr), New, M),
    option(error(ignore), Flags, highlight)
 => Flags = New,
    M = Expr.

math(Flags, omit_right(Expr), New, M),
    option(error(asis), Flags, highlight),
    Expr =.. [_Op, L, _R]
 => Flags = New,
    M = L.

math(Flags, omit_right(Expr), New, M),
    option(error(fix), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [L, box(list(space, [op(Op), R]))]).

math(Flags, omit_right(Expr), New, M),
    option(error(highlight), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [L, cancel(list(space, [op(Op), R]))]).

math(Flags, omit(Expr), New, M),
    option(error(ignore), Flags, highlight)
 => Flags = New,
    M = Expr.

math(Flags, omit(Expr), New, M),
    option(error(fix), Flags, highlight)
 => Flags = New,
    M = box(Expr).

math(Flags, omit(Expr), New, M),
    option(error(highlight), Flags, highlight)
 => Flags = New,
    M = cancel(Expr).

math(Flags, add_left(Expr), New, M),
    option(error(ignore), Flags, highlight),
    Expr =.. [_Op, _L, R]
 => Flags = New,
    M = R.

math(Flags, add_left(Expr), New, M),
    option(error(asis), Flags, highlight)
 => Flags = New,
    M = Expr.

math(Flags, add_left(Expr), New, M),
    option(error(fix), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [cancel(list(space, [L, op(Op)])), R]).

math(Flags, add_left(Expr), New, M),
    option(error(highlight), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [box(list(space, [L, op(Op)])), R]).

math(Flags, add_right(Expr), New, M),
    option(error(ignore), Flags, highlight),
    Expr =.. [_Op, L, _R]
 => Flags = New,
    M = L.

math(Flags, add_right(Expr), New, M),
    option(error(asis), Flags, highlight)
 => Flags = New,
    M = Expr.

math(Flags, add_right(Expr), New, M),
    option(error(fix), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [L, cancel(list(space, [op(Op), R]))]).

math(Flags, add_right(Expr), New, M),
    option(error(highlight), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [L, box(list(space, [op(Op), R]))]).

math(Flags, omit(Expr), New, M),
    option(error(ignore), Flags, highlight)
 => Flags = New,
    M = Expr.

math(Flags, omit(_Expr), New, M),
    option(error(asis), Flags, highlight)
 => Flags = New,
    M = "". % suppress at the next level, in the list

math(Flags, omit(Expr), New, M),
    option(error(fix), Flags, highlight)
 => Flags = New,
    M = box(Expr).

math(Flags, omit(Expr), New, M),
    option(error(highlight), Flags, highlight)
 => Flags = New,
    M = cancel(Expr).

math(Flags, add(_Expr), New, M),
    option(error(ignore), Flags, highlight)
 => Flags = New,
    M = "". % suppress at the next level, in the list

math(Flags, add(Expr), New, M),
    option(error(asis), Flags, highlight)
 => Flags = New,
    M = Expr.

math(Flags, add(Expr), New, M),
    option(error(fix), Flags, highlight)
 => Flags = New,
    M = cancel(Expr).

math(Flags, add(Expr), New, M),
    option(error(highlight), Flags, highlight)
 => Flags = New,
    M = box(Expr).

math(Flags, instead(_Wrong, Correct), New, M),
    option(error(ignore), Flags, highlight)
 => Flags = New,
    M = Correct.

math(Flags, instead(Wrong, _Correct), New, M),
    option(error(asis), Flags, highlight)
 => Flags = New,
    M = Wrong.

math(Flags, instead(_Wrong, Correct), New, M),
    option(error(fix), Flags, highlight)
 => Flags = New,
    M = box(Correct).

math(Flags, instead(Wrong, Correct), New, M),
    option(error(highlight), Flags, highlight)
 => Flags = New,
    M = underbrace(Wrong, list(space, ["instead of", Correct])).

%
% Expert and buggy rules
%
math(Flags, expert(Flags, _, B), New, X)
 => New = Flags,
    X = B.

math(Flags, buggy(Flags, _, B), New, X)
 => New = Flags,
    X = B.

%
% Probability distributions
%
math(dbinom(K, N, Pi), M)
 => M = fn(subscript('P', "Bi"), (['X' = K] ; [N, Pi])).

math(pbinom(K, N, Pi), M)
 => M = fn(subscript('P', "Bi"), (['X' =< K] ; [N, Pi])).

math(qbinom(Alpha, N, Pi), M)
 => M = fn(subscript("arg min", k),
          [fn(subscript('P', "Bi"), (['X' =< k] ; [N, Pi])) > Alpha]).

math(dnorm(Z), M)
 => M = fn(phi, [Z]).

math(dnorm(X, Mu, Sigma2), M)
 => M = fn(phi, ([X] ; [Mu, Sigma2])).

math(pnorm(Z), M)
 => M = fn('Phi', [Z]).

math(pnorm(X, Mu, Sigma2), M)
 => M = fn('Phi', ([X] ; [Mu, Sigma2])).

math(qnorm(Alpha), M)
 => M = fn('Phi' ^ -1, [Alpha]).

math(qnorm(Alpha, Mu, Sigma2), M)
 => M = fn('Phi' ^ -1, ([Alpha] ; [Mu, Sigma2])).

math(pchisq(X, Df), M)
 => M = fn(subscript('F', fn(chi^2, [list(space, [Df, "df"])])), [X]).

math(qchisq(Alpha, Df), M)
 => M = fn(subscript('F' ^ -1, fn(chi^2, [list(space, [Df, "df"])])), [Alpha]).

math(pt(T, Df), M)
 => M = fn('P', (['T' =< T] ; [list(space, [Df, "df"])])).

math(qt(Alpha, Df), M)
 => M = fn(subscript('T', Alpha), [list(space, [Df, "df"])]).

%
% Functions like f(x) and f(x; a, b)
%
ml(fn(Name, (Args ; Pars)), M, Flags)
 => ml(Name, F, Flags),
    ml(paren(list(op(;), [list(op(','), Args), list(op(','), Pars)])), X, Flags),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, (Args ; Pars)), M),
    string(Name)
 => jax(Flags, Name, F),
    jax(Flags, paren(list(op(';'), [list(op(','), Args), list(op(','), Pars)])), X),
    format(string(M), "{~w}\\,{~w}", [F, X]).

jax(Flags, fn(Name, (Args ; Pars)), M)
 => jax(Flags, Name, F),
    jax(Flags, paren(list(op(';'), [list(op(','), Args), list(op(','), Pars)])), X),
    format(string(M), "{~w}{~w}", [F, X]).

paren(Flags, fn(_Name, (Args ; Pars)), Paren)
 => paren(Flags, list(op(','), Args), X),
    paren(Flags, list(op(','), Pars), Y),
    Paren is max(X, Y) + 1.

prec(Flags, fn(_Name, (_Args ; _Pars)), Prec)
 => prec(Flags, a * b, P0),
    Prec is P0 - 1.

type(_Flags, fn(_Name, (_Args ; _Pars)), Type)
 => Type = paren.

ml(fn(Name, [Arg]), M, Flags),
    type(Flags, Arg, paren)
 => ml(Name, F, Flags),
    ml(Arg, X, Flags),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, [Arg]), M),
    type(Flags, Arg, paren)
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w}{~w}", [F, X]).

% Omit parenthesis in special functions
ml(fn(Name, [Arg]), M, Flags),
    (   type(Flags, Name, special)
    ;   type(Flags, Name, subscript(N, _)), type(Flags, N, special)
    ;   type(Flags, Name, superscript(N, _)), type(Flags, N, special)
    ;   type(Flags, Name, subsupscript(N, _, _)), type(Flags, N, special)
    ),
    prec(Flags, Name, P),
    prec(Flags, Arg, Prec),
    P > Prec
 => ml(Name, F, Flags),
    ml(Arg, X, Flags),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, [Arg]), M),
    (   type(Flags, Name, special)
    ;   type(Flags, Name, subscript(N, _)), type(Flags, N, special)
    ;   type(Flags, Name, superscript(N, _)), type(Flags, N, special)
    ;   type(Flags, Name, subsupscript(N, _, _)), type(Flags, N, special)
    ),
    prec(Flags, Name, P),
    prec(Flags, Arg, Prec),
    P > Prec
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w}{~w}", [F, X]).

ml(fn(Name, [Arg]), M, Flags),
    type(Flags, Name, Type),
    member(Type, [special, subscript(_), superscript(_)]),
    prec(Flags, Arg, 0)
 => ml(Name, F, Flags),
    ml(Arg, X, Flags),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, [Arg]), M),
    type(Flags, Name, Type),
    member(Type, [special, subscript(_), superscript(_)]),
    prec(Flags, Arg, 0)
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w}{~w}", [F, X]).

ml(fn(Name, Args), M, Flags)
 => ml(Name, F, Flags),
    ml(paren(list(op(','), Args)), X, Flags),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, Args), M)
 => jax(Flags, Name, F),
    jax(Flags, paren(list(op(','), Args)), X),
    format(string(M), "{~w}{~w}", [F, X]).

paren(Flags, fn(_Name, [Arg]), P),
    type(Flags, Arg, paren)
 => paren(Flags, Arg, P).

paren(Flags, fn(_Name, [Arg]), P),
    prec(Flags, Arg, P0),
    P0 = 0
 => paren(Flags, Arg, P).

paren(Flags, fn(_Name, Args), P)
 => paren(Flags, list(op(','), Args), P).

prec(Flags, fn(Name, _Args), Prec),
    prec(Flags, Name, P),
    P = 0
 => current(Prec0, yfx, *),
    Prec is Prec0 - 1.

prec(Flags, fn(Name, _Args), Prec)
 => prec(Flags, Name, Prec).

type(_Flags, fn(_Name, _Args), Type)
 => Type = function.

% Default compounds
%
% Can't use the macros here because of left recursion
ml(A, M, Flags),
    compound(A),
    compound_name_arguments(A, N, Args)
 => ml(fn(N, Args), M, Flags).

jax(Flags, A, M),
    compound(A),
    compound_name_arguments(A, N, Args)
 => jax(Flags, fn(N, Args), M).

type(Flags, A, M),
    compound(A),
    compound_name_arguments(A, N, Args)
 => type(Flags, fn(N, Args), M).

%
% Defaults
%
math(A, M)
 => M = A.

math(_Flags, A, M)
 => M = A.

math(Flags, A, New, M)
 => New = Flags,
    M = A.

paren(Flags, A, P),
    math(A, M),
    dif(A, M)
 => paren(Flags, M, P).

paren(Flags, A, P),
    math(Flags, A, M),
    dif(A, M)
 => paren(Flags, M, P).

paren(Flags, A, P),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => paren(New, M, P).

paren(_, _, P)
 => P = 0.

prec(Flags, A, Den),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => prec(New, M, Den).

prec(_, _, P)
 => P = 0.

type(Flags, A, Type),
    math(A, M),
    dif(A, M)
 => type(Flags, M, Type).

type(Flags, A, Type),
    math(Flags, A, M),
    dif(A, M)
 => type(Flags, M, Type).

type(Flags, A, Type),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => type(New, M, Type).

type(_Flags, A, Type),
    compound(A)
 => Type = compound.

denoting(Flags, A, Den),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => denoting(New, M, Den).

denoting(Flags, Expression, Den),
    compound(Expression)
 => compound_name_arguments(Expression, _, Arguments),
    maplist(denoting(Flags), Arguments, List),
    append(List, Den).

% If everything fails, there is no abbreviation
denoting(_Flags, _, Den)
 => Den = [].

% Precedence
current(Prec, Fix, Op) :-
    atom(Op),
    current_op(Prec, Fix, Op).

