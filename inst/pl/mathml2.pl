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
mathml(R, M, Flags)
 => ml(R, M0, Flags),
    denoting(R, Denoting, Flags),
    ml(with(Denoting), With, Flags),
    !, M = [math(M0), With].

mathjax(R, M, Flags)
 => jax(R, M0, Flags),
    denoting(R, Denoting, Flags),
    jax(with(Denoting), With, Flags),
    !, format(string(M), "$~w$~w", [M0, With]).

%
% Macros
%
% macro(Flags, R, New, M): translates the R expression to another R
% expression M, checking for Flags and eventually changing Flags to New
%
% Calls math/2,3,4 macros
%
macro(R, New, M, Flags) :-
    math_hook(R, M0),
    !, New = Flags,
    M = M0.

macro(R, New, M, Flags) :-
    math(R, New, M, Flags),
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
jax(R, M, Flags),
    macro(R, New, R0, Flags)
 => jax(New, R0, M).

% Same for precedence checks, parentheses and types
prec(R, Prec, Flags),
    macro(R, New, R0, Flags)
 => prec(New, R0, Prec).

paren(R, Paren, Flags),
    macro(R, New, R0, Flags)
 => paren(New, R0, Paren).

type(R, Paren, Flags),
    macro(R, New, R0, Flags)
 => type(New, R0, Paren).

%
% Check if the current element is to be replaced
%
math(R, New, M, Flags),
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
math(Sum, New, M, Flags),
    compound(Sum),
    compound_name_arguments(Sum, sum, [Arg]),
    select(subscript(From), Flags, New0),
    select(superscript(To), New0, New1)
 => New = New1,
    M = fn(subsupscript(sum, From, To), [Arg]).

% Summation sign from (only subscript)
math(Sum, New, M, Flags),
    compound(Sum),
    compound_name_arguments(Sum, sum, [Arg]),
    select(subscript(Idx), Flags, New0)
 => New = New0,
    M = fn(subscript(sum, Idx), [Arg]).

mathml :-
    mathml(sum('['(x, i))).

% Same for product sign
math(Prod, New, M, Flags),
    compound(Prod),
    compound_name_arguments(Prod, prod, Arg),
    select(subscript(Idx), Flags, New0),
    select(superscript(Pwr), New0, New1)
 => New = New1,
    M = fn(subsupscript(prod, Idx, Pwr), Arg).

math(Prod, New, M, Flags),
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
math(subsupscript(R, Idx, Pwr), New, M, Flags)
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

jax(R, M, Flags),
    select(subscript(Idx), Flags, New0),
    select(superscript(Pwr), New0, New)
 => jax(New, R, X),
    jax(New, Idx, Y),
    jax(New, Pwr, Z),
    format(string(M), "{~w}_{~w}^{~w}", [X, Y, Z]).

type(R, Type, Flags),
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

math(subscript(R, Idx), New, M, Flags)
 => New = [subscript(Idx) | Flags],
    M = R.

% Render
ml(R, M, Flags),
    select(subscript(Idx), Flags, New)
 => ml(New, R, X),
    ml(New, Idx, Y),
    M = msub([X, Y]).

jax(R, M, Flags),
    select(subscript(Idx), Flags, New)
 => jax(New, R, X),
    jax(New, Idx, Y),
    format(string(M), "{~w}_{~w}", [X, Y]).

type(R, Type, Flags),
    member(subscript(Idx), Flags)
 => Type = subscript(R, Idx).

mathml :-
    mathml(subscript(x, i)).

%
% Superscripts like s^D
%
math(superscript(R, Idx), New, M, Flags)
 => New = [superscript(Idx) | Flags],
    M = R.

% Render
ml(R, M, Flags),
    select(superscript(Pwr), Flags, New),
    prec(R, Prec, Flags),
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

jax(R, M, Flags),
    prec(R, Prec, Flags),
    current_op(P, xfy, ^),
    Prec >= P,
    select(superscript(Pwr), Flags, New)
 => jax(New, paren(R), X),
    jax(New, Pwr, Y),
    format(string(M), "{~w}^{~w}", [X, Y]).

jax(R, M, Flags),
    select(superscript(Pwr), Flags, New)
 => jax(New, R, X),
    jax(New, Pwr, Y),
    format(string(M), "{~w}^{~w}", [X, Y]).

type(R, Type, Flags),
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

jax(text(R), M, _Flags)
 => format(string(M), "\\mathrm{~w}", [R]).

type(text(_), T, _Flags)
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

jax(greek(R), M, _Flags)
 => format(string(M), "\\~w", [R]).

type(greek(_), T, _Flags)
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

jax(boolean(R), M, _Flags)
 => format(string(M), "~w", [R]).

type(boolean(_), T, _Flags)
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

jax(set(empty), M, _Flags)
 => M = "\\emptyset".

type(set(empty), T, _Flags)
 => T = atomic.

%
% sin^2(x) etc.
%
math(sin(A), New, M, Flags),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(sin^Pwr, [A]).

math(sinpi(A), New, M, Flags),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(sinpi^Pwr, [A]).

math(cos(A), New, M, Flags),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(cos^Pwr, [A]).

math(cospi(A), New, M, Flags),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(cospi^Pwr, [A]).

math(tan(A), New, M, Flags),
    select(superscript(Pwr), Flags, Flags1)
 => New = Flags1,
    M = fn(tan^Pwr, [A]).

math(tanpi(A), New, M, Flags),
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

prec(special(sum), Prec, _Flags)
 => current(P, yfx, *),
    Prec is P + 1.

ml(special(prod), M, _Flags)
 => M = mo(&(prod)).

prec(special(prod), Prec, _Flags)
 => current(P, yfx, *),
    Prec is P + 1.

ml(special(R), M, _Flags)
 => M = mi(R).

jax(special(sum), M, _Flags)
 => M = "{\\sum}".

jax(special(prod), M, _Flags)
 => M = "{\\prod}".

jax(special(sgn), M, _Flags)
 => M = "{\\mathrm{sgn}\\,}".

jax(special(R), M, _Flags)
 => format(string(M), "\\~w", [R]).

type(special(_), T, _Flags)
 => T = special.

prec(special(_), Prec, _Flags)
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

jax(space(thinmathspace), M, _Flags)
 => M = "\\,".

jax(space(_Width), M, _Flags)
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

jax(ident(R), M, Flags),
    member(mathvariant(calligraphy), Flags)
 => format(string(M), "\\mathcal{~w}", [R]).

jax(ident(R), M, _Flags)
 => format(string(M), "~w", [R]).

type(ident(_), T, _Flags)
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

jax(abs(R), M, Flags)
 => jax(R, X, Flags),
    format(string(M), "\\left|~w\\right|", [X]).

paren(abs(_), P, _Flags)
 => P = 0.

prec(abs(R), P, Flags)
 => prec(paren(R), P, Flags).

math(sign(R), M)
 => M = fn(sgn, [R]).

ml(sqrt(R), M, Flags)
 => ml(R, X, Flags),
    M = msqrt(X).

jax(sqrt(A), M, Flags)
 => jax(A, X, Flags),
    format(string(M), "\\sqrt{~w}", [X]).

paren(sqrt(_), P, _Flags)
 => P = 0.

prec(sqrt(_), P, _Flags)
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

jax(forall(A), M, Flags)
 => jax(A, X, Flags),
    format(string(M), "\\forall{~w}", [X]).

paren(forall(A), P, Flags)
 => paren(A, P, Flags).

prec(forall(_), P, _Flags)
 => current(P, yfx, *).

math(any(A), M)
 => M = exists(A).

ml(exists(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo(&('Exists')), mo(&(af)), X]).

jax(exists(A), M, Flags)
 => jax(A, X, Flags),
    format(string(M), "\\exists{~w}", [X]).

paren(exists(A), P, Flags)
 => paren(A, P, Flags).

prec(exists(_), P, _Flags)
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

jax(choose(N, K), M, Flags)
 => jax(N, X, Flags),
    jax(K, Y, Flags),
    format(string(M), "\\binom{~w}{~w}", [X, Y]).

paren(choose(_, _), P, _Flags)
 => P = 1.

prec(choose(_, _), P, _Flags)
 => P = 0.

type(choose(_, _), T, _Flags)
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

jax(ceiling(A), M, Flags)
 => jax(A, X, Flags),
    format(string(M), "\\lceil{~w}\\rceil", [X]).

paren(ceiling(_), P, _Flags)
 => P is 0.

ml(floor(A), M, Flags)
 => ml(A, X, Flags),
    M = mrow([mo(&(lfloor)), X, mo(&(rfloor))]).

jax(floor(A), M, Flags)
 => jax(A, X, Flags),
    format(string(M), "\\lfloor{~w}\\rfloor", [X]).

paren(floor(_), P, _Flags)
 => P is 0.

math((_F :- Body), M, _Flags)
 => M = Body.

math(Function, M, _Flags),
    compound(Function),
    compound_name_arguments(Function, function, Args)
 => M = lambda(Args).

math(lambda(Args), M, Flags),
    member(name-N, Flags)
 => M = fn(N, Args).

math(lambda(Args), M, Flags)
 => option(name(N), Flags),
    M = fn(N, Args).

math(lambda(Args), M, _Flags)
 => M = fn(lambda, Args).

type(lambda(_), T, _Flags)
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

jax(body([R]), M, Flags)
 => jax(R, M, Flags).

jax(body(Body), M, Flags)
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

ml_row(Row, M, Flags),
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

jax_row(Row, M, Flags),
    compound(Row),
    compound_name_arguments(Row, Name, Cells),
    member(Name, ['#', '$$', '%', '!'])
 => maplist(jax_cell(Flags), Cells, C),
    atomic_list_concat(C, ' & ', R),
    format(string(M), "~w\\\\\n", [R]).

jax_cell(C, M, Flags)
 => jax(C, X, Flags),
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

jax(ifelse(T, Y, N), M, Flags)
 => jax(T, Test, Flags),
    jax(Y, Yes, Flags),
    jax(N, No, Flags),
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

jax(if(T, Y), M, Flags)
 => jax(T, Test, Flags),
    jax(Y, Yes, Flags),
    format(string(M), "{~w},\\ \\mathrm{if}\\ {~w}", [Yes, Test]).

paren(if(_, _), P, _Flags)
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

jax(prod(A), M, Flags)
 => jax(A, X, Flags),
    format(string(M), "\\prod{~w}", [X]).

jax(Prod, M, Flags),
    compound(Prod),
    compound_name_arguments(Prod, prod, Args)
 => maplist(jax(Flags), Args, X),
    format(string(M), "\\prod{~w}", [X]).

paren(Prod, P, Flags),
    compound(Prod),
    compound_name_arguments(Prod, prod, Args)
 => maplist(paren(Flags), Args, PX),
    max_list(PX, P).

prec(Prod, P, _Flags),
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
math(cal(A), New, M, Flags)
 => New = [mathvariant(calligraphy) | Flags],
    M = A.

%
% Extract value from a result (e.g., integrate)
%
math($(Fn, "value"), New, M, Flags)
 => Flags = New,
    M = Fn.

%
% Integrate over range
%
% Case A: Fn is a function
math(integrate(Fn, Lower, Upper), New, M, Flags),
    Fn = (Head :- _Body),
    compound(Head),
    compound_name_arguments(Head, function, [DX | _]),
    member(name-Name, Flags)
 => Flags = New,
    M = integrate(fn(Name, [DX]), Lower, Upper, DX).

math(integrate(Fn, Lower, Upper), New, M, Flags),
    Fn = (Head :- _Body),
    compound(Head),
    compound_name_arguments(Head, function, [DX | _])
 => Flags = New,
    M = integrate(fn(lambda, [DX]), Lower, Upper, DX).

% Case B: Fn is an atom
math(integrate(Fn, Lower, Upper), New, M, Flags),
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

jax(integrate(Fn, From, To, DX), M, Flags)
 => jax(Fn, XFn, Flags),
    jax(From, XFrom, Flags),
    jax(To, XTo, Flags),
    jax(DX, XDX, Flags),
    format(string(M), "\\int_{~w}^{~w}{~w}\\,{d{~w}}", [XFrom, XTo, XFn, XDX]).

paren(integrate(_, _, _, A), Paren, Flags)
 => paren(A, Paren, Flags).

prec(integrate(_, _, _, _), Prec, _Flags)
 => current(Prec, yfx, *).

% hats
ml(hat(A), M, Flags)
 => ml(A, X, Flags),
    M = mover(accent(true), [X, mo(&('Hat'))]).

jax(hat(A), M, Flags)
 => jax(A, X, Flags),
    format(string(M), "\\hat{~w}", [X]).

paren(hat(A), Paren, Flags)
 => paren(A, Paren, Flags).

prec(hat(A), Prec, Flags)
 => prec(A, Prec, Flags).

type(hat(A), Type, Flags)
 => type(A, Type, Flags).

% tilde
ml(tilde(A), M, Flags)
 => ml(A, X, Flags),
	M = mover(accent(true), [X, mo(&(tilde))]).

paren(tilde(A), Paren, Flags)
 => paren(A, Paren, Flags).

prec(tilde(A), Prec, Flags)
 => prec(A, Prec, Flags).

type(tilde(A), Type, Flags)
 => type(A, Type, Flags).

%
% Mathematical signs
%
ml(op(le), M, _Flags)
 => M = mo(&(le)).

jax(op(le), M, _Flags)
 => M = "\\le".

ml(op(ge), M, _Flags)
 => M = mo(&(ge)).

jax(op(ge), M, _Flags)
 => M = "\\ge".

ml(op(ne), M, _Flags)
 => M = mo(&(ne)).

jax(op(ne), M, _Flags)
 => M = "\\ne".

ml(op(cdot), M, _Flags)
 => M = mo(&(sdot)).

jax(op(cdot), M, _Flags)
 => M = "\\cdot".

ml(op(pm), M, _Flags)
 => M = mo(&(pm)).

jax(op(pm), M, _Flags)
 => M = "\\pm".

ml(op(times), M, _Flags)
 => M = mo(&(times)).

jax(op(times), M, _Flags)
 => M = "\\times".

ml(op(sum), M, _Flags)
 => M = mo(&(sum)).

jax(op(sum), M, _Flags)
 => M = "\\sum".

ml(op(prod), M, _Flags)
 => M = mo(&(prod)).

jax(op(prod), M, _Flags)
 => M = "\\prod".

ml(op('#58'), M, _Flags)
 => M = mo(&('#58')).

jax(op('#58'), M, _Flags)
 => M = ":".

ml(op(','), M, _Flags)
 => M = mo(',').

jax(op(','), M, _Flags)
 => M = ",".

ml(op('CircleTimes'), M, _Flags)
 => M = mo(&('CircleTimes')).

jax(op('CircleTimes'), M, _Flags)
 => M = "\\otimes".

ml(op('#x2062'), M, _Flags)
 => M = mo(&('#x2062')).

jax(op('#x2062'), M, _Flags)
 => M = "{}".

ml(op('Tilde'), M, _Flags)
 => M = mo(&('Tilde')).

jax(op('Tilde'), M, _Flags)
 => M = "\\sim".

ml(op(leftrightarrow), M, _Flags)
 => M = mo(&(leftrightarrow)).

jax(op(leftrightarrow), M, _Flags)
 => M = "\\leftrightarrow".

ml(op(iff), M, _Flags)
 => M = mo(&(iff)).

jax(op(iff), M, _Flags)
 => M = "\\iff".

ml(op(rightarrow), M, _Flags)
 => M = mo(&(rightarrow)).

jax(op(rightarrow), M, _Flags)
 => M = "\\rightarrow".

ml(op(rArr), M, _Flags)
 => M = mo(&(rArr)).

jax(op(rArr), M, _Flags)
 => M = "\\Rightarrow".

ml(op(leftarrow), M, _Flags)
 => M = mo(&(leftarrow)).

jax(op(leftarrow), M, _Flags)
 => M = "\\leftarrow".

ml(op(lArr), M, _Flags)
 => M = mo(&(lArr)).

jax(op(lArr), M, _Flags)
 => M = "\\Leftarrow".

ml(op(uparrow), M, _Flags)
 => M = mo(&(uparrow)).

jax(op(uparrow), M, _Flags)
 => M = "\\uparrow".

ml(op(uArr), M, _Flags)
 => M = mo(&(uArr)).

jax(op(uArr), M, _Flags)
 => M = "\\Uparrow".

ml(op(downarrow), M, _Flags)
 => M = mo(&(downarrow)).

jax(op(downarrow), M, _Flags)
 => M = "\\downarrow".

ml(op(dArr), M, _Flags)
 => M = mo(&(dArr)).

jax(op(dArr), M, _Flags)
 => M = "\\Downarrow".

ml(op(approx), M, _Flags)
 => M = mo(&(approx)).

jax(op(approx), M, _Flags)
 => M = "\\approx".

ml(op(equiv), M, _Flags)
 => M = mo(&(equiv)).

jax(op(equiv), M, _Flags)
 => M = "\\equiv".

ml(op(cong), M, _Flags)
 => M = mo(&(cong)).

jax(op(cong), M, _Flags)
 => M = "\\cong".

ml(op(propto), M, _Flags)
 => M = mo(&(prop)).

jax(op(propto), M, _Flags)
 => M = "\\propto".

ml(op(and), M, _Flags)
 => M = mo(&(and)).

jax(op(and), M, _Flags)
 => M = "\\land".

ml(op(or), M, _Flags)
 => M = mo(&(or)).

jax(op(or), M, _Flags)
 => M = "\\lor".

ml(op(not), M, _Flags)
 => M = mo(&(not)).

jax(op(not), M, _Flags)
 => M = "\\lnot".

ml(op(veebar), M, _Flags)
 => M = mo(&(veebar)).

jax(op(veebar), M, _Flags)
 => M = "\\veebar".

ml(op(isin), M, _Flags)
 => M = mo(&(isin)).

jax(op(isin), M, _Flags)
 => M = "\\in".

ml(op(notin), M, _Flags)
 => M = mo(&(notin)).

jax(op(notin), M, _Flags)
 => M = "\\notin".

ml(op(cap), M, _Flags)
 => M = mo(&(cap)).

jax(op(cap), M, _Flags)
 => M = "\\cap".

ml(op(cup), M, _Flags)
 => M = mo(&(cup)).

jax(op(cup), M, _Flags)
 => M = "\\cup".

ml(op(A), M, _Flags)
 => M = mo(A).

jax(op(A), M, _Flags)
 => format(string(M), "~w", [A]).

prec(op(A), P, _Flags),
    current(P0, _Fix, A)
 => P = P0.

current(0, fy, op(sum)).

denoting(op(_), D, _Flags)
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

jax(posint(A), M, _Flags)
 => format(string(M), "~w", [A]).

jax(pos(1.0Inf), M, _Flags)
 => M = "\\infty".

jax(pos(A), M, Flags)
 => option(round(D), Flags, 2),
    format(atom(Mask), '~~~wf', [D]),
    format(string(M), Mask, [A]).

type(pos(_), Type, _Flags)
 => Type = atomic.

type(posint(_), Type, _Flags)
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
math(isin(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfx(Prec, isin, A, B).

math(notin(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfx(Prec, notin, A, B).

math(intersect(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, cap, A, B).

math(union(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, cup, A, B).

math(':'(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, '#58', A, B).

math(kronecker(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, 'CircleTimes', A, B).

math('=='(A, B), New, X, Flags)
 => math(A = B, New, X, Flags).

math(A = B, New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, =, A, B).

math(A \= B, New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, \=),
    X = xfx(Prec, ne, A, B).

math(A < B, New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, <),
    X = yfy(Prec, <, A, B).

math(A =< B, New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =<),
    X = yfy(Prec, le, A, B).

math(~(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, 'Tilde', A, B).

math(leftrightarrow(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, leftrightarrow, A, B).

math(iff(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, iff, A, B).

math(rightarrow(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, rightarrow, A, B).

math(rArr(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, rArr, A, B).

math(leftarrow(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, leftarrow, A, B).

math(lArr(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, lArr, A, B).

math(uparrow(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, uparrow, A, B).

math(uArr(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, uArr, A, B).

math(downarrow(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, downarrow, A, B).

math(dArr(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ->),
    X = yfy(Prec, dArr, A, B).

math(equiv(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, equiv, A, B).

math(cong(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, cong, A, B).

math(propto(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, propto, A, B).

math(A > B, New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, >),
    X = yfy(Prec, >, A, B).

math(A >= B, New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, >=),
    X = yfy(Prec, ge, A, B).

math(+A, New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, +),
    X = fy(Prec, +, A).

math(A + B, New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, +),
    X = yfy(Prec, +, A, B).

math(-A, New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, -),
    X = fy(Prec, -, A).

math(A - B, New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, -),
    X = yfy(Prec, -, A, B).

% Use dot or no dot instead of asterisk
math(A * B, New, X, Flags),
    type(A, TypeA, Flags),
    TypeA = atomic,
    type(B, TypeB, Flags),
    TypeB = atomic
 => New = Flags,
    X = nodot(A, B).

math(A * B, New, M, Flags)
 => New = Flags,
    M = cdot(A, B).

math('%*%'(A, B), M, _Flags)
 => M = times(A, B).

math(crossprod(A, B), M, _Flags)
 => M = '%*%'(t(A), B).

math(tcrossprod(A, B), M, _Flags)
 => M = '%*%'(A, t(B)).

math(approx(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, approx, A, B).

math(~(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, 'Tilde', A, B).

math(dot(A, B), New, X, Flags)
 => New = Flags,
    X = cdot(A, B).

math(cdot(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfy(Prec, cdot, A, B).

math(pm(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, +),
    X = yfy(Prec, pm, A, B).

math(nodot(A, B), New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfy(Prec, '#x2062', A, B).

math(times(A, B), M)
 => current_op(Prec, yfx, *),
    M = yfy(Prec, times, A, B).

math(A / B, New, X, Flags)
 => New = Flags,
    current_op(Prec, yfx, /),
    X = yfx(Prec, /, A, B).

math((A ; B), New, X, Flags)
 => New = Flags,
    current_op(Prec, xfy, ;),
    X = xfy(Prec, ;, A, B).

math(A^B, New, X, Flags)
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

jax(fy(Prec, Op, A), M, Flags)
 => jax(op(Op), S, Flags),
    jax(right(Prec, A), X, Flags),
    format(string(M), "{~w}{~w}", [S, X]).

jax(yf(Prec, Op, A), M, Flags)
 => jax(op(Op), S, Flags),
    jax(left(Prec, A), X, Flags),
    format(string(M), "{~w}{~w}", [X, S]).

jax(xfx(Prec, Op, A, B), M, Flags)
 => jax(left(Prec-1, A), X, Flags),
    jax(op(Op), S, Flags),
    jax(right(Prec-1, B), Y, Flags),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

jax(yfx(Prec, Op, A, B), M, Flags)
 => jax(left(Prec, A), X, Flags),
    jax(op(Op), S, Flags),
    jax(right(Prec-1, B), Y, Flags),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

jax(xfy(Prec, Op, A, B), M, Flags)
 => jax(left(Prec-1, A), X, Flags),
    jax(op(Op), S, Flags),
    jax(right(Prec, B), Y, Flags),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

jax(yfy(Prec, Op, A, B), M, Flags)
 => jax(left(Prec, A), X, Flags),
    jax(op(Op), S, Flags),
    jax(right(Prec, B), Y, Flags),
    format(string(M), "{~w}{~w}{~w}", [X, S, Y]).

denoting(fy(_, _, A), D, Flags)
 => denoting(A, D, Flags).

denoting(yf(_, _, A), D, Flags)
 => denoting(A, D, Flags).

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

