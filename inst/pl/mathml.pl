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
r2mathml(R, S) :-
    r2mathml([], R, S).

r2mathml(Flags, R, S) :-
    mathml(Flags, R, M),
    html(M, H, []),
    maplist(atom_string, H, S).

% Same for MathJax/LaTeX
r2mathjax(R, S) :-
    r2mathjax([], R, S).

r2mathjax(Flags, R, S) :-
    mathjax(Flags, R, S).

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

macro(Flags, R, Flags, M) :-
    math(R, M),
    dif(R, M).

macro(Flags, R, Flags, M) :-
    math(Flags, R, M),
    dif(R, M).

macro(Flags, R, New, M) :-
    math(Flags, R, New, M),
    dif(Flags-R, New-M).

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
ml(Flags, R, M),
    macro(Flags, R, New, R0)
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
    compound_name_arguments(Sum, sum, Arg),
    select(subscript(From), Flags, New0),
    select(superscript(To), New0, New1)
 => New = New1,
    M = fn(subsupscript(sum, From, To), Arg).

% Summation sign from (only subscript)
math(Flags, Sum, New, M),
    compound(Sum),
    compound_name_arguments(Sum, sum, Arg),
    select(subscript(Idx), Flags, New0)
 => New = New0,
    M = fn(subscript(sum, Idx), Arg).

math(Flags, Sum, New, M),
    compound(Sum),
    compound_name_arguments(Sum, sum, Arg)
 => New = Flags,
    M = fn(sum, Arg).

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

math(Flags, Prod, New, M),
    compound(Prod),
    compound_name_arguments(Prod, prod, Arg),
    select(superscript(Pwr), Flags, New1)
 => New = New1,
    M = fn(superscript(prod, Pwr), Arg).

math(Flags, Prod, New, M),
    compound(Prod),
    compound_name_arguments(Prod, prod, Arg)
 => New = Flags,
    M = fn(prod, Arg).

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
ml(Flags, R, M),
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
    format(string(M), "{~w_~w^~w}", [X, Y, Z]).

mathml :-
    mathml(subsupscript(x, i, 2)).

%
% Subscript like s_D
%
math('['(R, Idx), M)
 => M = subscript(R, Idx).

math(Flags, subscript(R, Idx), New, M)
 => New = [subscript(Idx) | Flags],
    M = R.

% Render
ml(Flags, R, M),
    select(subscript(Idx), Flags, New)
 => ml(New, R, X),
    ml(New, Idx, Y),
    M = msub([X, Y]).

jax(Flags, R, M),
    select(subscript(Idx), Flags, New)
 => jax(New, R, X),
    jax(New, Idx, Y),
    format(string(M), "{~w_~w}", [X, Y]).

mathml :-
    mathml(subscript(x, i)).

%
% Superscripts like s^D
%
math(Flags, superscript(R, Idx), New, M)
 => New = [superscript(Idx) | Flags],
    M = R.

% Render
ml(Flags, R, M),
    select(superscript(Pwr), Flags, New)
 => ml(New, R, X),
    ml(New, Pwr, Y),
    M = msup([X, Y]).

jax(Flags, R, M),
    select(superscript(Pwr), Flags, New)
 => jax(New, R, X),
    jax(New, Pwr, Y),
    format(string(M), "{~w^~w}", [X, Y]).

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

ml(_Flags, text(R), M)
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
        theta, vartheta, iota, kappa, lambda, mu, nu, xi, pi, rho, sigma, tau,
        upsilon, phi, varphi, chi, psi, omega, 'Gamma', 'Delta', 'Theta',
        'Lambda', 'Xi', 'Pi', 'Sigma', 'Upsilon', 'Phi', 'Psi', 'Omega'])
 => M = greek(R).

ml(_Flags, greek(R), M)
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

ml(_Flags, boolean(R), M)
 => M = mi(R).

jax(_Flags, boolean(R), M)
 => format(string(M), "{~w}", [R]).

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

ml(_Flags, set(empty), M)
 => M = mi(&(empty)).

jax(_Flags, set(empty), M)
 => M = "\\emptyset".

type(_Flags, set(empty), T)
 => T = atomic.

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

ml(_Flags, special(sum), M)
 => M = mo(&(sum)).

prec(_Flags, special(sum), Prec)
 => current(P, yfx, *),
    Prec is P + 1.

prec(_Flags, special(_), Prec)
 => current(Prec, yfx, *).

ml(_Flags, special(prod), M)
 => M = mo(&(prod)).

ml(_Flags, special(R), M)
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

mathml :-
    mathml(exp(x)),
    mathml(exp(x + y)).

%
% Space
%
math(space, M)
 => M = space(thinmathspace).

ml(_Flags, space(W), M)
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

ml(Flags, ident(R), M),
    member(mathvariant(calligraphy), Flags)
 => M = mi(mathvariant(script), R).

ml(_Flags, ident(R), M)
 => M = mi(R).

jax(Flags, ident(R), M),
    member(mathvariant(calligraphy), Flags)
 => format(string(M), "\\mathcal{~w}", [R]).

jax(_Flags, ident(R), M)
 => format(string(M), "{~w}", [R]).

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

ml(Flags, abs(R), M)
 => ml(Flags, R, X),
    M = mrow([mo(&(vert)), X, mo(&(vert))]).

jax(Flags, abs(R), M)
 => jax(Flags, R, X),
    format(string(M), "\\left|{~w}\\right|", [X]).

paren(_Flags, abs(_), P)
 => P = 0.

prec(Flags, abs(R), P)
 => prec(Flags, paren(R), P).

math(sign(R), M)
 => M = fn(sgn, [R]).

ml(Flags, sqrt(R), M)
 => ml(Flags, R, X),
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

ml(Flags, forall(A), M)
 => ml(Flags, A, X),
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

ml(Flags, exists(A), M)
 => ml(Flags, A, X),
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

ml(Flags, choose(N, K), M)
 => ml(Flags, N, X),
    ml(Flags, K, Y),
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

ml(Flags, ceiling(A), M)
 => ml(Flags, A, X),
    M = mrow([mo(&(lceil)), X, mo(&(rceil))]).

jax(Flags, ceiling(A), M)
 => jax(Flags, A, X),
    format(string(M), "{\\lceil~w\\rceil}", [X]).

paren(_Flags, ceiling(_), P)
 => P is 0.

ml(Flags, floor(A), M)
 => ml(Flags, A, X),
    M = mrow([mo(&(lfloor)), X, mo(&(rfloor))]).

jax(Flags, floor(A), M)
 => jax(Flags, A, X),
    format(string(M), "{\\lfloor~w\\rfloor}", [X]).

paren(_Flags, floor(_), P)
 => P is 0.

math(_Flags, (F :- Body), M)
 => M = (F == Body).

math(_Flags, Function, M),
    compound(Function),
    compound_name_arguments(Function, '$function', Args)
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

math(body([R]), M)
 => M = R.

math(body(Body), M)
 => findall([A], member(A, Body), List),
    M = rbind(List).

% Matrices
ml(Flags, rbind(R), M)
 => maplist(ml_row(Flags), R, Rows),
    M = mrow([mo('{'), mtable(columnalign(left), Rows)]).

ml_row(Flags, R, M)
 => maplist(ml_cell(Flags), R, Cells),
    M = mtr(Cells).

ml_cell(Flags, C, M)
 => ml(Flags, C, X),
    M = mtd(X).

jax(Flags, rbind([H | R]), M)
 => findall(l, member(_, H), Ls),
    atomic_list_concat(Ls, LLL),
    maplist(jax_row(Flags), [H | R], Rows),
    atomic_list_concat(Rows, Lines),
    format(string(M), "\\left\\{\\begin{array}{~w}~w\\end{array}\\right.", [LLL, Lines]).

jax_row(Flags, R, M)
 => maplist(jax_cell(Flags), R, Cells),
    atomic_list_concat(Cells, ' & ', Row),
    format(string(M), "~w\\\\\n", [Row]).

jax_cell(Flags, C, M)
 => jax(Flags, C, X),
    format(string(M), "~w", [X]).

math(Identical, M),
    compound(Identical),
    compound_name_arguments(Identical, identical, [X, Y])
 => M = (X == Y).

ml(Flags, ifelse(T, Y, N), M)
 => ml(Flags, T, Test),
    ml(Flags, Y, Yes),
    ml(Flags, N, No),
    ml(Flags, space, S),
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
      "\\left\\{\\begin{array}{ll} ~w & \\mathrm{if}~~~w\\\\ ~w & \\mathrm{otherwise}\\end{array}\\right.",
      [Yes, Test, No]).

paren(_Flags, ifelse(_, _, _), P)
 => P is 0.

ml(Flags, if(T, Y), M)
 => ml(Flags, T, Test),
    ml(Flags, Y, Yes),
    ml(Flags, space, S),
    M = mrow([Yes, mtext(","), S, mtext("if"), S, Test]).

jax(Flags, if(T, Y), M)
 => jax(Flags, T, Test),
    jax(Flags, Y, Yes),
    format(string(M), "~w,\\ \\mathrm{if}\\ ~w", [Yes, Test]).

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

ml(Flags, Prod, M),
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
% Integrate over range
%
% extract value
math(Flags, $(integrate(Fn, Lower, Upper), value), New, M)
 => Flags = New,
    M = integrate(Fn, Lower, Upper).

% No argument names
math(Flags, integrate(Fn, Lower, Upper), New, M)
 => Flags = New,
    r_eval('['(formalArgs(args(Fn)), 1), Arg1),
    atom_string(DX, Arg1),
    M = integrate(fn(Fn, [DX]), Lower, Upper, DX).

% Internal
ml(Flags, integrate(Fn, From, To, DX), M)
 => ml(Flags, Fn, XFn),
    ml(Flags, From, XFrom),
    ml(Flags, To, XTo),
    ml(Flags, DX, XDX),
    ml(Flags, space, Space),
    M = mrow([munderover([mo(&(int)), XFrom, XTo]), XFn, Space, mi(d), XDX]).

paren(Flags, integrate(_, _, _, A), Paren)
 => paren(Flags, A, Paren).

prec(_Flags, integrate(_, _, _, _), Prec)
 => current(Prec, yfx, *).

% hats
ml(Flags, hat(A), M)
 => ml(Flags, A, X),
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
ml(Flags, tilde(A), M)
 => ml(Flags, A, X),
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
ml(_Flags, op(le), M)
 => M = mo(&(le)).

jax(_Flags, op(le), M)
 => M = "\\le".

ml(_Flags, op(ge), M)
 => M = mo(&(ge)).

jax(_Flags, op(ge), M)
 => M = "\\ge".

ml(_Flags, op(ne), M)
 => M = mo(&(ne)).

jax(_Flags, op(ne), M)
 => M = "\\ne".

ml(_Flags, op(cdot), M)
 => M = mo(&(sdot)).

jax(_Flags, op(cdot), M)
 => M = "\\cdot".

ml(_Flags, op(pm), M)
 => M = mo(&(pm)).

jax(_Flags, op(pm), M)
 => M = "\\pm".

ml(_Flags, op(times), M)
 => M = mo(&(times)).

jax(_Flags, op(times), M)
 => M = "\\times".

ml(_Flags, op(sum), M)
 => M = mo(&(sum)).

jax(_Flags, op(sum), M)
 => M = "\\sum".

ml(_Flags, op(prod), M)
 => M = mo(&(prod)).

jax(_Flags, op(prod), M)
 => M = "\\prod".

ml(_Flags, op('#58'), M)
 => M = mo(&('#58')).

jax(_Flags, op('#58'), M)
 => M = ":".

ml(_Flags, op(','), M)
 => M = mo(',').

jax(_Flags, op(','), M)
 => M = ",".

ml(_Flags, op('CircleTimes'), M)
 => M = mo(&('CircleTimes')).

jax(_Flags, op('CircleTimes'), M)
 => M = "\\otimes".

ml(_Flags, op('#x2062'), M)
 => M = mo(&('#x2062')).

jax(_Flags, op('#x2062'), M)
 => M = "{}".

ml(_Flags, op('Tilde'), M)
 => M = mo(&('Tilde')).

jax(_Flags, op('Tilde'), M)
 => M = "\\sim".

ml(_Flags, op(leftrightarrow), M)
 => M = mo(&(leftrightarrow)).

jax(_Flags, op(leftrightarrow), M)
 => M = "\\leftrightarrow".

ml(_Flags, op(iff), M)
 => M = mo(&(iff)).

jax(_Flags, op(iff), M)
 => M = "\\iff".

ml(_Flags, op(rightarrow), M)
 => M = mo(&(rightarrow)).

jax(_Flags, op(rightarrow), M)
 => M = "\\rightarrow".

ml(_Flags, op(rArr), M)
 => M = mo(&(rArr)).

jax(_Flags, op(rArr), M)
 => M = "\\Rightarrow".

ml(_Flags, op(leftarrow), M)
 => M = mo(&(leftarrow)).

jax(_Flags, op(leftarrow), M)
 => M = "\\leftarrow".

ml(_Flags, op(lArr), M)
 => M = mo(&(lArr)).

jax(_Flags, op(lArr), M)
 => M = "\\Leftarrow".

ml(_Flags, op(uparrow), M)
 => M = mo(&(uparrow)).

jax(_Flags, op(uparrow), M)
 => M = "\\uparrow".

ml(_Flags, op(uArr), M)
 => M = mo(&(uArr)).

jax(_Flags, op(uArr), M)
 => M = "\\Uparrow".

ml(_Flags, op(downarrow), M)
 => M = mo(&(downarrow)).

jax(_Flags, op(downarrow), M)
 => M = "\\downarrow".

ml(_Flags, op(dArr), M)
 => M = mo(&(dArr)).

jax(_Flags, op(dArr), M)
 => M = "\\Downarrow".

ml(_Flags, op(approx), M)
 => M = mo(&(approx)).

jax(_Flags, op(approx), M)
 => M = "\\approx".

ml(_Flags, op(equiv), M)
 => M = mo(&(equiv)).

jax(_Flags, op(equiv), M)
 => M = "\\equiv".

ml(_Flags, op(cong), M)
 => M = mo(&(cong)).

jax(_Flags, op(cong), M)
 => M = "\\cong".

ml(_Flags, op(propto), M)
 => M = mo(&(prop)).

jax(_Flags, op(propto), M)
 => M = "\\propto".

ml(_Flags, op(and), M)
 => M = mo(&(and)).

jax(_Flags, op(and), M)
 => M = "\\land".

ml(_Flags, op(or), M)
 => M = mo(&(or)).

jax(_Flags, op(or), M)
 => M = "\\lor".

ml(_Flags, op(not), M)
 => M = mo(&(not)).

jax(_Flags, op(not), M)
 => M = "\\lnot".

ml(_Flags, op(veebar), M)
 => M = mo(&(veebar)).

jax(_Flags, op(veebar), M)
 => M = "\\veebar".

ml(_Flags, op(isin), M)
 => M = mo(&(isin)).

jax(_Flags, op(isin), M)
 => M = "\\in".

ml(_Flags, op(notin), M)
 => M = mo(&(notin)).

jax(_Flags, op(notin), M)
 => M = "\\notin".

ml(_Flags, op(cap), M)
 => M = mo(&(cap)).

jax(_Flags, op(cap), M)
 => M = "\\cap".

ml(_Flags, op(cup), M)
 => M = mo(&(cup)).

jax(_Flags, op(cup), M)
 => M = "\\cup".

ml(_Flags, op(A), M)
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

ml(_Flags, posint(A), M)
 => M = mn(A).

ml(_Flags, pos(1.0Inf), M)
 => M = mi(&('#x221E')).

ml(Flags, pos(A), M)
 => option(round(D), Flags, 2),
    format(atom(Mask), '~~~wf', [D]),
    format(string(X), Mask, [A]),
    M = mn(X).

jax(_Flags, posint(A), M)
 => format(string(M), "{~w}", [A]).

jax(_Flags, pos(1.0Inf), M)
 => M = "\\infty".

jax(Flags, pos(A), M)
 => option(round(D), Flags, 2),
    format(atom(Mask), '{~~~wf}', [D]),
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
ml(Flags, fy(Prec, Op, A), M)
 => ml(Flags, op(Op), S),
    ml(Flags, right(Prec, A), X),
    M = mrow([S, X]).

ml(Flags, yf(Prec, Op, A), M)
 => ml(Flags, op(Op), S),
    ml(Flags, left(Prec, A), X),
    M = mrow([X, S]).

ml(Flags, xfx(Prec, Op, A, B), M)
 => ml(Flags, left(Prec-1, A), X),
    ml(Flags, op(Op), S),
    ml(Flags, right(Prec-1, B), Y),
    M = mrow([X, S, Y]).

ml(Flags, yfx(Prec, Op, A, B), M)
 => ml(Flags, left(Prec, A), X),
    ml(Flags, op(Op), S),
    ml(Flags, right(Prec-1, B), Y),
    M = mrow([X, S, Y]).

ml(Flags, xfy(Prec, Op, A, B), M)
 => ml(Flags, left(Prec-1, A), X),
    ml(Flags, op(Op), S),
    ml(Flags, right(Prec, B), Y),
    M = mrow([X, S, Y]).

ml(Flags, yfy(Prec, Op, A, B), M)
 => ml(Flags, left(Prec, A), X),
    ml(Flags, op(Op), S),
    ml(Flags, right(Prec, B), Y),
    M = mrow([X, S, Y]).

jax(Flags, fy(Prec, Op, A), M)
 => jax(Flags, op(Op), S),
    jax(Flags, right(Prec, A), X),
    format(string(M), "{~w~w}", [S, X]).

jax(Flags, yf(Prec, Op, A), M)
 => jax(Flags, op(Op), S),
    jax(Flags, left(Prec, A), X),
    format(string(M), "{~w~w}", [X, S]).

jax(Flags, xfx(Prec, Op, A, B), M)
 => jax(Flags, left(Prec-1, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec-1, B), Y),
    format(string(M), "{~w~w~w}", [X, S, Y]).

jax(Flags, yfx(Prec, Op, A, B), M)
 => jax(Flags, left(Prec, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec-1, B), Y),
    format(string(M), "{~w~w~w}", [X, S, Y]).

jax(Flags, xfy(Prec, Op, A, B), M)
 => jax(Flags, left(Prec-1, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec, B), Y),
    format(string(M), "{~w~w~w}", [X, S, Y]).

jax(Flags, yfy(Prec, Op, A, B), M)
 => jax(Flags, left(Prec, A), X),
    jax(Flags, op(Op), S),
    jax(Flags, right(Prec, B), Y),
    format(string(M), "{~w~w~w}", [X, S, Y]).

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
 => M = left(Prec+1, A).

denoting(Flags, left(_, A), D)
 => denoting(Flags, A, D).

denoting(Flags, right(_, A), D)
 => denoting(Flags, A, D).

%
% Abbreviations
%
% with s^2_pool denoting the pooled variance
%
ml(Flags, with(A, _, _), X)
 => ml(Flags, A, X).

jax(Flags, with(A, _, _), X)
 => jax(Flags, A, X).

paren(Flags, with(A, _, _), Paren)
 => paren(Flags, A, Paren).

prec(Flags, with(A, _, _), Prec)
 => prec(Flags, A, Prec).

type(Flags, with(A, _, _), Type)
 => type(Flags, A, Type).

denoting(Flags, with(A, Expr, Info), Den)
 => denoting(Flags, Expr, T),
    Den = [denoting(A, Expr, Info) | T].

%
% Expand abbreviations
%
ml(Flags, denoting(A, Expr, Info), X)
 => ml(Flags, list(space, [A = Expr, "denoting", Info]), X).

type(Flags, denoting(A, _, _), Type)
 => type(Flags, A, Type).

denoting(_Flags, denoting(_, _, _), Den)
 => Den = [].

%
% Collect abbreviations
%
ml(Flags, with(Abbreviations), X)
 => sort(Abbreviations, Sorted), % remove duplicates
    ml(Flags, with_(Sorted), X).

ml(_Flags, with_([]), W)
 => W = "".

ml(Flags, with_([A]), W)
 => ml(Flags, A, X),
    W = span([", with", &(nbsp), math(X)]).

ml(Flags, with_([A, B | T]), W)
 => ml(Flags, A, X),
    ml(Flags, and([B | T]), Y),
    W = span([", with", &(nbsp), math(X) | Y]).

ml(_Flags, and([]), W)
 => W = ".".

ml(Flags, and([A | T]), W)
 => ml(Flags, A, X),
    ml(Flags, and(T), Y),
    W = span([", and", &(nbsp), math(X) | Y]).

jax(Flags, with(Abbreviations), X)
 => sort(Abbreviations, Sorted), % remove duplicates
    jax(Flags, with_(Sorted), X).

jax(_Flags, with_([]), W)
 => W = "".

jax(Flags, with_([A]), W)
 => jax(Flags, A, X),
    format(string(W), ", with~~$~w$", [X]).

jax(Flags, with_([A, B | T]), W)
 => jax(Flags, A, X),
    jax(Flags, and([B | T]), Y),
    format(string(W), ", with~~$~w$~w", [X, Y]).

jax(_Flags, and([]), W)
 => W = ".".

jax(Flags, and([A | T]), W)
 => jax(Flags, A, X),
    jax(Flags, and(T), Y),
    format(string(W), ", and~~$~w$~w", [X, Y]).

%
% Parentheses
%
math('('(A), M)
 => M = paren(A).

ml(Flags, paren(A), M),
    paren(Flags, A, P),
    2 is P mod 3
 => ml(Flags, braces(A), M).

ml(Flags, paren(A), M),
    paren(Flags, A, P),
    1 is P mod 3
 => ml(Flags, brackets(A), M).

ml(Flags, paren(A), M)
 => ml(Flags, parentheses(A), M).

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

ml(Flags, parentheses(A), M)
 => ml(Flags, A, X),
    M = mrow([mo('('), X, mo(')')]).

jax(Flags, parentheses(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\left(~w\\right)", [X]).

paren(_Flags, parentheses(_), P)
 => P = 1.

type(_Flags, parentheses(_), T)
 => T = paren.

ml(Flags, brackets(A), M)
 => ml(Flags, A, X),
    M = mrow([mo('['), X, mo(']')]).

jax(Flags, brackets(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\left[~w\\right]", [X]).

paren(_Flags, brackets(_), P)
 => P = 2.

type(_Flags, brackets(_), T)
 => T = paren.

ml(Flags, braces(A), M)
 => ml(Flags, A, X),
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

ml(Flags, list(_, [A]), M)
 => ml(Flags, A, M).

ml(Flags, list(Sep, [A, B | T]), M)
 => ml(Flags, A, X),
    ml(Flags, tail(Sep, [B | T]), Y),
    M = mrow([X | Y]).

ml(Flags, tail(Sep, [A]), M)
 => ml(Flags, Sep, S),
    ml(Flags, A, X),
    M = [S, X].

ml(Flags, tail(Sep, [A, B | T]), M)
 => ml(Flags, Sep, S),
    ml(Flags, A, X),
    ml(Flags, tail(Sep, [B | T]), Y),
    M = [S, X | Y].

jax(Flags, list(_, [A]), M)
 => jax(Flags, A, M).

jax(Flags, list(Sep, [A, B | T]), M)
 => jax(Flags, A, X),
    jax(Flags, tail(Sep, [B | T]), Y),
    format(string(M), "{~w~w}", [X, Y]).

jax(Flags, tail(Sep, [A]), M)
 => jax(Flags, Sep, S),
    jax(Flags, A, X),
    format(string(M), "{~w~w}", [S, X]).

jax(Flags, tail(Sep, [A, B | T]), M)
 => jax(Flags, Sep, S),
    jax(Flags, A, X),
    jax(Flags, tail(Sep, [B | T]), Y),
    format(string(M), "{~w~w~w}", [S, X, Y]).

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
ml(Flags, frac(N, D), M)
 => ml(Flags, N, X),
    ml(Flags, D, Y),
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
ml(Flags, integrate(Fn, From, To, DX), M)
 => ml(Flags, Fn, XFn),
    ml(Flags, From, XFrom),
    ml(Flags, To, XTo),
    ml(Flags, DX, XDX),
    ml(Flags, space, Space),
    M = mrow([munderover([mo(&(int)), XFrom, XTo]), XFn, Space, mi(d), XDX]).

paren(Flags, integrate(_, _, _, A), Paren)
 => paren(Flags, A, Paren).

prec(_Flags, integrate(_, _, _, _), Prec)
 => current(Prec, yfx, +).

%
% Large font ("displaystyle")
%
ml(Flags, display(A), M)
 => ml(Flags, A, X),
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
ml(Flags, overline(A), M)
 => ml(Flags, A, X),
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
ml(Flags, cancel(A), M)
 => ml(Flags, A, X),
    M = menclose(notation(updiagonalstrike), X).

paren(Flags, cancel(A), Paren)
 => paren(Flags, A, Paren).

prec(Flags, cancel(A), Prec)
 => prec(Flags, A, Prec).

type(Flags, cancel(A), Type)
 => type(Flags, A, Type).

%
% Box
%
ml(Flags, box(A), M)
 => ml(Flags, A, X),
    M = menclose(notation(roundedbox), X).

paren(Flags, box(A), Paren)
 => paren(Flags, A, Paren).

prec(Flags, box(A), Prec)
 => prec(Flags, A, Prec).

type(Flags, box(A), Type)
 => type(Flags, A, Type).

%
% Underbrace
%
ml(Flags, underbrace(A, U), M)
 => ml(Flags, A, X),
    ml(Flags, U, Y),
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
    option(error(fix), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [L, box(list(space, [op(Op), R]))]).

math(Flags, omit_right(Expr), New, M),
    option(error(highlight), Flags, highlight),
    Expr =.. [Op, L, R]
 => Flags = New,
    M = list(space, [L, cancel(list(space, [op(Op), R]))]).

math(Flags, instead(_Wrong, Correct), New, M),
    option(error(ignore), Flags, highlight)
 => Flags = New,
    M = Correct.

math(Flags, instead(_Wrong, Correct), New, M),
    option(error(fix), Flags, highlight)
 => Flags = New,
    M = box(Correct).

math(Flags, instead(Wrong, Correct), New, M),
    option(error(highlight), Flags, highlight)
 => Flags = New,
    M = underbrace(list(space, ["instead of", Correct]), Wrong).

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
% Binomial distribution
%
% Density, distribution etc.
math(Flags, dbinom(K, N, Pi), New, X)
 => New = Flags,
    X = fn(subscript('P', "Bi"), (['X' = K] ; [N, Pi])).

math(Flags, pbinom(K, N, Pi), New, X)
 => New = Flags,
    X = fn(subscript('P', "Bi"), (['X' =< K] ; [N, Pi])).

math(Flags, upbinom(K, N, Pi), New, X)
 => New = Flags,
    X = fn(subscript('P', "Bi"), (['X' >= K] ; [N, Pi])).

math(Flags, cbinom(Alpha, N, Pi, Tail, Dist), New, X)
 => New = Flags,
    X = fn(Tail, [fn(subscript('P', "Bi"), ([Dist] ; [N, Pi])) =< Alpha]).

math(Flags, tail("upper"), New, X)
 => New = Flags,
    X = subscript("argmin", k).

math(Flags, tail("lower"), New, X)
 => New = Flags,
    X = subscript("argmax", k).

math(Flags, tail("upperdens"), New, X)
 => New = Flags,
    X = subscript("argmin", k > 'N' * pi).

math(Flags, tail("lowerdens"), New, X)
 => New = Flags,
    X = subscript("argmax", k < 'N' * pi).

math(Flags, dist("upper"), New, X)
 => New = Flags,
    X = ('X' >= k).

math(Flags, dist("lower"), New, X)
 => New = Flags,
    X = ('X' =< k).

math(Flags, dist("density"), New, X)
 => New = Flags,
    X = ('X' = k).

%
% Normal distribution
%
math(pnorm(X, Mu, Sigma2), M)
 => M = fn('Phi', ([X] ; [Mu, Sigma2])).

math(pnorm(Z), M)
 => M = fn('Phi', [Z]).

%
% Functions like f(x) and f(x; a, b)
%
ml(Flags, fn(Name, (Args ; Params)), M)
 => ml(Flags, Name, F),
    ml(Flags, paren(list(op(';'), [list(op(','), Args), list(op(','), Params)])), X),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, (Args ; Params)), M),
    string(Name)
 => jax(Flags, Name, F),
    jax(Flags, paren(list(op(';'), [list(op(','), Args), list(op(','), Params)])), X),
    format(string(M), "{~w\\,~w}", [F, X]).

jax(Flags, fn(Name, (Args ; Params)), M)
 => jax(Flags, Name, F),
    jax(Flags, paren(list(op(';'), [list(op(','), Args), list(op(','), Params)])), X),
    format(string(M), "{~w~w}", [F, X]).

paren(Flags, fn(_Name, (Args ; Params)), Paren)
 => paren(Flags, list(op(','), Args), X),
    paren(Flags, list(op(','), Params), Y),
    Paren is max(X, Y) + 1.

prec(Flags, fn(_Name, (_Args ; _Params)), Prec)
 => prec(Flags, a * b, P0),
    Prec is P0 - 1.

type(_Flags, fn(_Name, (_Args ; _Params)), Type)
 => Type = paren.

ml(Flags, fn(Name, [Arg]), M),
    type(Flags, Arg, paren)
 => ml(Flags, Name, F),
    ml(Flags, Arg, X),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, [Arg]), M),
    type(Flags, Arg, paren)
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w~w}", [F, X]).

% Omit parenthesis in special functions
ml(Flags, fn(Name, [Arg]), M),
%    type(Flags, Arg, function),
    type(Flags, Name, special),
    prec(Flags, Name, P),
    prec(Flags, Arg, Prec),
%    writeln(P-Prec),
    P > Prec
 => ml(Flags, Name, F),
    ml(Flags, Arg, X),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, [Arg]), M),
%    type(Flags, Arg, function)
    type(Flags, Name, special),
    prec(Flags, Name, P),
    prec(Flags, Arg, Prec),
%    writeln(P-Prec),
    P > Prec
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w~w}", [F, X]).

ml(Flags, fn(Name, [Arg]), M),
    type(Flags, Name, special),
    prec(Flags, Arg, 0)
 => ml(Flags, Name, F),
    ml(Flags, Arg, X),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, [Arg]), M),
    type(Flags, Name, special),
    prec(Flags, Arg, 0)
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w~w}", [F, X]).

ml(Flags, fn(Name, Args), M)
 => ml(Flags, Name, F),
    ml(Flags, paren(list(op(','), Args)), X),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, Args), M)
 => jax(Flags, Name, F),
    jax(Flags, paren(list(op(','), Args)), X),
    format(string(M), "{~w~w}", [F, X]).

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
ml(Flags, A, M),
    compound(A),
    compound_name_arguments(A, N, Args)
 => ml(Flags, fn(N, Args), M).

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

