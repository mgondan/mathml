:- discontiguous test/0, math/2, math/3, math/4, current/3, paren/3, prec/3, type/3, denoting/3, ml/3, jax/3.
:- use_module(library(http/html_write)).

%
% R interface
%
r2mathml(A, X) :-
    r2mathml([], A, X).

r2mathml(Flags, A, S) :-
    mathml(Flags, A, M),
    html(M, X, []),
    maplist(atom_string, X, S).

r2mathjax(A, X) :-
    r2mathjax([], A, X).

r2mathjax(Flags, A, X) :-
    mathjax(Flags, A, X).

mathml(Flags, A, X) :-
    ml(Flags, A, M),
    denoting(Flags, A, Denoting),
    ml(Flags, with(Denoting), With),
    !, X = [math(M), With].

mathml(Flags, A, X) :-
    ml(Flags, "Conversion failed", X),
    writeln(A).

mathjax(Flags, A, X) :-
    jax(Flags, A, M),
%    denoting(Flags, A, Denoting),
%    ml(Flags, with(Denoting), With),
    !,
    format(string(X), "$~w$", [M]).

mathjax(Flags, A, X) :-
    jax(Flags, "Conversion failed", X),
    writeln(A).

%
% Macros
%
ml(Flags, A, X),
    math(A, M),
    dif(A, M)
 => ml(Flags, M, X).

ml(Flags, A, X),
    math(Flags, A, M),
    dif(A, M)
 => ml(Flags, M, X).

ml(Flags, A, X),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => ml(New, M, X).

jax(Flags, A, X),
    math(A, M),
    dif(A, M)
 => jax(Flags, M, X).

jax(Flags, A, X),
    math(Flags, A, M),
    dif(A, M)
 => jax(Flags, M, X).

jax(Flags, A, X),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => jax(New, M, X).

prec(Flags, A, P),
    math(A, M),
    dif(A, M)
 => prec(Flags, M, P).

prec(Flags, A, P),
    math(Flags, A, M),
    dif(A, M)
 => prec(Flags, M, P).

prec(Flags, A, P),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => prec(New, M, P).

math(Flags, A, New, X),
    member(replace(A, _), Flags)
 => select(replace(A, X), Flags, New).

%
% Upright text
%
math(A, M),
    string(A)
 => M = text(A).

ml(_Flags, text(A), M)
 => M = mtext(A).

jax(_Flags, text(A), M)
 => format(string(M), "\\mathrm{~w}", [A]).

type(_Flags, text(_), T)
 => T = atomic.

%
% Greek letters
%
math(A, M),
    atom(A),
    memberchk(A, [alpha, beta, gamma, delta, epsilon, varepsilon, zeta, eta,
        theta, vartheta, iota, kappa, lambda, mu, nu, xi, pi, rho, sigma, tau,
        upsilon, phi, varphi, chi, psi, omega, 'Gamma', 'Delta', 'Theta',
        'Lambda', 'Xi', 'Pi', 'Sigma', 'Upsilon', 'Phi', 'Psi', 'Omega'])
 => M = greek(A).

ml(_Flags, greek(A), M) =>
    M = mi(&(A)).

jax(_Flags, greek(A), M) =>
    format(string(M), "\\~w", [A]).

type(_Flags, greek(_), T)
 => T = atomic.

%
% Booleans
%
math('TRUE', M)
 => M = boolean('T').

math('FALSE', M)
 => M = boolean('F').

ml(_Flags, boolean(A), M)
 => M = mi(A).

jax(_Flags, boolean(A), M)
 => format(string(M), "{~w}", [A]).

type(_Flags, boolean(_), T)
 => T = atomic.

%
% Set
%
math('is.null'(A), M)
 => M = (A = 'NULL').

math('NULL', M)
 => M = set(empty).

ml(_Flags, set(empty), M)
 => M = mi(&(empty)).

jax(_Flags, set(empty), M)
 => M = "\\emptyset".

type(_Flags, set(_), T)
 => T = atomic.

%
% Function names
%
math(A, M),
    atom(A),
    member(A, [sgn, sin, cos, tan, asin, acos, atan, arctan2, sinh, cosh, tanh,
               arsinh, arcosh, artanh, log, exp])
 => M = func(A).

ml(_Flags, func(A), X)
 => X = mi(A).

jax(_Flags, func(sgn), M)
 => M = "{\\mathrm{sgn}\\,}".

jax(_Flags, func(A), M)
 => format(string(M), "\\~w", [A]).

%
% Space
%
math(space, X)
 => X = space(thinmathspace).

ml(_Flags, space(Width), M)
 => M = mspace(width(Width), []).

jax(_Flags, space(thinmathspace), M)
 => M = "\\,".

jax(_Flags, space(_Width), M)
 => M = "\\ ".

%
% Symbols/Identifiers
%
math(A, M),
    atom(A)
 => M = ident(A).

ml(Flags, ident(A), M),
    member(mathvariant(calligraphy), Flags)
 => M = mi(mathvariant(script), A).

ml(_Flags, ident(A), M)
 => M = mi(A).

jax(Flags, ident(A), M),
    member(mathvariant(calligraphy), Flags)
 => format(string(M), "\\mathcal{~w}", [A]).

jax(_Flags, ident(A), M)
 => format(string(M), "{~w}", [A]).

type(_Flags, ident(_), T)
 => T = atomic.

%
% R package base
%
math(length(A), M)
 => M = abs(A).

ml(Flags, abs(A), M)
 => ml(Flags, A, X),
    M = mrow([mo(&(vert)), X, mo(&(vert))]).

jax(Flags, abs(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\left|{~w}\\right|", [X]).

paren(_Flags, abs(_), P)
 => P = 0.

prec(Flags, abs(A), P)
 => prec(Flags, paren(A), P).

math(sign(A), M)
 => M = fn(sgn, [A]).

ml(Flags, sqrt(A), M)
 => ml(Flags, A, X),
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
 => M = fn(sup(sin, -1), [A]).

math(acos(A), M)
 => M = fn(sup(cos, -1), [A]).

math(atan(A), M)
 => M = fn(sup(tan, -1), [A]).

math(atan2(y=A, x=B), M)
 => M = atan2(A, B).

math(atan2(A, B), M)
 => M = fn(sup(tan, -1), [A, B]).

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
 => M = fn(sup(sinh, -1), [A]).

math(acosh(A), M)
 => M = fn(sup(cosh, -1), [A]).

math(atanh(A), M)
 => M = fn(sup(tanh, -1), [A]).

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
 => M = fn(sub('I', Nu), [paren(X)]).

% todo: expon.scaled
math(besselK(x=X, nu=Nu), M)
 => M = besselK(X, Nu).

math(besselK(X, Nu), M)
 => M = fn(sub('K', Nu), [paren(X)]).

math(besselJ(x=X, nu=Nu), M)
 => M = besselJ(X, Nu).

math(besselJ(X, Nu), M)
 => M = fn(sub('J', Nu), [paren(X)]).

math(besselY(x=X, nu=Nu), M)
 => M = besselY(X, Nu).

math(besselY(X, Nu), M)
 => M = fn(sub('Y', Nu), [paren(X)]).

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

math(lchoose(n=N, k=K), M)
 => M = lchoose(N, K).

math(lchoose(N, K), M)
 => M = log(choose(N, K)).

math(factorial(x=N), M)
 => M = factorial(N).

math(factorial(N), M)
 => current(Prec, xfy, ^),
    M = yf(Prec, !, N).

math(lfactorial(x=N), M)
 => M = lfactorial(N).

math(lfactorial(N), M)
 => M = log(factorial(N)).

math(and(A, B), M)
 => current(Prec, xfy, ','),
    M = xfy(Prec, &(and), A, B).

math(or(A, B), M)
 => current(Prec, xfy, ';'),
    M = xfy(Prec, &(or), A, B).

math(!(A), M)
 => current(Prec, xfy, ^),
    M = fy(Prec, &(not), A).

math(xor(x=A, y=B), M)
 => M = xor(A, B).

math(xor(A, B), M)
 => current(Prec, xfy, ';'),
    M = xfy(Prec, &(veebar), A, B).

math(exp(A), M)
 => M = fn(exp, [A]).

math(expm1(A), M)
 => M = exp(A) - 1.

math(log(X), M)
 => M = fn(log, [X]).

math(log10(X), M)
 => M = logb(x=X, base=10).

math(log2(X), M)
 => M = logb(x=X, base=2).

% unclear why the cases need to be handled separately
math(logb(x=X, base=B), M)
 => M = logb(X, B).

math(logb(X, B), M)
 => M = fn(sub(log, B), [X]).

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

math((_F :- Body), M)
 => M = Body.

math(function(na, Body, _), M),
    compound_name_arguments(Body, '{', [Arg])
 => M = Arg.

math(function(na, Body, _), M),
    compound_name_arguments(Body, '{', Args)
 => M = Args.

math(function(na, Body, _), M)
 => M = Body.

math(Identical, M),
    compound(Identical),
    compound_name_arguments(Identical, identical, Args)
 => member(x=X, Args),
    member(y=Y, Args),
    M = (X == Y).

ml(Flags, ifelse(test=T, yes=Y, no=N), M)
 => ml(Flags, T, Test),
    ml(Flags, Y, Yes),
    ml(Flags, N, No),
    ml(Flags, space, S),
    M = mrow([mo('{'),
      mtable(columnalign(left),
      [ mtr([Yes, mrow([mtext("if"), S, Test])]),
        mtr([No, mtext("otherwise")])
      ])]).

jax(Flags, ifelse(test=T, yes=Y, no=N), M)
 => jax(Flags, T, Test),
    jax(Flags, Y, Yes),
    jax(Flags, N, No),
    format(string(M),
      "\\left\\{\\begin{array}{ll}~w & \\mathrm{if}\\ ~w\\\\ ~w & \\mathrm{otherwise}\\end{array}\\right.",
      [Yes, Test, No]).

paren(_Flags, ifelse(_), P)
 => P is 0.

math('%in%'(x=X, 'table'=Y), M)
 => M = isin(X, Y).

math('%noin%'(x=X, 'table'=Y), M)
 => M = notin(X, Y).

math(intersect(x=X, y=Y), M)
 => M = intersect(X, Y).

math(union(x=X, y=Y), M)
 => M = union(X, Y).

math(setdiff(x=X, y=Y), M)
 => M = setdiff(X, Y).

math(setdiff(X, Y), M)
 => M = X - Y.

math('%x%'('X'=X, 'Y'=Y), M)
 => M = kronecker(X, Y).

math('&'(A, B), M)
 => M = and(A, B).

math('|'(A, B), M)
 => M = or(A, B).

math(xor(x=A, y=B), M)
 => M = xor(A, B).

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

math(mean(x=A), M)
 => M = mean(A).

math(mean(A), M)
 => M = overline(A).

math(Min, M),
    compound(Min),
    compound_name_arguments(Min, min, Args)
 => M = fn("min", Args).

math(Max, M),
    compound(Max),
    compound_name_arguments(Max, max, Args)
 => M = fn("max", Args).

math(t(x=A), M)
 => M = t(A).

math(t(A), M)
 => M = A^"T".

math(Which, M),
    compound(Which),
    compound_name_arguments(Which, which, Args)
 => M = sub("I", Args).

math('which.max'(x=A), M)
 => M = 'which.max'(A).

math('which.max'(A), M)
 => M = fn("argmax", [A]).

math('which.min'(x=A), M)
 => M = 'which.min'(A).

math('which.min'(A), M)
 => M = fn("argmin", [A]).

% calligraphic letters
math(cal(x=A), M)
 => M = cal(A).

math(Flags, cal(A), New, M)
 => New = [mathvariant(calligraphy) | Flags],
    M = A.

% Sum over index
math(Flags, over(index=I, from=F, to=T, fun=Fn), New, M)
 => New = [index(I), from(F), to(T) | Flags],
    M = Fn.

math(Flags, Sum, New, M),
    compound(Sum),
    compound_name_arguments(Sum, sum, Args),
    select(index(I), Flags, F0),
    select(from(F), F0, F1),
    select(to(T), F1, F2)
 => New=F2,
    M = fn(subsup(op(&(sum)), I=F, T), Args).

ml(Flags, Sum, M),
    compound(Sum),
    compound_name_arguments(Sum, sum, Args)
 => maplist(ml(Flags), Args, X),
    M = mrow([mo(&(sum)), mrow(X)]).

jax(Flags, sum(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\sum{~w}", [X]).

jax(Flags, Sum, M),
    compound(Sum),
    compound_name_arguments(Sum, sum, Args)
 => maplist(jax(Flags), Args, X),
    format(string(M), "\\sum{~w}", [X]).

paren(Flags, Sum, P),
    compound(Sum),
    compound_name_arguments(Sum, sum, Args)
 => maplist(paren(Flags), Args, PX),
    max_list(PX, P).

prec(_Flags, Sum, P),
    compound(Sum),
    compound_name_arity(Sum, sum, _)
 => current(P, yfx, +).

%
%
ml(Flags, sum(I, From, To, A), M)
 => ml(Flags, I = From, XFrom),
    ml(Flags, To, XTo),
    ml(Flags, A, X),
    M = mrow([munderover([mo(&(sum)), XFrom, XTo]), X]).

paren(Flags, sum(_, _, _, A), Paren)
 => paren(Flags, A, Paren).

prec(_Flags, sum(_, _, _, _), Prec)
 => current(Prec, yfx, +).

test :- test(sum(i, 1, 10, i)).

%
% Integrate over range
%
% extract value
math(Flags, $(integrate(Fn, Lower, Upper), value), New, M)
 => Flags = New,
    M = integrate(Fn, Lower, Upper).

% with named arguments
math(Flags, integrate(f=Fn, lower=Lower, upper=Upper), New, M)
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

test :- test(integrate(sin, 0, pi)).

% hats
ml(Flags, hat(A), M)
 => ml(Flags, A, X),
	M = mover(accent(true), [X, mo(&('Hat'))]).

paren(Flags, hat(A), Paren)
 => paren(Flags, A, Paren).

prec(Flags, hat(A), Prec)
 => prec(Flags, A, Prec).

type(Flags, hat(A), Type)
 => type(Flags, A, Type).

test :- test(hat('K')).

test :- test(hat('K'^2)).
test :- test(hat('K')^2).
test :- test(hat('sigma')^2).

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
ml(_Flags, op(A), M)
 => M = mo(A).

jax(_Flags, op(&(le)), M)
 => M = "\\le".

jax(_Flags, op(&(ge)), M)
 => M = "\\ge".

jax(_Flags, op(&(ne)), M)
 => M = "\\ne".

jax(_Flags, op(&(sdot)), M)
 => M = "\\cdot".

jax(_Flags, op(&(times)), M)
 => M = "\\times".

jax(_Flags, op(&(sum)), M)
 => M = "\\sum".

jax(_Flags, op(&(prod)), M)
 => M = "\\prod".

jax(_Flags, op(&('#58')), M)
 => M = ":".

jax(_Flags, op(','), M)
 => M = ",".

jax(_Flags, op(&('CircleTimes')), M)
 => M = "\\otimes".

jax(_Flags, op(&('#x2062')), M)
 => M = "{}".

jax(_Flags, op(&('Tilde')), M)
 => M = "\\sim".

jax(_Flags, op(&(and)), M)
 => M = "\\land".

jax(_Flags, op(&(or)), M)
 => M = "\\lor".

jax(_Flags, op(&(not)), M)
 => M = "\\lnot".

jax(_Flags, op(&(veebar)), M)
 => M = "\\veebar".

jax(_Flags, op(&(isin)), M)
 => M = "\\in".

jax(_Flags, op(&(notin)), M)
 => M = "\\notin".

jax(_Flags, op(&(cap)), M)
 => M = "\\cap".

jax(_Flags, op(&(cup)), M)
 => M = "\\cup".

jax(_Flags, op(A), M)
 => format(string(M), "~w", [A]).

prec(_Flags, op(A), P),
    current(P0, _Fix, A)
 => P = P0.

current(0, fy, &(sum)).

current(Prec, yfx, &(sdot)) :-
    current_op(Prec, yfx, *).

current(Prec, yfy, '<=') :-
    current_op(Prec, xfx, =<).

denoting(_Flags, op(_), D)
 => D = [].

%
% Indices like s_D
%
math(Flags, '['(A, Idx), New, X)
 => math(Flags, sub(A, Idx), New, X).

%
% Check for sub(sup(A, Power), Index)
%
math(Flags, sub(A, Idx), New, M),
    type(Flags, A, sup(Bas, Pwr))
 => New = [replace(sup(Bas, Pwr), subsup(Bas, Idx, Pwr)) | Flags],
    M = A.

%
% Render
%

%math(Flags, sub(A, Idx), M),
%    prec(Flags, sub(A, Idx), Outer),
%    prec(Flags, A, Inner),
%    Outer < Inner
% => M = sub(paren(A), Idx).

ml(Flags, sub(A, B), M)
 => ml(Flags, A, X),
    ml(Flags, B, Y),
    M = msub([X, Y]).

jax(Flags, sub(A, B), M)
 => jax(Flags, A, X),
    jax(Flags, B, Y),
    format(string(M), "{~w_~w}", [X, Y]).

paren(Flags, sub(A, _), P)
 => paren(Flags, A, P).

prec(Flags, sub(A, _), P)
 => prec(Flags, A, P).

type(Flags, sub(A, _), T)
 => type(Flags, A, T).

%
% Powers like s^D
%
% Check for sup(sub(A, Index), Power)
%
math(Flags, sup(A, Pwr), New, M),
    type(Flags, A, sub(Bas, Idx))
 => New = [replace(sub(Bas, Idx), subsup(Bas, Idx, Pwr)) | Flags],
    M = A.

%
% Render
%

%math(Flags, sup(A, Pwr), M),
%    prec(Flags, sup(A, Pwr), Outer),
%    prec(Flags, A, Inner),
%    Outer < Inner
% => M = sup(paren(A), Pwr).

ml(Flags, sup(A, B), M)
 => ml(Flags, A, X),
    ml(Flags, B, Y),
    M = msup([X, Y]).

jax(Flags, sup(A, B), M)
 => jax(Flags, A, X),
    jax(Flags, B, Y),
    format(string(M), "{~w^~w}", [X, Y]).

paren(Flags, sup(A, _), P)
 => paren(Flags, A, P).

prec(_Flags, sup(_, _), P)
 => current(P, xfy, ^).

type(_Flags, sup(A, B), T)
 => T = sup(A, B).

%
% Index and Exponent: s_D^2
%

%math(Flags, subsup(A, Idx, Pwr), X),
%    prec(Flags, subsup(A, Idx, Pwr), Outer),
%    prec(Flags, A, Inner),
%    Outer < Inner
% => X = subsup(paren(A), Idx, Pwr).

ml(Flags, subsup(A, B, C), M)
 => ml(Flags, A, X),
    ml(Flags, B, Y),
    ml(Flags, C, Z),
    M = msubsup([X, Y, Z]).

jax(Flags, subsup(A, B, C), M)
 => jax(Flags, A, X),
    jax(Flags, B, Y),
    jax(Flags, C, Z),
    format(string(M), "{~w_~w^~w}", [X, Y, Z]).

paren(Flags, subsup(A, _, _), P)
 => paren(Flags, A, P).

prec(Flags, subsup(A, _, C), P)
 => prec(Flags, sup(A, C), P).

type(_Flags, subsup(A, B, C), T)
 => T = subsup(A, B, C).

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
    D =< 99,
    format(codes(X), '~99f', [A]),
    nth1(Dot, X, 46),
    N is Dot + D,
    findall(E, (nth1(I, X, E), I =< N), Round),
    string_codes(S, Round),
    M = mn(S).

jax(_Flags, posint(A), M)
 => format(string(M), "{~w}", [A]).

jax(_Flags, pos(1.0Inf), M)
 => M = "\\infty".

jax(Flags, pos(A), M)
 => option(round(D), Flags, 2),
    D =< 99,
    format(codes(X), '~99f', [A]),
    nth1(Dot, X, 46),
    N is Dot + D,
    findall(E, (nth1(I, X, E), I =< N), Round),
    string_codes(M, Round).

math(number(A), M),
    A < 0
 => Abs is abs(A),
    M = -pos(Abs).

math(number(A), M)
 => M = pos(A).

%
% Operators
%
math(Flags, isin(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfx(Prec, &(isin), A, B).

math(Flags, notin(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfx(Prec, &(notin), A, B).

math(Flags, intersect(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, &(cap), A, B).

math(Flags, union(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, &(cup), A, B).

math(Flags, ':'(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, &('#58'), A, B).

math(Flags, kronecker(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfx(Prec, &('CircleTimes'), A, B).

math(Flags, '=='(A, B), New, X)
 => math(Flags, A = B, New, X).

math(Flags, A = B, New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, =, A, B).

math(Flags, A \= B, New, X)
 => New = Flags,
    current_op(Prec, xfx, \=),
    X = xfx(Prec, &(ne), A, B).

math(Flags, A < B, New, X)
 => New = Flags,
    current_op(Prec, xfx, <),
    X = yfy(Prec, <, A, B).

math(_Flags, '<='(A, B), M)
 => M = (A =< B).

math(Flags, A =< B, New, X)
 => New = Flags,
    current_op(Prec, xfx, =<),
    X = yfy(Prec, &(le), A, B).

math(Flags, ~(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, &('Tilde'), A, B).

math(Flags, A > B, New, X)
 => New = Flags,
    current_op(Prec, xfx, >),
    X = yfy(Prec, >, A, B).

math(Flags, A >= B, New, X)
 => New = Flags,
    current_op(Prec, xfx, >=),
    X = yfy(Prec, &(ge), A, B).

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
    type(Flags, B, Type),
    Type = atomic
 => New = Flags,
    X = nodot(A, B).

math(Flags, A * B, New, M)
 => New = Flags,
    M = dot(A, B).

math(_Flags, '%*%'(A, B), M)
 => M = times(A, B).

math(_Flags, crossprod(x=A, y=B), M)
 => M = crossprod(A, B).

math(_Flags, crossprod(A, B), M)
 => M = '%*%'(t(A), B).

math(_Flags, tcrossprod(x=A, y=B), M)
 => M = tcrossprod(A, B).

math(_Flags, tcrossprod(A, B), M)
 => M = '%*%'(A, t(B)).

math(Flags, ~(A, B), New, X)
 => New = Flags,
    current_op(Prec, xfx, =),
    X = yfy(Prec, &('Tilde'), A, B).

test :- test(a * b).
test :- test(a * (b * c)).
test :- test((a * b) * c).
test :- test((2 * b) * c).
test :- test((-2 * b) * c).
test :- test((-2 * 2) * c).
test :- test((-2 * -2) * c).

math(Flags, dot(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfy(Prec, &(sdot), A, B).

math(Flags, nodot(A, B), New, X)
 => New = Flags,
    current_op(Prec, yfx, *),
    X = yfy(Prec, &('#x2062'), A, B).

math(times(A, B), M)
 => current_op(Prec, yfx, *),
    M = yfy(Prec, &(times), A, B).

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
    X = sup(A, B).

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

math(Flags, left(Prec, A), M),
    prec(Flags, A, P),
    P > Prec
 => M = paren(A).

math(left(_, A), M)
 => M = A.

math(right(Prec, A), M)
 => M = left(Prec, A).

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

paren(Flags, with(A, _, _), Paren)
 => paren(Flags, A, Paren).

prec(Flags, with(A, _, _), Prec)
 => prec(Flags, A, Prec).

type(Flags, with(A, _, _), Type)
 => type(Flags, A, Type).

denoting(Flags, with(A, Expr, Info), Den)
 => denoting(Flags, Expr, T),
    Den = [denoting(A, Expr, Info) | T].

test :-
    S2P = with(sub(s, "pool")^2,
                   frac((sub('N', "A") - 1) * sub(s, "A")^2 +
                        (sub('N', "B") - 1) * sub(s, "B")^2,
                        sub('N', "A") + sub('N', "B") - 2),
                   "the pooled variance"),
    test(frac(sub(overline('X'), "A") - sub(overline('X'), "B"),
                  sqrt(nodot(S2P, 1/sub('N', "A") + 1/sub('N', "B"))))).

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
 => W = " ".

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
    P is P0 - 1.

%
% Large fraction
%
math(dfrac(Num, Den), M)
 => M = display(frac(Num, Den)).

%
% Sum over index
%
ml(Flags, sum(I, From, To, A), M)
 => ml(Flags, I = From, XFrom),
    ml(Flags, To, XTo),
    ml(Flags, A, X),
    M = mrow([munderover([mo(&(sum)), XFrom, XTo]), X]).

paren(Flags, sum(_, _, _, A), Paren)
 => paren(Flags, A, Paren).

prec(_Flags, sum(_, _, _, _), Prec)
 => current(Prec, yfx, +).

test :- test(sum(i, 1, 10, i)).

%
% Intgrate over range
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

test :- test(integrate(fn(f, [paren(x)]), 1, 10, x)).

%
% Large font ("displaystyle")
%
ml(Flags, display(A), M)
 => ml(Flags, A, X),
    M = mstyle(displaystyle(true), X).

jax(Flags, display(A), M)
 => jax(Flags, A, X),
    format(string(M), "\\displaystyle{~w}", [X]).

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
 => current(Prec, yfx, *).

type(Flags, overline(A), Type)
 => type(Flags, A, Type).

test :- test(overline('D')).

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

test :- test(cancel('D')).

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

test :- test(box('D')).

%
% Underbrace
%
ml(Flags, underbrace(A, U), M)
 => ml(Flags, A, X),
    ml(Flags, U, Y),
    M = munder([munder(accentunder(true),
                  [Y, mo(stretchy(true), &('UnderBrace'))]), X]).

paren(Flags, underbrace(A, _), Paren)
 => paren(Flags, A, Paren).

prec(Flags, underbrace(A, _), Prec)
 => prec(Flags, A, Prec).

type(Flags, underbrace(A, _), Type)
 => type(Flags, A, Type).

test :- test(underbrace('D', u)).

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

test :- test(dfrac(omit_right(overline('D') - mu),
                   sub(s, 'D') / sqrt('N'))).

test :- test(dfrac(overline('D') - mu,
                   sub(s, 'D') / instead('N', sqrt('N')))).

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
    X = fn(sub('P', "Bi"), (['X' = K] ; [N, Pi])).

math(_Flags, pbinom(q=K, size=N, prob=Pi), M)
 => M = fn(sub('P', "Bi"), (['X' =< K] ; [N, Pi])).

math(Flags, pbinom(K, N, Pi), New, X)
 => New = Flags,
    X = fn(sub('P', "Bi"), (['X' =< K] ; [N, Pi])).

math(Flags, upbinom(K, N, Pi), New, X)
 => New = Flags,
    X = fn(sub('P', "Bi"), (['X' >= K] ; [N, Pi])).

math(Flags, cbinom(Alpha, N, Pi, Tail, Dist), New, X)
 => New = Flags,
    X = fn(Tail, [fn(sub('P', "Bi"), ([Dist] ; [N, Pi])) =< Alpha]).

math(Flags, tail("upper"), New, X)
 => New = Flags,
    X = sub("argmin", k).

math(Flags, tail("lower"), New, X)
 => New = Flags,
    X = sub("argmax", k).

math(Flags, tail("upperdens"), New, X)
 => New = Flags,
    X = sub("argmin", k > 'N' * pi).

math(Flags, tail("lowerdens"), New, X)
 => New = Flags,
    X = sub("argmax", k < 'N' * pi).

math(Flags, dist("upper"), New, X)
 => New = Flags,
    X = ('X' >= k).

math(Flags, dist("lower"), New, X)
 => New = Flags,
    X = ('X' =< k).

math(Flags, dist("density"), New, X)
 => New = Flags,
    X = ('X' = k).

test :- test(dbinom(k, 'N', pi)).
test :- test(pbinom(k, 'N', pi)).
test :- test(upbinom(k, 'N', pi)).
test :- test(cbinom(alpha, 'N', pi, tail("lower"), dist("lower"))).
test :- test(cbinom(alpha, 'N', pi, tail("upper"), dist("upper"))).
test :- test(cbinom(alpha, 'N', pi, tail("lowerdens"), dist("density"))).
test :- test(cbinom(alpha, 'N', pi, tail("upperdens"), dist("density"))).

%
% Functions like f(x) and f(x; a, b)
%
ml(Flags, fn(Name, (Args ; Params)), M)
 => ml(Flags, Name, F),
    ml(Flags, paren(list(op(';'), [list(op(','), Args), list(op(','), Params)])), X),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, (Args ; Params)), M)
 => jax(Flags, Name, F),
    jax(Flags, paren(list(op(';'), [list(op(','), Args), list(op(','), Params)])), X),
    format(string(M), "{~w~w}", [F, X]).

paren(Flags, fn(_Name, (Args ; Params)), Paren)
 => paren(Flags, list(op(','), Args), X),
    paren(Flags, list(op(','), Params), Y),
    Paren is max(X, Y) + 1.

prec(Flags, fn(_Name, (_Args ; _Params)), Prec)
 => prec(Flags, a * b, Prec).

type(_Flags, fn(_Name, (_Args ; _Params)), Type)
 => Type = paren.

ml(Flags, fn(Name, [Arg]), M),
    type(Flags, Arg, paren)
 => ml(Flags, Name, F),
    ml(Flags, Arg, X),
    M = mrow([F, mo(&(af)), X]).

ml(Flags, fn(Name, [Arg]), M),
    type(Flags, Arg, function)
 => ml(Flags, Name, F),
    ml(Flags, Arg, X),
    M = mrow([F, mo(&(af)), X]).

ml(Flags, fn(Name, [Arg]), M),
    prec(Flags, Arg, P),
    P = 0
 => ml(Flags, Name, F),
    ml(Flags, Arg, X),
    M = mrow([F, mo(&(af)), X]).

ml(Flags, fn(Name, Args), M)
 => ml(Flags, Name, F),
    ml(Flags, paren(list(op(','), Args)), X),
    M = mrow([F, mo(&(af)), X]).

jax(Flags, fn(Name, [Arg]), M),
    type(Flags, Arg, paren)
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w~w}", [F, X]).

jax(Flags, fn(Name, [Arg]), M),
    type(Flags, Arg, function)
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w~w}", [F, X]).

jax(Flags, fn(Name, [Arg]), M),
    prec(Flags, Arg, P),
    P = 0
 => jax(Flags, Name, F),
    jax(Flags, Arg, X),
    format(string(M), "{~w~w}", [F, X]).

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
 => paren(Flags, Arg, P1),
    succ(P1, P).

paren(Flags, fn(_Name, Args), P)
 => paren(Flags, list(op(','), Args), P0),
    succ(P0, P).

prec(_Flags, fn(_Name, _Args), Prec)
 => current(Prec, yfx, *).

type(_Flags, fn(_Name, _Args), Type)
 => Type = function.

% Default compounds
%math(_Flags, A, M),
%    compound(A),
%    compound_name_arguments(A, N, Args)
% => M = fn(N, Args).

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

paren(_, _, P) =>
    P = 0.

prec(Flags, A, Den),
    math(Flags, A, New, M),
    dif(Flags-A, New-M)
 => prec(New, M, Den).

prec(_, _, P) =>
    P = 0.

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

%
% Tests
%
test(A) :-
    writeln(A),
    mathml([], A, M),
    html(math(M)).
