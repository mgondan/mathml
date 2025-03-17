:- module(pval, []).
:- reexport(library(mathml)).

:- multifile math_hook/2, math_hook/3, math_hook/4.

% p-value: if numeric, show e.g. as < 0.001
mathml:math_hook(pval(A), M, Flags, Flags1) :-
    M = A,
    Flags1 = [pval(.) | Flags].

% p-value with p symbol
mathml:math_hook(pval(A, P), M, Flags, Flags1) :-
    M = A,
    Flags1 = [pval(P) | Flags].

% Round t-statistic to two digits
mathml:math_hook(tstat(A), M, Flags, Flags1) :-
    M = A,
    Flags1 = [digits(2) | Flags].

% Render 0.05 as 5%
mathml:math_hook(percent(A), M, Flags, Flags1) :-
    option(digits(D), Flags, 2),
    D1 is D - 2,
    Flags1 = [digits(D1), mult(100) | Flags],
    M = list("", [A, '%']).

