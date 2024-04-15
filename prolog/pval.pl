:- module(pval, []).
:- reexport(library(mathml)).

:- multifile math_hook/2, math_hook/3, math_hook/4.

% p-value
mathml:math_hook(pval(A), M, Flags, Flags1) :-
    type(A, T, Flags),
    member(numeric(N), T),
    N =< 1,
    N >= 0.1,
    !,
    M = A,
    Flags1 = [round(2) | Flags].

mathml:math_hook(pval(A), M, Flags, Flags1) :-
    type(A, T, Flags),
    member(numeric(_N), T),
    !,
    M = A,
    Flags1 = [round(3) | Flags].

mathml:math_hook(pval(A), M, Flags, Flags1) :-
    !,
    M = A,
    Flags1 = Flags.

mathml:math_hook(pval(A, P), M, Flags) :-
    type(A, T, Flags),
    member(numeric(N), T),
    N < 0.001,
    !,
    M = (P < pval(0.001)).

mathml:math_hook(pval(A, P), M, _Flags) :-
    !,
    M = (P == pval(A)).
