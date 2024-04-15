:- module(pval, []).
:- reexport(library(mathml)).

:- multifile math_hook/2, math_hook/3, math_hook/4.

% p-value
math_hook(pval(A), M, Flags, Flags1) :-
    mathml:type(A, T, Flags),
    member(numeric(N), T),
    N =< 1,
    N >= 0.1,
    !,
    M = A,
    Flags1 = [round(2) | Flags].

math_hook(pval(A), M, Flags, Flags1) :-
    mathml:type(A, T, Flags),
    member(numeric(_N), T),
    !,
    M = A,
    Flags1 = [round(3) | Flags].

math_hook(pval(A), M, Flags, Flags1) :-
    !,
    M = A,
    Flags1 = Flags.

math_hook(pval(A, P), M, Flags) :-
    mathml:type(A, T, Flags),
    member(numeric(N), T),
    N < 0.001,
    !,
    M = (P < pval(0.001)).

math_hook(pval(A, P), M, _Flags) :-
    !,
    M = (P == pval(A)).
