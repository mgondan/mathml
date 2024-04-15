% p-value
math_hook(pval(A), M, Flags) :-
    type(A, T, Flags),
    member(numeric(N), T),
    N =< 1,
    N >= 0.1,
    !,
    M = round(A, 2).

math_hook(pval(A), M, Flags) :-
    type(A, T, Flags),
    member(numeric(_N), T),
    !,
    M = round(A, 3).

math_hook(pval(A), M, _Flags) :-
    !,
    M = round(A, 4).

math_hook(pval(A, P), M, Flags) :-
    type(A, T, Flags),
    member(numeric(N), T),
    N < 0.001,
    !,
    M = (P < pval(0.001)).

math_hook(pval(A, P), M, _Flags) :-
    !,
    M = (P == pval(A)).
