% Render in MathML
mlx(nthroot(X, N), M, Flags) :-
    ml(X, X1, Flags),
    ml(N, N1, Flags),
    M = mroot([X1, N1]).

% Render in MathJax
jaxx(nthroot(X, N), M, Flags) :-
    jax(X, X1, Flags),
    jax(N, N1, Flags),
    format(string(M), "\\sqrt[~w]{~w}", [N1, X1]).

% Precedence above power to ensure parenthesis around (nthroot)^2
precx(nthroot(_X, _N), P, _Flags) :-
    current_op(P0, xfy, ^),
    P is P0 + 1.

% Continue counting parentheses below root
parenx(nthroot(X, _N), P, Flags) :-
    paren(X, P, Flags).

% Show x^(1/n) as nthroot(x, n)
math_hook(X ^ '('(1/N), nthroot(X, N)).
