:- module(leanseq, []).

:- reexport(library(mathml)).

:- multifile jaxx/3.


:- op( 500, fy, ~).     % negation
:- op(1000, xfy, &).    % conjunction
:- op(1100, xfy, '|').  % disjunction
:- op(1110, xfy, =>).   % conditional
:- op(1120, xfy, <=>).  % biconditional
:- op( 500, fy, !).     % universal quantifier:  ![X]:
:- op( 500, fy, ?).     % existential quantifier:  ?[X]:
:- op( 500,xfy, :).


atom_split(X, D, L) :-
    atomic_list_concat(L, D, X).

jaxx(prover(H, J), M, Flags) :-
    jax(latex_root(H, J), M0, Flags),
    atom_split(M, '', M0),
    write(M), nl.

jaxx(latex_root(P, J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    format(string(M), "\\begin{prooftree} {~w} \\end{prooftree}", [M0]).
 
jaxx(latex_proof(ax(S, U), J), M, Flags) :-
    jax(latex_sequent(S, ax(S,U), J), M0, Flags),
    format(string(M), "\\RightLabel{$Ax.$} \\AxiomC{} \\UnaryInfC{~w} ", [M0]).

jaxx(latex_proof(rcond(S, P), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_sequent(S, rcond(S, P), J), M1, Flags),
    format(string(M), "{~w}\\RightLabel{$R\\to$} \\UnaryInfC{~w}", [M0, M1]).

jaxx(latex_proof(land(S, P), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_sequent(S, land(S, P), J), M1, Flags),
    format(string(M), "{~w}\\RightLabel{$L\\land$} \\UnaryInfC{~w}", [M0, M1]).

jaxx(latex_proof(ror(S, P), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_sequent(S, ror(S, P), J), M1, Flags),
    format(string(M), "{~w}\\RightLabel{$R\\lor$} \\UnaryInfC{~w}", [M0, M1]).

jaxx(latex_proof(lneg(S, P), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_sequent(S, lneg(S, P), J), M1, Flags),
    format(string(M), "{~w}\\RightLabel{$L\\neg$} \\UnaryInfC{~w}", [M0, M1]).

jaxx(latex_proof(rneg(S, P), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_sequent(S, rneg(S, P), J), M1, Flags),
    format(string(M), "{~w}\\RightLabel{$R\\neg$} \\UnaryInfC{~w}", [M0, M1]).

jaxx(latex_proof(rand(S, P, Q), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_proof(Q, J), M1, Flags),
    jax(latex_sequent(S, rand(S, P, Q), J), M2, Flags),
    format(string(M), "{~w}{~w} \\RightLabel{$R\\land$} \\BinaryInfC{~w}", [M0, M1, M2]).

jaxx(latex_proof(lor(S, P, Q), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_proof(Q, J), M1, Flags),
    jax(latex_sequent(S, lor(S, P, Q), J), M2, Flags),
    format(string(M), "{~w}{~w} \\RightLabel{$L\\lor$} \\BinaryInfC{~w}", [M0, M1, M2]).

jaxx(latex_proof(lcond(S, P, Q), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_proof(Q, J), M1, Flags),
    jax(latex_sequent(S, lcond(S, P, Q), J), M2, Flags),
    format(string(M), "{~w}{~w} \\RightLabel{$L\\to$} \\BinaryInfC{~w}", [M0, M1, M2]).

jaxx(latex_proof(lbicond(S, P, Q), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_proof(Q, J), M1, Flags),
    jax(latex_sequent(S, lbicond(S, P, Q), J), M2, Flags),
    format(string(M), "{~w}{~w} \\RightLabel{$R\\leftrightarrow$} \\BinaryInfC{~w}", [M0, M1, M2]).

jaxx(latex_proof(rbicond(S, P, Q), J), M, Flags) :-
    jax(latex_proof(P, J), M0, Flags),
    jax(latex_proof(Q, J), M1, Flags),
    jax(latex_sequent(S, rbicond(S, P, Q), J), M2, Flags),
    format(string(M), "{~w}{~w} \\RightLabel{$R\\leftrightarrow$} \\BinaryInfC{~w}", [M0, M1, M2]).

jaxx(latex_proof(asq(S, U), J), M, Flags) :-
    jax(latex_antisequent(S, asq(S, U), J), M0, Flags),
    format(string(M), "\\RightLabel{$R\\Asq.$} \\AxiomC{} \\UnaryInfC{~w}", [M0]).

jaxx(latex_sequent(G > D, P, J), M, Flags) :-
    jax(latex_list(G, P, left, 0, J), M0, Flags),
    jax(latex_list(D, P, right, 0, J), M1, Flags),
    format(string(M), "$ {~w} \\vdash {~w} $", [M0, M1]).

jaxx(latex_antisequent(G > D, P, J), M, Flags) :-
    jax(latex_list(G, P, left, 0, J), M0, Flags),
    jax(latex_list(D, P, right, 0, J), M1, Flags),
    format(string(M), "$ {~w} \\nvdash {~w} $", [M0, M1]).

jaxx(latex_list([], _, _, _, _), M, _Flags) :-
    format(string(M), "").

/* jaxx(latex_list([X|L], P, F, N, J), M, Flags) :-
    has_usage(P, F, N),
    !,
    jax(latex_formula(X, J), M0, Flags),
    I is N+1,
    jax(latex_rest(L, P, F, I, J), M1, Flags),
    format(string(M), "{~w}{~w}", [M0, M1]). */

jaxx(latex_list([X|L], P, F, N, J), M, Flags) :-
    jax(latex_formula(X, J), M0, Flags),
    I is N+1,
    jax(latex_rest(L, P, F, I, J), M1, Flags),
    format(string(M), "{~w}{~w}", [M0, M1]).

jaxx(latex_list([_|L], P, F, N, J), M, Flags) :-
    I is N+1,
    jax(latex_list(L, P, F, I, J), M0, Flags),
    format(string(M), "{~w}", [M0]).

jaxx(latex_rest([], _, _, _, _), M, _Flags) :-
    format(string(M), "").

/* jaxx(latex_rest([X|L], P, F, N, J), M, Flags) :-
    has_usage(P, F, N),
    !,
    jax(latex_formula(X, J), M0, Flags),
    I is N+1, 
    jax(latex_rest(L, P, F, I, J), M1, Flags),
    format(string(M), ", {~w}{~w}", [M0, M1]). */

jaxx(latex_rest([X|L], P, F, N, J), M, Flags) :-
    jax(latex_formula(X, J), M0, Flags),
    I is N+1,
    jax(latex_rest(L, P, F, I, J), M1, Flags),
    format(string(M), ", {~w}{~w}", [M0, M1]).

jaxx(latex_rest([_|L], P, F, N, J), M, Flags) :-
    I is N+1,
    jax(latex_rest(L, P, F, I, J), M0, Flags),
    format(string(M), "{~w}", [M0]).

jaxx(latex_formula(~A, J), M, Flags) :-
    !,
    jax(latex_term(A, J), M0, Flags),
    format(string(M), "\\neg {~w}", [M0]).

jaxx(latex_formula((A&B), J), M, Flags) :-
    !,
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "{~w} \\land {~w}", [M0, M1]).

jaxx(latex_formula((A|B), J), M, Flags) :-
    !,
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "{~w} \\lor {~w}", [M0, M1]).

jaxx(latex_formula((A=>B), J), M, Flags) :-
    !,
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "{~w} \\to {~w}", [M0, M1]).

jaxx(latex_formula((A<=>B), J), M, Flags) :-
    !,
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "{~w} \\leftrightarrow {~w} ", [M0, M1]).

/* jaxx(latex_formula(![X], J), M, Flags) :-
    !,
    jax(latex_term(X, J), M0, Flags),
    format(string(M), "\\forall {~w}", [M0]).

jaxx(latex_formula(?[Y], J), M, Flags) :-
    !,
    jax(latex_term(Y, J), M0, Flags),
    format(string(M), "\\exists {~w}", [M0]). */

jaxx(latex_formula(X, J), M, Flags) :-
    jax(latex_factor(X, J), M, Flags).

jaxx(latex_term(~A, J), M, Flags) :-
    !, 
    jax(latex_term(A, J), M0, Flags),
    format(string(M), "\\neg {~w} ", [M0]).
 
jaxx(latex_term((A&B), J), M, Flags) :-
    !, 
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "( {~w} \\land {~w} )", [M0, M1]).

jaxx(latex_term((A|B), J), M, Flags) :-
    !, 
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "( {~w} \\lor {~w} )", [M0, M1]).

jaxx(latex_term((A=>B), J), M, Flags) :-
    !, 
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "( {~w} \\to {~w} )", [M0, M1]).

jaxx(latex_term((A<=>B), J), M, Flags) :-
    !, 
    jax(latex_term(A, J), M0, Flags),
    jax(latex_term(B, J), M1, Flags),
    format(string(M), "( {~w} \\leftrightarrow {~w} )", [M0, M1]).

/* jaxx(latex_term(![X], J), M, Flags) :-
    !,
    jax(latex_term(X, J), M0, Flags),
    format(string(M), "\\forall {~w}", [M0]).

jaxx(latex_term(?[Y], J), M, Flags) :-
    !,
    jax(latex_term(Y, J), M0, Flags),
    format(string(M), "\\exists {~w}", [M0]). */

jaxx(latex_term(X, J), M, Flags) :-
    jax(latex_factor(X, J), M, Flags).

jaxx(latex_factor(X, J), M, _Flags) :-
    var(X),
    member(W = N, J),
    W == X,
    !,
    format(string(M), "{~w}", [N]).

jaxx(latex_factor(X, _), M, _Flags) :-
    var(X),
    !,
    format(string(M), "?").

jaxx(latex_factor(X, J), M, Flags) :-
    X=..[F,Y|L],
    !,
    jax(latex_factor(Y, J), M0, Flags),
    jax(latex_args(L, J), M1, Flags),
    format(string(M), "{~w} ( {~w} {w} )", [F, M0, M1]).

jaxx(latex_factor(X, _), M, _Flags) :-
    format(string(M), "{~w}", [X]).  

jaxx(latex_args([], _), M, _Flags) :-
    format(string(M), "").  

jaxx(latex_args([X|L], J), M, Flags) :-
    jax(latex_factor(X, J), M0, Flags),
    jax(latex_args(L, J), M1, Flags),
    format(string(M), ", {~w}{~w}", [M0, M1]). 


