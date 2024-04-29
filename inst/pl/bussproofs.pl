:- module(bussproofs,
  [ op(500, fy, ~),     % negation
    op(1000, xfy, &),   % conjunction
    op(1100, xfy, '|'), % disjunction
    op(1110, xfy, =>),  % conditional
    op(1120, xfy, <=>)  % biconditional
  ]).

:- reexport(library(mathml)).

:- multifile mlx/3.

% Render in MathML
%
%      <mpadded height="+.5em" width="+.5em" voffset="-.15em" semantics="bspr_prooflabel:right">
%        <mstyle displaystyle="false" scriptlevel="0">
%          <mrow data-mjx-texclass="ORD">
%            <mstyle mathsize="0.7em">
%              <mrow data-mjx-texclass="ORD">
%                <mi>R</mi>
%                <mo accent="false" stretchy="false">&#x2192;</mo>
%              </mrow>
%            </mstyle>
%          </mrow>
%        </mstyle>
%      </mpadded>
%
% math_hook(rcond(A, B), [frac(B, A), '%->%'('R', xx)]).

mathml:math_hook([], '').
mathml:math_hook(&(A, B), and(A, B)).
mathml:math_hook('|'(A, B), or(A, B)).
mathml:math_hook(=>(A, B), '%->%'(A, B)).
mathml:math_hook(<=>(A, B), '%<->%'(A, B)).

% Render in MathML
mathml:mlx(proof(Denominator, Numerator), M, Flags) :-
    X1 =.. ['###1', '##'(Numerator)],
    X =.. ['###2', '##'(X1), '##1'(Denominator)],
    mathml:ml(proof_tree(X), M, Flags).

mathml:mlx(proof(Denominator, Numerator1, Numerator2), M, Flags) :-
    X1 =.. ['###1', '##'(Numerator1, '', Numerator2)],
    X =.. ['###2', '##'(X1), '##1'(Denominator)],
    mathml:ml(proof_tree(X), M, Flags).

mathml:mlx(format_operator(Op), R, Flags) :-
    mathml:ml(Op, R1, Flags),
    mathml:ml(size(R1, '0.7em'), R2, Flags),
    mathml:ml(mstyle_right(R2), R3, Flags),
    mathml:ml(mpadded_right(R3), R, Flags).

mathml:mlx(size(A, Size), M, _Flags) :-
    M = mrow([mstyle([mathsize(Size)], A)]).

mathml:mlx(mstyle_right(A), M, _Flags) :-
    M = mstyle([displaystyle('false'), scriptlevel('0')], A).

mathml:mlx(mpadded_right(A), M, _Flags) :-
    M = mpadded([height('.5em'), width('.5em'), voffset('.9em'), semantics('bspr_prooflabel:right')], A).

mathml:mlx(rbicond(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator('%<->%'('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(rbicond(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator('%<->%'('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(rcond(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator('%->%'('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(rcond(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator('%->%'('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(rand(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator(and('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(rand(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator(and('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(ror(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator(or('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(ror(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator(or('R', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(rneg(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator(!('R','')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(rneg(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator(!('R','')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lbicond(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator('%<->%'('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lbicond(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator('%<->%'('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lcond(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator('%->%'('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lcond(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator('%->%'('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(land(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator('&'('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(land(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator('&'('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lor(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator(or('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lor(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator(or('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lneg(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator(!('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(lneg(A, B, C), M, Flags) :-
    mathml:ml(proof(A, B, C), F, Flags),
    mathml:ml(format_operator(!('L', '')), R, Flags),
    M = mrow([F, R]).

mathml:mlx(ax(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator(ident('Ax.')), R, [mathvariant(italic)]),
    M = mrow([F, R]).

mathml:mlx(asq(A, B), M, Flags) :-
    mathml:ml(proof(A, B), F, Flags),
    mathml:ml(format_operator(ident('Asq.')), R, [mathvariant(italic)]),
    M = mrow([mathcolor(red)], [F, R]).
