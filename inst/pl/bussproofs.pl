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

% Render in MathML
mlx(rcond(A, B), M, Flags) :-
    build_table([B], [##(A)], M1),
    M2 = [##('') | M1],
    M3 =.. [### | M2],
    ml(proof_tree(M3), M4, Flags),
    M = mrow([M4]).
 
build_table([], Acc, Acc).

build_table([B], Acc, Matrix) :-
    B =.. L,
    nth0(0, L, Func, L1),
    nth0(0, L1, Arg, L2),
    B1 =.. [Func, Arg],
    append([##(B1)], Acc, List),
    build_table(L2, List, Matrix). 

/* 
mlx(rcond(A, B), M, Flags) :-
    ml(frac(B, A), F1, Flags),
    ml('%->%'('R', ''), R1, Flags),
    M = mrow([F1, mstyle([mathsize('0.7em')], R1)]).

mlx(land(A, B), M, Flags) :-
    ml(frac(B, A), F1, Flags),
    ml('&'('L', ''), R1, Flags),
    M = mrow([F1, mstyle([mathsize('0.7em')], R1)]).

mlx(ax(A, B), M, Flags) :-
    ml(frac(B, A), F1, Flags),
    ml(ident('Ax.'), R1, [mathvariant(italic)]),
    M = mrow([F1, mstyle([mathsize('0.7em')], R1)]). */

% Render in MathJax
% jaxx(rcond(A, B), M, Flags) :-
%     jax(A, A1, Flags),
%     jax(B, B1, Flags),
%     format(string(M), "\\frac{~w}{~w}", [B1, A1]).

