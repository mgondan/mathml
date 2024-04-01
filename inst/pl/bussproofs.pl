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
mlx(proof(Denominator, Numerator), X, _Flags) :-
    X1 =.. ['###1', '##'(Numerator)],
    X =.. ['###2', '##'(X1), '##1'(Denominator)].

mlx(size(A, Size), M, _Flags) :-
    M = mrow([mstyle([mathsize(Size)], A)]).

mlx(mstyle_right(A), M, _Flags) :-
    M = mstyle([displaystyle('false'), scriptlevel('0')], A).

mlx(mpadded_right(A), M, _Flags) :-
    M = mpadded([height('.5em'), width('.5em'), voffset('.9em'), semantics('bspr_prooflabel:right')], A).

mlx(rcond(A, B), M, Flags) :-
    ml(proof(A, B), F1, Flags),
    ml(proof_tree(F1), F2, Flags),
    ml('%->%'('R', ''), R1, Flags),
    ml(size(R1, '0.7em'), R2, Flags),
    ml(mstyle_right(R2), R3, Flags),
    ml(mpadded_right(R3), R4, Flags),
    M = mrow([F2, R4]).  

mlx(land(A, B), M, Flags) :-
    ml(proof(A, B), F1, Flags),
    ml(proof_tree(F1), F2, Flags),
    ml('&'('L', ''), R1, Flags),
    ml(size(R1, '0.7em'), R2, Flags),
    ml(mstyle_right(R2), R3, Flags),
    ml(mpadded_right(R3), R4, Flags),
    M = mrow([F2, R4]).  

mlx(ax(A, B), M, Flags) :-
    ml(proof(A, B), F1, Flags),
    ml(proof_tree(F1), F2, Flags),
    ml(ident('Ax.'), R1, [mathvariant(italic)]),
    ml(size(R1, '0.7em'), R2, Flags),
    ml(mstyle_right(R2), R3, Flags),
    ml(mpadded_right(R3), R4, Flags),
    M = mrow([F2, R4]).

% Render in MathJax
% jaxx(rcond(A, B), M, Flags) :-
%     jax(A, A1, Flags),
%     jax(B, B1, Flags),
%     format(string(M), "\\frac{~w}{~w}", [B1, A1]).

