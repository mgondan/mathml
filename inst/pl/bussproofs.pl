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
mlx(table(A, B), M, _Flags) :-
    M1 =.. ['###1', '##'(B)],
    M =.. ['###', '##'(M1), '##1'(A)].

mlx(size(A, S), M, _Flags) :-
    M = mrow([mstyle([mathsize(S)], A)]).

mlx(rcond(A, B), M, Flags) :-
    ml(table(A, B), F1, Flags),
    ml(proof_tree(F1), F2, Flags),
    ml('%->%'('R', ''), R1, Flags),
    ml(size(R1, '0.7em'), R2, Flags),
    M = mrow([F2, mpadded([height('.5em'), width('.5em'), voffset('.9em'), semantics('bspr_prooflabel:right')], mstyle([displaystyle('false'), scriptlevel('0')], R2))]).  

mlx(land(A, B), M, Flags) :-
    ml(table(A, B), F1, Flags),
    ml(proof_tree(F1), F2, Flags),
    ml('&'('L', ''), R1, Flags),
    ml(size(R1, '0.7em'), R2, Flags),
    M = mrow([F2, mpadded([height('.5em'), width('.5em'), voffset('.9em'), semantics('bspr_prooflabel:right')], mstyle([displaystyle('false'), scriptlevel('0')], R2))]).  

mlx(ax(A, B), M, Flags) :-
    ml(table(A, B), F1, Flags),
    ml(proof_tree(F1), F2, Flags),
    ml(ident('Ax.'), R1, [mathvariant(italic)]),
    ml(size(R1, '0.7em'), R2, Flags),
    M = mrow([F2, mpadded([height('.5em'), width('.5em'), voffset('.9em'), semantics('bspr_prooflabel:right')], mstyle([displaystyle('false'), scriptlevel('0')], R2))]). 

% Render in MathJax
% jaxx(rcond(A, B), M, Flags) :-
%     jax(A, A1, Flags),
%     jax(B, B1, Flags),
%     format(string(M), "\\frac{~w}{~w}", [B1, A1]).

