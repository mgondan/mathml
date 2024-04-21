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
mlx(proof(Denominator, Numerator), M, Flags) :-
    X1 =.. ['###1', '##'(Numerator)],
    X =.. ['###2', '##'(X1), '##1'(Denominator)],
    ml(proof_tree(X), M, Flags).

mlx(proof(Denominator, Numerator1, Numerator2), M, Flags) :-
    X1 =.. ['###1', '##'(Numerator1, '', Numerator2)],
    X =.. ['###2', '##'(X1), '##1'(Denominator)],
    ml(proof_tree(X), M, Flags).

mlx(format_operator(Op), R, Flags) :-
    ml(Op, R1, Flags),
    ml(size(R1, '0.7em'), R2, Flags),
    ml(mstyle_right(R2), R3, Flags),
    ml(mpadded_right(R3), R, Flags).

mlx(size(A, Size), M, _Flags) :-
    M = mrow([mstyle([mathsize(Size)], A)]).

mlx(mstyle_right(A), M, _Flags) :-
    M = mstyle([displaystyle('false'), scriptlevel('0')], A).

mlx(mpadded_right(A), M, _Flags) :-
    M = mpadded([height('.5em'), width('.5em'), voffset('.9em'), semantics('bspr_prooflabel:right')], A).

mlx(rbicond(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator('%<->%'('R', '')), R, Flags),
    M = mrow([F, R]).  

mlx(rbicond(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator('%<->%'('R', '')), R, Flags),
    M = mrow([F, R]).  

mlx(rcond(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator('%->%'('R', '')), R, Flags),
    M = mrow([F, R]).  

mlx(rcond(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator('%->%'('R', '')), R, Flags),
    M = mrow([F, R]).  

mlx(rand(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator(and('R', '')), R, Flags),
    M = mrow([F, R]). 

mlx(rand(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator(and('R', '')), R, Flags),
    M = mrow([F, R]). 
    
mlx(ror(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator(or('R', '')), R, Flags),
    M = mrow([F, R]). 

mlx(ror(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator(or('R', '')), R, Flags),
    M = mrow([F, R]). 

mlx(rneg(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator(!('R','')), R, Flags),
    M = mrow([F, R]). 

mlx(rneg(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator(!('R','')), R, Flags),
    M = mrow([F, R]). 

mlx(lbicond(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator('%<->%'('L', '')), R, Flags),
    M = mrow([F, R]). 

mlx(lbicond(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator('%<->%'('L', '')), R, Flags),
    M = mrow([F, R]).  

mlx(lcond(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator('%->%'('L', '')), R, Flags),
    M = mrow([F, R]).  

mlx(lcond(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator('%->%'('L', '')), R, Flags),
    M = mrow([F, R]).  

mlx(land(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator('&'('L', '')), R, Flags),
    M = mrow([F, R]).  

mlx(land(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator('&'('L', '')), R, Flags),
    M = mrow([F, R]).  

mlx(lor(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator(or('L', '')), R, Flags),
    M = mrow([F, R]). 

mlx(lor(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator(or('L', '')), R, Flags),
    M = mrow([F, R]). 

mlx(lneg(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator(!('L', '')), R, Flags),
    M = mrow([F, R]). 

mlx(lneg(A, B, C), M, Flags) :-
    ml(proof(A, B, C), F, Flags),
    ml(format_operator(!('L', '')), R, Flags),
    M = mrow([F, R]). 

mlx(ax(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator(ident('Ax.')), R, [mathvariant(italic)]),
    M = mrow([F, R]).

mlx(asq(A, B), M, Flags) :-
    ml(proof(A, B), F, Flags),
    ml(format_operator(ident('Asq.')), R, [mathvariant(italic)]),
    M = mrow([mathcolor(red)], [F, R]).

% Render in MathJax
% jaxx(rcond(A, B), M, Flags) :-
%     jax(A, A1, Flags),
%     jax(B, B1, Flags),
%     format(string(M), "\\frac{~w}{~w}", [B1, A1]).

