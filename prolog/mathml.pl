:- module(mathml, [pl_mathml/2, pl_mathml/3, pl_mathjax/2, pl_mathjax/3]).

:- discontiguous mathml/0, math/2, math/3, math/4, current/3, paren/3, prec/3.
:- discontiguous type/3, denoting/3, ml/3, jax/3.

:- use_module(library(http/html_write)).
:- use_module(library(rolog)).
:- consult(['../inst/pl/lib/core.pl', '../inst/pl/lib/op.pl']).

% Hook to defined own macros
%
% Example
% assert(math_hook(t0, subscript(t, 0))).
%
% From R, the hook is installed by
% mathml:hook(t0, subscript(t, 0))
%
:- dynamic math_hook/2, math_hook/3, math_hook/4.
:- multifile math_hook/2, math_hook/3, math_hook/4.

% Low-level functions (see, e.g. nthroot.pl)
%
% Example
% see nthroot.pl
%
:- multifile mlx/3.    % translate term to mathml
:- multifile jaxx/3.   % translate to LaTeX
:- multifile precx/3.  % operator precedence
:- multifile parenx/3. % count parentheses
:- multifile typex/3.  % some type information

% Translate prolog expression to MathML string
%
% Example
% pl_mathml(sin(pi/2), M).
%
pl_mathml(R, S)
=> pl_mathml(R, S, []).

% The flags allow for context-dependent translation
%
% Examples
% see vignette of R package mathml
%
pl_mathml(R, S, Flags0)
 => digits_(Flags0, Flags1),
    mathml(R, M, Flags1),
    html(M, H, []),
    maplist(atom_string, H, S).

% R interface: Translate R expression to MathJax string
pl_mathjax(R, S)
 => pl_mathjax(R, S, []).

pl_mathjax(R, S, Flags0)
 => digits_(Flags0, Flags1),
    mathjax(R, S, Flags1).

 % Default digits if not defined in flags
 digits_(Flags0, Flags1),
    option(digits(_), Flags0)
 => Flags1 = Flags0.

 digits_(Flags0, Flags1)
  => Flags1 = [digits(2) | Flags0].
