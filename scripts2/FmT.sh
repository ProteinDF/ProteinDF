#!/usr/bin/octave -qf

#printf ("execute>>>> %s", program_name);
#for i = 1:nargin
#  printf (" %s", argv{i});
#endfor
#printf ("\n");

global m = str2num(argv{1})
global T = str2num(argv{2})

function y = f(t)
global m;
global T;
y = (t^(2.0 * m)) * exp(-T * t * t);
endfunction

[answer, iner, nfun, err] = quad("f", 0, 1)


