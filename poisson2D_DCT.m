
% Execute the Discrete Cosine Transform of the right-hand-side f, and solve
% the linear system resulting from 2nd order centered differences.

function u = poisson2D_DCT(Mw,f)

% 2D DCT of the right hand side
fhat = dct2(f);
uhat = fhat./Mw;

% The solution is not unique (uhat(0,0) = inf);
uhat(1,1) = 0;

% Inverse 2D DCT
u = idct2(uhat);

end