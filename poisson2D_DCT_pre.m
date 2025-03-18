
% Compute the wavenumber matrix needed for the discrete cosine transform

function Mw = poisson2D_DCT_pre(Lx,Ly,nx,ny)

% grid spacing
dx = Lx/(nx-1); 
dy = Ly/(ny-1);

% modified wavenumbers
kx = 0:nx-2;
ky = 0:ny-2;
mwx = 2.*(cos(pi*kx./(nx-1)) - 1)./dx^2;
mwy = 2.*(cos(pi*ky./(ny-1)) - 1)./dy^2;
[Mwx, Mwy] = meshgrid(mwx,mwy);
Mw = Mwx + Mwy;

end