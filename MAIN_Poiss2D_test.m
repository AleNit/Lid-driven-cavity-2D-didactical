
% Solve the 2D Poisson equation with homogeneous Neumann boundary condition
% over a uniform staggered grid with a Discrete Cosine Transform. 
% Derivatives are discretized by 2nd-order accurate centered differences.
% A. Nitti, Polytechnic University of Bari (2025)

clc
clear 
close all


%% input parameters
Lx=pi;                                          % domain length [m]
Ly=pi;                                          % domain length [m]
nx=41;                                          % number of nodes in x dir.
ny=41;                                          % number of nodes in y dir.
ua=@(x,y) cos(2.*x).*cos(2.*y);                 % analytical solution
rhs=@(x,y) -8.*cos(2.*x).*cos(2.*y);            % right-hand-side function

dx = Lx/(nx-1); 
dy = Ly/(ny-1);
xp = dx/2:dx:Lx-dx/2;
yp = dy/2:dy:Ly-dy/2;
[Xp,Yp] = meshgrid(xp,yp); 

f = rhs(Xp,Yp);

Mw = poisson2D_DCT_pre(Lx,Ly,nx,ny);

u = poisson2D_DCT(Mw,f);


%% plot result
figure(1)
[Xp,Yp]=meshgrid(xp,yp);
ua2=ua(Xp,Yp);
surf(Xp,Yp,ua2)
shading interp
hold on
surf(Xp,Yp,u,'FaceColor','none')
xlabel('x');
ylabel('y');
zlabel('\phi');
legend('analytical','numerical')


%% compute error
err=(ua2-u)./(ua2+1.0e-6);
disp(['relative error rms: ',num2str(rms(err,'all'))])



%% run convergence analysis
nn = [21,41,81];

Lx=pi;                                          % domain length [m]
Ly=Lx;                                          % domain length [m]

for ii=1:length(nn)

    nx = nn(ii);
    ny = nx;

    dx = Lx/(nx-1); 
    dy = dx;
    xp = dx/2:dx:Lx-dx/2;
    yp = xp;
    [Xp,Yp]=meshgrid(xp,yp);
    
    f = rhs(Xp,Yp);
    
    Mw = poisson2D_DCT_pre(Lx,Ly,nx,ny);
    
    u = poisson2D_DCT(Mw,f);
    
    % compute error
    ua2=ua(Xp,Yp);
    err = (ua2-u)./(ua2+1.0e-6);
    erri(ii) = rms(err,'all');

end

figure(2)
loglog(nn,erri,'-o')
hold on
o2 = 4.*nn.^-2;
loglog(nn,o2,'--')
xlabel('n');
ylabel('\epsilon');
legend('numerical','2nd order')
