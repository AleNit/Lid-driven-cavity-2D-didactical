
% Simulate the 2d Lid-Driven-Cavity test by solving incompressible 
% Navier-Stokes equations by means a fractional step method. Spatial derivatives 
% are discretized by 2nd-order accurate centered differences on a staggered
% grid. Time integration is carried out by a 2nd order Adams-Bashforth method.
% The steady solution at Re=400 is varified against the reference data from
% Ghia et al. (1982).
% A. Nitti, Polytechnic University of Bari (2025)

clc
clear 
close all
kf=1;


%% input parameters
Re = 400;                       % Reynolds number
T = 30;                         % simulation duration (dim.less time units)
Lx = 1;                         % domain length, x dir. 
Ly = 1;                         % domain length, y dir.
Ulid = 1;                       % lid velocity (dim.less)
nx = 101;                        % number of nodes in the x dir. 
ny = 101;                        % number of nodes in the y dir. 
CFLi = 0.5;                      % CFL control value
np=100;                          % step interval for visualization


%% pre-processing operations
dx = Lx/(nx-1);
dy = Ly/(ny-1);
dt = CFLi*min(dx,dy)/Ulid;      % choose time step size from CFL condition and lid velocity
nt = ceil(T/dt);                % number of time steps

% Coordinate of nodes and cell centers
xco = 0:dx:Lx;
xce = ( xco(2:end) + xco(1:end-1) ).*0.5;
yco = 0:dy:Ly;
yce = ( yco(2:end) + yco(1:end-1) ).*0.5;

% graphic variables
[Xce,Yce] = meshgrid(xce, yce);
mx=1:ceil(nx/15):nx;
my=1:ceil(ny/15):ny;

% initialization
time=0;
u = zeros(nx,ny+1);           % x velocity component
v = zeros(nx+1,ny);           % y velocity component
p = zeros(nx-1,ny-1);         % pressure 
Hxm = zeros(nx-2,ny-1);
Hym = zeros(nx-1,ny-2);

% compute wavenumber matrix
Mw = poisson2D_DCT_pre(Lx,Ly,nx,ny);


%% start time loop
tbet=1.0;
tgam=0.0;
talp=tbet+tgam;

for n = 1:nt

    % set Dirichelet boundary conditions on borders
    u(:,1) = -u(:,2);                       % south
    v(:,1) = 0;
    u(:,end) = 2*Ulid-u(:,end-1);           % north
    v(:,end) = 0;
    u(1,:) = 0;                             % west
    v(1,:) = -v(2,:);
    u(end,:) = 0;                           % east
    v(end,:) = -v(end-1,:);

    % viscous terms
    Vu = ( u(1:end-2,2:end-1) - 2*u(2:end-1,2:end-1) + u(3:end,2:end-1) )./dx^2 + ...  % nx-1 * ny
         ( u(2:end-1,1:end-2) - 2*u(2:end-1,2:end-1) + u(2:end-1,3:end) )./dy^2; 
    Vv = ( v(1:end-2,2:end-1) - 2*v(2:end-1,2:end-1) + v(3:end,2:end-1) )./dx^2 + ...  % nx * ny-1
         ( v(2:end-1,1:end-2) - 2*v(2:end-1,2:end-1) + v(2:end-1,3:end) )./dy^2; 

    % interpolate velocity at cell center/node
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))./2;
    uco = (u(:,1:end-1)+u(:,2:end))./2;
    vco = (v(1:end-1,:)+v(2:end,:))./2;
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))./2;
        
    % Convective terms
    uvco = uco.*vco;
    Cu = ( uce(2:end,:).^2 - uce(1:end-1,:).^2 )./dx;
    Cu = Cu + ( uvco(2:end-1,2:end) - uvco(2:end-1,1:end-1) )./dy;    
    Cv = ( vce(:,2:end).^2 - vce(:,1:end-1).^2 )./dy;
    Cv = Cv + ( uvco(2:end,2:end-1) - uvco(1:end-1,2:end-1) )./dx;

    % store convective and viscous terms
    Hx = -Cu + Vu./Re;
    Hy = -Cv + Vv./Re;

    % pressure gradient
    Gx = ( p(2:end,:)-p(1:end-1,:) )./dx;
    Gy = ( p(:,2:end)-p(:,1:end-1) )./dy;

    % advance momentum equation (provisional velocity field)
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + dt.*( -talp.*Gx + tbet.*Hx + tgam.*Hxm );
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) + dt.*( -talp.*Gy + tbet.*Hy + tgam.*Hym );

    % velocity divergence 
    div = ( u(2:end,2:end-1) - u(1:end-1,2:end-1) )./dx + ...
          ( v(2:end-1,2:end) - v(2:end-1,1:end-1) )./dy;
    
    % Solve for p (using cosine transform, faster)
    div = div./(dt*talp);
    phi = poisson2D_DCT(Mw,div);
    
    % compute phi gradient
    Gx = ( phi(2:end,:)-phi(1:end-1,:) )./dx;
    Gy = ( phi(:,2:end)-phi(:,1:end-1) )./dy;

    % update velocity field
    u(2:end-1,2:end-1) = u(2:end-1,2:end-1) - talp*dt.*Gx;
    v(2:end-1,2:end-1) = v(2:end-1,2:end-1) - talp*dt.*Gy;

    % update pressure field
    p = p + phi;
    
    % perform some checks
    div = ( u(2:end,2:end-1) - u(1:end-1,2:end-1) )./dx + ...
          ( v(2:end-1,2:end) - v(2:end-1,1:end-1) )./dy;
    dmax = max(max(div));
    maxudx = max(max(abs(u)))/dx;
    maxvdy = max(max(abs(v)))/dy;
    CFL = max(maxudx,maxvdy)*dt;

    % update variables
    time = time + dt;
    Hxm  = Hx;
    Hym  = Hy;    
    tbet = 1.5;
    tgam = -0.5;
    talp = tbet + tgam;
        
    % write header
    disp([  'step ',num2str(n),9, ...
            't = ',num2str(time),9, ...            
            'd_max = ',num2str(dmax),9, ...
            'CFL = ',num2str(CFL),9, ...
            ])
    

    % plot velocity magnitude and direction
    if ( mod(n,np)==0 )

        uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
        vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
        umod = sqrt( uce.^2 + vce.^2 );
        surf(Xce,Yce,zeros(size(Xce)),umod')
        shading interp
        hold on
        colormap(turbo)
        colorbar
        quiver(Xce(mx,my),Yce(mx,my),uce(mx,my)',vce(mx,my)','k')
        axis equal
        axis([0,Lx,0,Ly])
        xlabel('$x/L$','interpreter','latex');    
        ylabel('$y/L$','interpreter','latex');                    
        tit=strcat('$|\mathbf{u}|, \, t \, U/L=',num2str(time,'%2.2f'),'$');
        title(tit,'interpreter','latex')
        view(2)
        drawnow
        hold off

    end

end



%% compare result with reference solution
fileID = fopen('./Ghia_YU.txt','r');
uref = fscanf(fileID,'%f',[2 Inf])';
fclose(fileID);

fileID = fopen('./Ghia_XV.txt','r');
vref = fscanf(fileID,'%f',[2 Inf])';
fclose(fileID);

up = u(ceil((ny-1)/2),1:end-1);
vp = v(1:end-1,ceil((nx-1)/2));

figure(2)
plot(yco,up,'-k')
hold on
plot(uref(:,1),uref(:,2),'sr')
xlabel('$y/L$','Interpreter','latex','FontSize',14)
ylabel('$u_x/U$','Interpreter','latex','FontSize',14)
legend('present work','Ghia et al. (1982) - $Re=400$','Interpreter','latex','FontSize',13)

figure(3)
plot(xco,vp,'-k')
hold on
plot(vref(:,1),vref(:,2),'sr')
xlabel('$x/L$','Interpreter','latex','FontSize',14)
ylabel('$u_y/U$','Interpreter','latex','FontSize',14)
legend('present work','Ghia et al. (1982) - $Re=400$','Interpreter','latex','FontSize',13)
