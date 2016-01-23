close all;
clear all;
format short g;
clc;
%-----------------------------------------------------------------------
Re = 1e2;     % Reynolds number
dt = 1e-2;    % time step
tf = 5.0;    % final time
xStart = 0.0;
xEnd = 1.0;
yStart = 0.0;
yEnd = 1.0;
nx = 400;      % number of x-gridpoints
ny = 400;      % number of y-gridpoints
nsteps = 10;  % number of steps with graphic output
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(xStart,xEnd,nx+1); hx = (xEnd - xStart)/nx;
y = linspace(yStart,yEnd,ny+1); hy = (yEnd - yStart)/ny;
[Y,X] = meshgrid(y,x);
%-----------------------------------------------------------------------
% initial conditions
U = zeros(nx-1,ny) + eps; V = zeros(nx,ny-1) + eps; Pressure = zeros(nx,ny) + eps; P = Pressure + eps;
% boundary conditions
uN = x*0+0;    vN = avg(x)*0;
uS = x*0+0;      vS = avg(x)*0;
uW = avg(y)*0+1; vW = y*0+0;
uE = avg(y)*0+1; vE = y*0+0;

%-----------------------------------------------------------------------
% Free-slip boundary condition on top and bottom walls
uN(2:end-1) = U(:,end);
uS(2:end-1) = U(:,1);

% Zero gradient for velocity on east wall
uE = (mean(uW) ./ mean(U(end,:))) .* U(end,:);
vE(2:end-1) = (mean(vW) ./ mean(V(end,:))) .* V(end,:);
%-----------------------------------------------------------------------
Ubc = dt/Re*...
    (...
    [2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)'] / hx^2 + ...
    [uW;zeros(nx-3,ny);uE] / hy^2 ...
    );

Vbc = dt/Re*...
    (...
    [vS' zeros(nx,ny-3) vN'] / hx^2 + ...
    [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)] / hy^2 ...
    );

fprintf('initialization (governing equations solution) \n')
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));
Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp); Rp = chol(Lp(perp,perp)); Rpt = Rp';

Lu = speye((nx-1)*ny)+dt/Re*(kron(speye(ny),K1(nx-1,hx,2)) + ...
           kron(K1(ny,hy,3),speye(nx-1)));
peru = symamd(Lu); Ru = chol(Lu(peru,peru)); Rut = Ru';

Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3)) + ...
           kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv); Rv = chol(Lv(perv,perv)); Rvt = Rv';

Lq = kron(speye(ny-1),K1(nx-1,hx,2))+kron(K1(ny-1,hy,2),speye(nx-1));
perq = symamd(Lq); Rq = chol(Lq(perq,perq)); Rqt = Rq';

fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')

Uhist = sparse((nx-1)*ny,nt); % Time history of U
Vhist = sparse(nx*(ny-1),nt); % Time history of V

% Immersed boundary
% N = 100; % Number of Lagrange points
% xImm = 0.5 + 0.12 * cos(linspace(0,2*pi,N));
% yImm = 0.5 + 0.12 * sin(linspace(0,2*pi,N));
% dxImm_db = cos(linspace(0,2*pi,N));
% dyImm_db = sin(linspace(0,2*pi,N));
pointCloud = dlmread('point_cloud.txt');
pointCloudSensitivity = dlmread('point_cloud_sensitivity.txt');
xImm = pointCloud(:,1);
yImm = pointCloud(:,2);
dxImm_db = pointCloudSensitivity(:,1);
dyImm_db = pointCloudSensitivity(:,2);

[deltaMatU,deltaMatUSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'U');
[deltaMatV,deltaMatVSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'V');
[deltaMatP,deltaMatPSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'P');

deltaMatU = sparse(deltaMatU);
deltaMatUSensitivity = sparse(deltaMatUSensitivity);
deltaMatV = sparse(deltaMatV);
deltaMatVSensitivity = sparse(deltaMatVSensitivity);
deltaMatP = sparse(deltaMatP);
deltaMatPSensitivity = sparse(deltaMatPSensitivity);

alpha = -1e4;
beta = -5e1;
% break
T = linspace(0,tf,nt);

fxHist = zeros(size(xImm,1),length(T));
fyHist = zeros(size(yImm,1),length(T));
fpxHist = zeros(size(xImm,1),length(T));
fpyHist = zeros(size(yImm,1),length(T));

timeIntegralx = 0;
timeIntegraly = 0;
for k = 1:nt
   % treat nonlinear terms
%    gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   gamma = 0.0;
   
   % Free-slip boundary condition on top and bottom walls
   uN(2:end-1) = U(:,end);
   uS(2:end-1) = U(:,1);

   % Zero gradient for velocity on east wall
   uE = (mean(uW) ./ mean(U(end,:))) .* U(end,:);
   vE(2:end-1) = (mean(vW) ./ mean(V(end,:))) .* V(end,:);

   Ubc = dt/Re*...
    (...
    [2*uS(2:end-1)' zeros(nx-1,ny-2) 2*uN(2:end-1)'] / hx^2 + ...
    [uW;zeros(nx-3,ny);uE] / hy^2 ...
    );

   Vbc = dt/Re*...
    (...
    [vS' zeros(nx,ny-3) vN'] / hx^2 + ...
    [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)] / hy^2 ...
    );

   Ue = [uW;U;uE]; Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
   Ve = [vS' V vN']; Ve = [2*vW-Ve(1,:);Ve;2*vE-Ve(end,:)];
   
   % Calculate convective terms based on upwind approximation
   % upwindMatrix = upwindFiniteDifference(upwindMatrix,flowDirection,originalDirection,neededDirection,gamma)
   U_xMOMdx = upwindFiniteDifference(Ue,Ue,'x','x',gamma);
   U_xMOMdy = upwindFiniteDifference(Ue,Ve,'x','y',gamma);
   V_xMOMdy = upwindFiniteDifference(Ve,Ue,'y','x',gamma);
   
   U_yMOMdx = upwindFiniteDifference(Ue,Ve,'x','y',gamma);
   V_yMOMdx = upwindFiniteDifference(Ve,Ue,'y','x',gamma);
   V_yMOMdy = upwindFiniteDifference(Ve,Ve,'y','y',gamma);
   
   dUUdX = diff(U_xMOMdx(:,2:end-1) .* U_xMOMdx(:,2:end-1)) / hx;
   dUVdY = diff([U_xMOMdy(2:end-1,:) .* V_xMOMdy(2:end-1,:)]')' / hy;
   dVUdX = diff(U_yMOMdx(:,2:end-1) .* V_yMOMdx(:,2:end-1)) / hx;
   dVVdY = diff([V_yMOMdy(2:end-1,:) .* V_yMOMdy(2:end-1,:)]')' / hy;
   % Calculate force terms
   if k > 1
%        fx = alpha * deltaMatU' * trapz(T(1:k),(deltaMatU * Uhist(:,1:k) - 0)')' + ...
%            beta * deltaMatU' * (deltaMatU * (reshape(U,[],1) - 0));
%        fy = alpha * deltaMatV' * trapz(T(1:k),(deltaMatV * Vhist(:,1:k) - 0)')' + ...
%            beta * deltaMatV' * (deltaMatV * (reshape(V,[],1) - 0));
       timeIntegralx = timeIntegralx + alpha * transpose(deltaMatU) * deltaMatU * (reshape(U,[],1) + Uhist(:,k-1)) * dt / 2;
       fx = timeIntegralx + ...
           beta * transpose(deltaMatU) * (deltaMatU * (reshape(U,[],1) - 0));
       timeIntegraly = timeIntegraly + alpha * transpose(deltaMatV) * deltaMatV * (reshape(V,[],1) + Vhist(:,k-1)) * dt / 2;
       fy = timeIntegraly + ...
           beta * transpose(deltaMatV) * (deltaMatV * (reshape(V,[],1) - 0));
       fxHist(:,k) = deltaMatU * fx;
       fyHist(:,k) = deltaMatV * fy;
   else
       fx = 0.0;
       fy = 0.0;
   end
   
%    rhs = reshape(U - dt * (dUUdX + dUVdY) + Ubc - diff(P) / hx,[],1) + dt * fx;
   rhs = reshape(U - dt * (dUUdX + dUVdY) + Ubc,[],1) + dt * fx;
   u(peru) = Ru\(Rut\rhs(peru));
   U = reshape(u,nx-1,ny);

%    rhs = reshape(V - dt * (dVUdX + dVVdY) + Vbc  - diff(P')' / hy,[],1) + dt * fy;
   rhs = reshape(V - dt * (dVUdX + dVVdY) + Vbc,[],1) + dt * fy;
   v(perv) = Rv\(Rvt\rhs(perv));
   V = reshape(v,nx,ny-1);
   
   % pressure correction
   rhs = reshape(diff([uW;U;uE])/hx+diff([vS' V vN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   P = reshape(p,nx,ny);
   U = U-diff(P)/hx;
   V = V-diff(P')'/hy;
   Pressure = Pressure + P;
   P = Pressure;
   % Save time history of U and V
   Uhist(:,k) = reshape(U,[],1);
   Vhist(:,k) = reshape(V,[],1);
   
   % visualization
   if floor(25*k/nt)>floor(25*(k-1)/nt)
       fprintf('.')
   end
   if max(max(isnan(U))) || max(max(isnan(V)))
       fprintf('\n')
       disp('NaN detected!');
       fprintf('\n')
       disp(['time: ' num2str(k*dt)]);
       break
   end
end
fprintf('\n')
%=======================================================================
xu = linspace(xStart,xEnd,size(U,1));
yu = linspace(yStart,yEnd,size(U,2));
[Yu,Xu] = meshgrid(yu,xu);
figure,
contourf(Xu,Yu,U,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal

xp = linspace(xStart,xEnd,size(P,1));
yp = linspace(yStart,yEnd,size(P,2));
[Yp,Xp] = meshgrid(yp,xp);

figure,
contourf(Xp,Yp,P,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal

uBoundary = deltaMatU * reshape(U,[],1);
vBoundary = deltaMatV * reshape(V,[],1);
pBoundary = deltaMatP * reshape(P,[],1);

figure,plot(linspace(0,2*pi,length(xImm)),uBoundary);
figure,plot(linspace(0,2*pi,length(xImm)),vBoundary);
% Uavg = avg(U')';
% Vavg = avg(V);
% figure,
% contourf(sqrt(Uavg.^2 + Vavg.^2))
% break
% ----------------------------------------------------------------------- %
% %{
% Sensitivity Analysis
% initial conditions
Up = zeros(nx-1,ny) + eps; Vp = zeros(nx,ny-1) + eps; Pressurep = zeros(nx,ny) + eps; Pp = Pressurep;
% boundary conditions
UpN = x*0+0;    VpN = avg(x)*0;
UpS = x*0+0;      VpS = avg(x)*0;
UpW = avg(y)*0+0; VpW = y*0+0;
UpE = avg(y)*0+0; VpE = y*0+0;
%-----------------------------------------------------------------------
% Free-slip boundary condition on top and bottom walls
UpN(2:end-1) = Up(:,end);
UpS(2:end-1) = Up(:,1);

% Zero gradient for velocity on east wall
UpE = (mean(UpW) ./ mean(Up(end,:))) .* Up(end,:);
VpE(2:end-1) = (mean(VpW) ./ mean(Vp(end,:))) .* Vp(end,:);
%-----------------------------------------------------------------------
Upbc = dt/Re*...
    (...
    [2*UpS(2:end-1)' zeros(nx-1,ny-2) 2*UpN(2:end-1)'] / hx^2 + ...
    [UpW;zeros(nx-3,ny);UpE] / hy^2 ...
    );

Vpbc = dt/Re*...
    (...
    [VpS' zeros(nx,ny-3) VpN'] / hx^2 + ...
    [2*VpW(2:end-1);zeros(nx-2,ny-1);2*VpE(2:end-1)] / hy^2 ...
    );

fprintf('initialization (sensitivity analysis) \n')

% Convective operator for U and V
% upwindMatrix = upwindFiniteDifference(upwindMatrix,flowDirection,originalDirection,neededDirection,gamma)
Upe = [UpW;Up;UpE]; Upe = [2*UpS'-Upe(:,1) Upe 2*UpN'-Upe(:,end)];
Vpe = [VpS' Vp VpN']; Vpe = [2*VpW-Vpe(1,:);Vpe;2*VpE-Vpe(end,:)];

U_xMOMdx = upwindFiniteDifference(Ue,Ue,'x','x',gamma);
U_xMOMdy = upwindFiniteDifference(Ue,Ve,'x','y',gamma);
V_xMOMdy = upwindFiniteDifference(Ve,Ue,'y','x',gamma);

U_yMOMdx = upwindFiniteDifference(Ue,Ve,'x','y',gamma);
V_yMOMdx = upwindFiniteDifference(Ve,Ue,'y','x',gamma);
V_yMOMdy = upwindFiniteDifference(Ve,Ve,'y','y',gamma);

% dUUdX = diff(U_xMOMdx(:,2:end-1) .* U_xMOMdx(:,2:end-1)) / hx;
% dUVdY = diff([U_xMOMdy(2:end-1,:) .* V_xMOMdy(2:end-1,:)]')' / hy;
% dVUdX = diff(U_yMOMdx(:,2:end-1) .* V_yMOMdx(:,2:end-1)) / hx;
% dVVdY = diff([V_yMOMdy(2:end-1,:) .* V_yMOMdy(2:end-1,:)]')' / hy;

fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')

Uphist = sparse((nx-1)*ny,nt); % Time history of dU/dB
Vphist = sparse(nx*(ny-1),nt); % Time history of dV/dB

timeIntegralx = 0;
timeIntegraly = 0;
for k = 1:nt
%    k
   % treat nonlinear terms
%    gamma = min(1.2*dt*max(max(max(abs(Up)))/hx,max(max(abs(Vp)))/hy),1);
   gamma = 0.0;
%    tic

       % ADDED THIS ----------------------------------------------------------
   % Free-slip boundary condition on top and bottom walls
    UpN(2:end-1) = Up(:,end);
    UpS(2:end-1) = Up(:,1);

    % Zero gradient for velocity on east wall
%     UpE = (mean(UpW) ./ mean(Up(end,:))) .* Up(end,:);
%     VpE(2:end-1) = (mean(VpW) ./ mean(Vp(end,:))) .* Vp(end,:);
    UpE = 1 .* Up(end,:);
    VpE(2:end-1) = 1 .* Vp(end,:);
    % ----------------------------------------------------------------- %
   Upe = [UpW;Up;UpE]; Upe = [2*UpS'-Upe(:,1) Upe 2*UpN'-Upe(:,end)];
   Vpe = [VpS' Vp VpN']; Vpe = [2*VpW-Vpe(1,:);Vpe;2*VpE-Vpe(end,:)];
   
   Up_xMOMdx = upwindFiniteDifference(Upe,Upe,'x','x',gamma);
   Up_xMOMdy = upwindFiniteDifference(Upe,Vpe,'x','y',gamma);
   Vp_xMOMdy = upwindFiniteDifference(Vpe,Upe,'y','x',gamma);

   Up_yMOMdx = upwindFiniteDifference(Upe,Vpe,'x','y',gamma);
   Vp_yMOMdx = upwindFiniteDifference(Vpe,Upe,'y','x',gamma);
   Vp_yMOMdy = upwindFiniteDifference(Vpe,Vpe,'y','y',gamma);

   dUUpdX = diff(U_xMOMdx(:,2:end-1) .* Up_xMOMdx(:,2:end-1)) / hx;
   dUVpdY = diff([U_xMOMdy(2:end-1,:) .* Vp_xMOMdy(2:end-1,:)]')' / hy;
   dUpVdY = diff([Up_xMOMdy(2:end-1,:) .* V_xMOMdy(2:end-1,:)]')' / hy;
   
   dVUpdX = diff(Up_yMOMdx(:,2:end-1) .* V_yMOMdx(:,2:end-1)) / hx;
   dVpUdX = diff(U_yMOMdx(:,2:end-1) .* Vp_yMOMdx(:,2:end-1)) / hx;
   dVVpdY = diff([V_yMOMdy(2:end-1,:) .* Vp_yMOMdy(2:end-1,:)]')' / hy;
%    toc
   % Calculate force terms
%    tic
   if k > 1
%    fpx = alpha * trapz(T(1:k),(deltaMatUSensitivity' * deltaMatU * Uhist(:,1:k) + ...
%                           deltaMatU' * deltaMatUSensitivity * Uhist(:,1:k) + ...
%                           deltaMatU' * deltaMatU * Uphist(:,1:k))')' + ...
%          beta * (deltaMatUSensitivity' * deltaMatU * reshape(U,[],1) + ...
%                  deltaMatU' * deltaMatUSensitivity * reshape(U,[],1) + ...
%                  deltaMatU' * deltaMatU * reshape(Up,[],1));
%              
%    fpy = alpha * trapz(T(1:k),(deltaMatVSensitivity' * deltaMatV * Vhist(:,1:k) + ...
%                           deltaMatV' * deltaMatVSensitivity * Vhist(:,1:k) + ...
%                           deltaMatV' * deltaMatV * Vphist(:,1:k))')' + ...
%          beta * (deltaMatVSensitivity' * deltaMatV * reshape(V,[],1) + ...
%                  deltaMatV' * deltaMatVSensitivity * reshape(V,[],1) + ...
%                  deltaMatV' * deltaMatV * reshape(Vp,[],1));
       timeIntegralx = timeIntegralx + ...
           alpha * transpose(deltaMatUSensitivity) * deltaMatU * (reshape(U,[],1) + Uhist(:,k-1)) * dt + ...
           alpha * transpose(deltaMatU) * deltaMatUSensitivity * (reshape(U,[],1) + Uhist(:,k-1)) * dt + ...
           alpha * transpose(deltaMatU) * deltaMatU * (reshape(Up,[],1) + Uphist) * dt;
       
       fpx = timeIntegralx + ...
           beta * transpose(deltaMatUSensitivity) * (deltaMatU * (reshape(U,[],1) - 0)) + ...
           beta * transpose(deltaMatU) * (deltaMatUSensitivity * (reshape(U,[],1) - 0)) + ...
           beta * transpose(deltaMatU) * (deltaMatU * (reshape(Up,[],1) - 0));
       
       timeIntegraly = timeIntegraly + ...
           alpha * transpose(deltaMatVSensitivity) * deltaMatV * (reshape(V,[],1) + Vhist(:,k-1)) * dt + ...
           alpha * transpose(deltaMatV) * deltaMatVSensitivity * (reshape(V,[],1) + Vhist(:,k-1)) * dt + ...
           alpha * transpose(deltaMatV) * deltaMatV * (reshape(Vp,[],1) + Vphist) * dt;
       
       fpy = timeIntegraly + ...
           beta * transpose(deltaMatVSensitivity) * (deltaMatV * (reshape(V,[],1) - 0)) + ...
           beta * transpose(deltaMatV) * (deltaMatVSensitivity * (reshape(V,[],1) - 0)) + ...
           beta * transpose(deltaMatV) * (deltaMatV * (reshape(Vp,[],1) - 0));
       
       fpxHist(:,k) = deltaMatU * fpx;
       fpyHist(:,k) = deltaMatV * fpy;
   else
       fpx = 0;
       fpy = 0;
   end
%    toc
%    tic
%    rhs = reshape(Up - dt * (2 * dUUpdX + dUVpdY + dUpVdY) + Upbc - diff(Pp) / hx,[],1) + dt * fpx;
   rhs = reshape(Up - dt * (2 * dUUpdX + dUVpdY + dUpVdY) + Upbc,[],1) + dt * fpx;
   u(peru) = Ru\(Rut\rhs(peru));
   Up = reshape(u,nx-1,ny);
   
%    rhs = reshape(Vp - dt * (dVUpdX + dVpUdX + 2 * dVVpdY) + Vpbc - diff(Pp')' / hy,[],1) + dt * fpy;
   rhs = reshape(Vp - dt * (dVUpdX + dVpUdX + 2 * dVVpdY) + Vpbc,[],1) + dt * fpy;
   v(perv) = Rv\(Rvt\rhs(perv));
   Vp = reshape(v,nx,ny-1);
%    toc
   % pressure correction
%    tic
   rhs = reshape(diff([UpW;Up;UpE])/hx+diff([VpS' Vp VpN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   Pp = reshape(p,nx,ny);
   Up = Up - diff(Pp)/hx;
   Vp = Vp - diff(Pp')'/hy;
   Pressurep = Pressurep + Pp;
   Pp = Pressurep;
%    toc
   % Save time history of U and V
%    tic
   Uphist = reshape(Up,[],1);
   Vphist = reshape(Vp,[],1);
%    Uphist(:,k) = reshape(Up,[],1);
%    Vphist(:,k) = reshape(Vp,[],1);
%    toc
   % visualization
   if floor(25*k/nt)>floor(25*(k-1)/nt)
       fprintf('.')
%        try
%            clf(fig)
%        catch
%        end
%        figure(4),
%        contourf(Xu,Yu,Up,50,'linestyle','none')
%        hold all
%        plot(xImm,yImm,'w.','linewidth',2)
%        set(gca,'XTickLabel','','YTickLabel','')
%        axis equal
%        fig = gcf;
   end
   if max(max(isnan(Up))) || max(max(isnan(Vp))) || max(max(isnan(Pp)))
       fprintf('\n')
       disp('NaN detected!');
       fprintf('\n')
       disp(['time: ' num2str(k*dt)]);
       break
   end
end
fprintf('\n')

xu = linspace(xStart,xEnd,size(U,1));
yu = linspace(yStart,yEnd,size(U,2));
[Yu,Xu] = meshgrid(yu,xu);
figure,
contourf(Xu,Yu,Up,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal

xv = linspace(xStart,xEnd,size(V,1));
yv = linspace(yStart,yEnd,size(V,2));
[Yv,Xv] = meshgrid(yv,xv);
figure,
contourf(Xv,Yv,Vp,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal

xp = linspace(xStart,xEnd,size(P,1));
yp = linspace(yStart,yEnd,size(P,2));
[Yp,Xp] = meshgrid(yp,xp);
figure,
contourf(Xp,Yp,Pp,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal

upBoundary = deltaMatU * reshape(Up,[],1);
vpBoundary = deltaMatV * reshape(Vp,[],1);
ppBoundary = deltaMatP * reshape(Pp,[],1);

figure,
subplot(2,1,1)
plot(xImm,pBoundary)
subplot(2,1,2)
plot(xImm,ppBoundary)
title('pressure')
trapz(xImm,pBoundary)

dlmwrite('mesh_convergence/ppBoundary.txt', ppBoundary)
dlmwrite('mesh_convergence/Up.txt', Up)
dlmwrite('mesh_convergence/Vp.txt', Vp)
dlmwrite('mesh_convergence/Pp.txt', Pp)

% Save data
% save('current_simulation_results/U_0','U')
% save('current_simulation_results/V_0','V')
% save('current_simulation_results/P_0','P')
% save('current_simulation_results/uBoundary_0','uBoundary')
% save('current_simulation_results/vBoundary_0','vBoundary')
% save('current_simulation_results/pBoundary_0','pBoundary')
% save('current_simulation_results/fyHist_0','fyHist')
% save('current_simulation_results/fxHist_0','fxHist')
% save('current_simulation_results/deltaMatU','deltaMatU')
% save('current_simulation_results/deltaMatV','deltaMatV')
% save('current_simulation_results/deltaMatP','deltaMatP')
% 
% save('current_simulation_results/Up_CSA001_','Up')
% save('current_simulation_results/Vp_CSA001_','Vp')
% save('current_simulation_results/Pp_CSA001_','Pp')
% save('current_simulation_results/ppBoundary_CSA001_','ppBoundary')
% save('current_simulation_results/upBoundary_CSA001_','upBoundary')
% save('current_simulation_results/vpBoundary_CSA001_','vpBoundary')
% save('current_simulation_results/fpyHist_CSA001_','fpyHist')
% save('current_simulation_results/fpxHist_CSA001_','fpxHist')
% figure,plot(linspace(0,2*pi,length(xImm)),upBoundary);
% figure,plot(linspace(0,2*pi,length(xImm)),vpBoundary);
%}
