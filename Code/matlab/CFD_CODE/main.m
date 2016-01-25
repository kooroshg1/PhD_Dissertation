close all;
clear all;
format short g;
clc;
%-----------------------------------------------------------------------
fontsize = 40;
linewidth = 5;
markersize = 20;
%-----------------------------------------------------------------------
Re = 1e2;     % Reynolds number
dt = 1e-2;    % time step
tf = 5;    % final time
xStart = 0.0;
xEnd = 4.0;
yStart = -1.5;
yEnd = 1.5;
nx = 400;      % number of x-gridpoints
ny = 300;      % number of y-gridpoints
nsteps = 10;  % number of steps with graphic output
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;
x = linspace(xStart,xEnd,nx+1); hx = (xEnd - xStart)/nx;
y = linspace(yStart,yEnd,ny+1); hy = (yEnd - yStart)/ny;
[Y,X] = meshgrid(y,x);
%-----------------------------------------------------------------------
% initial conditions
U = zeros(nx-1,ny) + eps; V = zeros(nx,ny-1) + eps; Pressure = ones(nx,ny); P = Pressure;
% boundary conditions
uN = x*0+0;      vN = avg(x)*0;
uS = x*0+0;      vS = avg(x)*0;
uW = avg(y)*0+1; vW = y*0+0;
uE = avg(y)*0+0; vE = y*0+0;

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

fprintf('initialization \n')
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

% Uhist = zeros((nx-1)*ny,nt); % Time history of U
% Vhist = zeros(nx*(ny-1),nt); % Time history of V

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

dlmwrite('point_cloud_velocity.txt',zeros(size(pointCloud,1),size(pointCloud,2)));

% [deltaMatU,deltaMatUSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'U');
% [deltaMatV,deltaMatVSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'V');
% [deltaMatP,deltaMatPSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'P');
% 
% deltaMatU = sparse(deltaMatU);
% deltaMatUSensitivity = sparse(deltaMatUSensitivity);
% deltaMatV = sparse(deltaMatV);
% deltaMatVSensitivity = sparse(deltaMatVSensitivity);
% deltaMatP = sparse(deltaMatP);
% deltaMatPSensitivity = sparse(deltaMatPSensitivity);

maxSubIteration = 50;

alpha = -1e4;
beta = -5e1;
% break
% T = linspace(0,tf,nt);
% fx = 0;
% fy = 0;
% fxHist = zeros(size(xImm,1),length(T) * 10);
% fyHist = zeros(size(yImm,1),length(T) * 10);

timeIntegralx = 0;
timeIntegraly = 0;

pointCloudVelocity = dlmread('point_cloud_velocity.txt');

UBoundary = pointCloudVelocity(:,1); % Velocity of boundary
VBoundary = pointCloudVelocity(:,2); % Velocity of boundary

[nodeEffectiveArea,nodeUnitNormalVector] = calcNodeAreaOrientation(pointCloud);

centerLocation = sparse(nt,1);
figureName = 1.0;
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

pointCloudVelocity = dlmread('point_cloud_velocity.txt');


UBoundary = pointCloudVelocity(:,1); % Velocity of boundary
VBoundary = pointCloudVelocity(:,2); % Velocity of boundary

movePointcloud(deltaMatP * reshape(P,[],1),dt,nodeEffectiveArea,nodeUnitNormalVector)

fx = 0.0;
fy = 0.0;
timeIntegralx = 0;
timeIntegraly = 0;

fxHist = zeros(length(xImm), nt);
fyHist = zeros(length(yImm), nt);

% centerLocation(timeIndex) = mean(pointCloud(:,1));
xu = linspace(xStart,xEnd,size(U,1));
yu = linspace(yStart,yEnd,size(U,2));
[Yu,Xu] = meshgrid(yu,xu);
[Uin, Uon] = inpolygon(Xu, Yu, xImm, yImm); 

for k = 1:nt
%         continue
%     for subIteration = 1:maxSubIteration
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
            timeIntegralx = timeIntegralx + ...
               alpha * transpose(deltaMatU) * ( ...
                   deltaMatU * reshape(U,[],1)  - UBoundary + ...
                   deltaMatU * Ulast - UBoundary ...
                   ) * dt / 2;
           fx = timeIntegralx + ...
               beta * transpose(deltaMatU) * ( ...
                   deltaMatU * reshape(U,[],1) - UBoundary...
                   );

           timeIntegraly = timeIntegraly + ...
               alpha * transpose(deltaMatV) * ( ...
                   deltaMatV * reshape(V,[],1)  - VBoundary + ...
                   deltaMatV * Vlast - VBoundary ...
                   ) * dt / 2;
           fy = timeIntegraly + ...
               beta * transpose(deltaMatV) * ( ...
                   deltaMatV * reshape(V,[],1) - VBoundary...
                   );    
               
           fxHist(:,k) = deltaMatU * fx;
           fyHist(:,k) = deltaMatV * fy;
       else
           fx = 0.0;
           fy = 0.0;
       end

       
%        rhs = reshape(U - dt * (dUUdX + dUVdY) + Ubc - diff(P) / hx,[],1) + dt * fx;
       rhs = reshape(U - dt * (dUUdX + dUVdY) + Ubc,[],1) + dt * fx;
       u(peru) = Ru\(Rut\rhs(peru));
       U = reshape(u,nx-1,ny);

%        rhs = reshape(V - dt * (dVUdX + dVVdY) + Vbc  - diff(P')' / hy,[],1) + dt * fy;
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
%        Uhist(:,k) = reshape(U,[],1);
%        Vhist(:,k) = reshape(V,[],1);
       Ulast = reshape(U,[],1);
       Vlast = reshape(V,[],1);
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
%        if rem(k,50) == 0
%            Uplot = U;
%            Uplot(Uin) = nan;
%            figure(1),
%            contourf(Xu,Yu,Uplot,50,'linestyle','none')
%            colorbar;
%            hold all
%            plot(xImm,yImm,'w.','linewidth',2)
%            set(gca,'XTickLabel','','YTickLabel','')
%            axis equal
%            title(['Time:' num2str(k)])
%            print(['fig/' num2str(figureName)],'-dpng')
%            figureName = figureName + 1;
%            pause(0.01)
%            clf
% %     %        mean(pointCloud(:,1))
% %            dlmwrite(['surfacePressureHistory/' num2str(figureName)],deltaMatP * reshape(P,[],1));
%        end
end
fprintf('\n')
Ux = diff(U);
Vy = diff(V')';
diver = (Ux(:,2:end-1)/hx + Vy(2:end-1,:)/hy);
max(max(abs(diver)))
%=======================================================================
xv = linspace(xStart,xEnd,size(V,1));
yv = linspace(yStart,yEnd,size(V,2));
[Yv,Xv] = meshgrid(yv,xv);

figure,
contourf(Xv,Yv,V,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal
title('v-velocity')

xu = linspace(xStart,xEnd,size(U,1));
yu = linspace(yStart,yEnd,size(U,2));
[Yu,Xu] = meshgrid(yu,xu);

figure,
contourf(Xu,Yu,U,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal
title('u-velocity')

xp = linspace(xStart,xEnd,size(P,1));
yp = linspace(yStart,yEnd,size(P,2));
[Yp,Xp] = meshgrid(yp,xp);

figure,
contourf(Xp,Yp,Pressure,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal
title('pressure')

% Uavg = avg(U);
% Vavg = avg(V')';
% Umag = (Uavg(:,2:end-1).^2 + Vavg(2:end-1,:).^2);
% xavg = linspace(xStart+hx/2,xEnd-hx/2,size(Umag,1));
% yavg = linspace(yStart+hy/2,yEnd-hy/2,size(Umag,2));
% [Yavg,Xavg] = meshgrid(yavg,xavg);
% 
% figure,
% contourf(Xavg,Yavg,Umag,50,'linestyle','none')
% hold all
% plot(xImm,yImm,'w.','linewidth',2)
% set(gca,'XTickLabel','','YTickLabel','')
% axis equal
% title('velocity magnitude')

uBoundary = deltaMatU * reshape(U,[],1);
vBoundary = deltaMatV * reshape(V,[],1);
pBoundary = deltaMatP * reshape(P,[],1);
xpBoundary = deltaMatP * reshape(Xp,[],1);
ypBoundary = deltaMatP * reshape(Yp,[],1);

% figure,plot(linspace(0,360,length(xImm)),sgolayfilt(uBoundary,3,51),'k','linewidth',5);
% xlabel('Angular location on cylinder','fontsize',fontsize)
% ylabel('U velocity','fontsize',fontsize)
% set(gca,'fontsize',fontsize)
% 
% figure,plot(linspace(0,360,length(xImm)),sgolayfilt(vBoundary,3,51),'k','linewidth',5);
% xlabel('Angular location on cylinder','fontsize',fontsize)
% ylabel('V velocity','fontsize',fontsize)
% set(gca,'fontsize',fontsize)
% 
% figure,plot(linspace(0,360,length(xImm)),sgolayfilt(pBoundary,3,51),'k','linewidth',5);
% xlabel('Angular location on cylinder','fontsize',fontsize)
% ylabel('Pressure','fontsize',fontsize)
% set(gca,'fontsize',fontsize)

filterDegree = 17;
filterSize = 41;
figure,scatter3(xImm,yImm,uBoundary),title('u')
% figure,plot(xImm,pBoundary),title('pressure')
figure,
% plot(xImm,sgolayfilt(pBoundary,filterDegree,filterSize),...
%      xImm,pBoundary)
plot(xImm,sgolayfilt(pBoundary,filterDegree,filterSize))
title('pressure')
trapz(xImm,pBoundary)

figure,
plot(xImm,uBoundary,'k','linewidth',linewidth)
xlabel('X','fontsize',fontsize)
ylabel('U-velocity','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
plot(xImm,vBoundary,'k','linewidth',linewidth)
xlabel('X','fontsize',fontsize)
ylabel('V-velocity','fontsize',fontsize)
set(gca,'fontsize',fontsize)


% Save data
% save('current_simulation_results/U_0000001','U')
% save('current_simulation_results/V_0000001','V')
% save('current_simulation_results/P_0000001','P')
% save('current_simulation_results/uBoundary_0000001','uBoundary')
% save('current_simulation_results/vBoundary_0000001','vBoundary')
% save('current_simulation_results/pBoundary_0000001','pBoundary')
% save('current_simulation_results/fyHist_0000001','fyHist')
% save('current_simulation_results/fxHist_0000001','fxHist')
% save P
% Uavg = avg(U')';
% Vavg = avg(V);
% figure,
% contourf(sqrt(Uavg.^2 + Vavg.^2))

dlmwrite('cylinder/xImm', xImm)
dlmwrite('cylinder/yImm', yImm)
dlmwrite('cylinder/Xu.txt', Xu)
dlmwrite('cylinder/Yu.txt', Yu)
dlmwrite('cylinder/Xv.txt', Xv)
dlmwrite('cylinder/Yv.txt', Yv)
dlmwrite('cylinder/Xp.txt', Xp)
dlmwrite('cylinder/Yp.txt', Yp)
dlmwrite('cylinder/U.txt', U)
dlmwrite('cylinder/V.txt', V)
dlmwrite('cylinder/P.txt', P)
dlmwrite('cylinder/pBoundary.txt', pBoundary)
dlmwrite('cylinder/uBoundary.txt', uBoundary)
dlmwrite('cylinder/vBoundary.txt', vBoundary)
% dlmwrite('cylinder/U0.txt', U)
% dlmwrite('cylinder/V0.txt', V)
% dlmwrite('cylinder/P0.txt', P)
% dlmwrite('cylinder/p0Boundary.txt', pBoundary)

% dlmwrite('cylinder/U1.txt', U)
% dlmwrite('cylinder/V1.txt', V)
% dlmwrite('cylinder/P1.txt', P)
% dlmwrite('cylinder/p1Boundary.txt', pBoundary)

% dlmwrite('cylinder/U_CS.txt', imag(U) / -0.0001)
% dlmwrite('cylinder/V_CS.txt', imag(V) / -0.0001)
% dlmwrite('cylinder/P_CS.txt', imag(P) / -0.0001)
% dlmwrite('cylinder/pBoundary_CS.txt', imag(pBoundary) / -0.0001)

