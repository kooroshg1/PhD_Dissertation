clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
multiplier = 2;
N1 = 10 * multiplier;
N2 = 50 * multiplier;
N3 = 10 * multiplier;
N4 = 50 * multiplier;
N5 = 20 * multiplier;
% Calculate point cloud
yi = 0.3;
xi = 0.0;
yo = 0.3;
xo = 1.0-xi;
xt = 0.3;
yt = 0.1 + 0.0000075;

A1 = [xt^2,xt,1;2*xt,1,0;xi^2,xi,1];
coef = A1 \ [yt;0;yi];
a1 = coef(1);
b1 = coef(2);
c1 = coef(3);

A2 = [xt^2,xt,1;2*xt,1,0;xo^2,xo,1];
coef = A2 \ [yt;0;yo];
a2 = coef(1);
b2 = coef(2);
c2 = coef(3);

f1 = @(x) a1 * x.^2 + b1 * x + c1 + 0.5;
f2 = @(x) a2 * x.^2 + b2 * x + c2 + 0.5;
f3 = @(x) -a1 * x.^2 - b1 * x - c1 + 0.5;
f4 = @(x) -a2 * x.^2 - b2 * x - c2 + 0.5;

x1 = linspace(0,xt,10)'; x3 = x1;
x2 = linspace(xt,1.0,10)'; x4 = x2;

pointCloud1 = [xi*ones(N1,1) linspace(yi+0.5,1,N1)';...
                   linspace(xi,xo,N2)' ones(N2,1);...
                   xo*ones(N3,1) linspace(1,yo+0.5,N3)';
                   linspace(xo,xt,N4)' f2(linspace(xo,xt,N4)');...
                   linspace(xt,xi,N5)' f1(linspace(xt,xi,N5)')];
pointCloud2 = [xi*ones(N1,1) linspace(0.5-yi,0,N1)';...
                   linspace(xi,xo,N2)' 0*ones(N2,1);...
                   xo*ones(N3,1) linspace(0,0.5-yo,N3)';
                   linspace(xo,xt,N4)' f4(linspace(xo,xt,N4)');...
                   linspace(xt,xi,N5)' f3(linspace(xt,xi,N5)')];

% dlmwrite('topBack.txt',[linspace(xo,xt,N4)' 0.5-f2(linspace(xo,xt,N4)')]);
% dlmwrite('topFront.txt',[linspace(xt,xi,N5)' 0.5-f1(linspace(xt,xi,N5)')]);
% dlmwrite('bottomBack.txt',[linspace(xo,xt,N4)' 0.5-f4(linspace(xo,xt,N4)')]);
% dlmwrite('bottomFront.txt',[linspace(xt,xi,N5)' 0.5-f3(linspace(xt,xi,N5)')]);

findSame = pointCloud1(1:end-1,:) == pointCloud1(2:end,:);
findSame = findSame(:,1) & findSame(:,2);
findSame = [N1,N1+N2,N1+N2+N3,N1+N2+N3+N4];
pointCloud1(findSame,:) = [];

% findSame = pointCloud2(1:end-1,:) == pointCloud2(2:end,:);
% findSame = findSame(:,1) & findSame(:,2);
pointCloud2(findSame,:) = [];

pointCloud = [pointCloud1;pointCloud2];

% Calculate point cloud perturbation
dyt = 0.02;
yt = 0.1 + dyt;


A1 = [xt^2,xt,1;2*xt,1,0;xi^2,xi,1];
coef = A1 \ [yt;0;yi];
a1 = coef(1);
b1 = coef(2);
c1 = coef(3);

A2 = [xt^2,xt,1;2*xt,1,0;xo^2,xo,1];
coef = A2 \ [yt;0;yo];
a2 = coef(1);
b2 = coef(2);
c2 = coef(3);

f1 = @(x) a1 * x.^2 + b1 * x + c1 + 0.5;
f2 = @(x) a2 * x.^2 + b2 * x + c2 + 0.5;
f3 = @(x) -a1 * x.^2 - b1 * x - c1 + 0.5;
f4 = @(x) -a2 * x.^2 - b2 * x - c2 + 0.5;

x1 = linspace(0,xt,10)'; x3 = x1;
x2 = linspace(xt,1.0,10)'; x4 = x2;

pointCloudPerturb1 = [xi*ones(N1,1) linspace(yi+0.5,1,N1)';...
                   linspace(xi,xo,N2)' ones(N2,1);...
                   xo*ones(N3,1) linspace(1,yo+0.5,N3)';
                   linspace(xo,xt,N4)' f2(linspace(xo,xt,N4)');...
                   linspace(xt,xi,N5)' f1(linspace(xt,xi,N5)')];
pointCloudPerturb2 = [xi*ones(N1,1) linspace(0.5-yi,0,N1)';...
                   linspace(xi,xo,N2)' 0*ones(N2,1);...
                   xo*ones(N3,1) linspace(0,0.5-yo,N3)';
                   linspace(xo,xt,N4)' f4(linspace(xo,xt,N4)');...
                   linspace(xt,xi,N5)' f3(linspace(xt,xi,N5)')];


% findSame = pointCloudPerturb1(1:end-1,:) == pointCloudPerturb1(2:end,:);
% findSame = findSame(:,1) & findSame(:,2);
pointCloudPerturb1(findSame,:) = [];

% findSame = pointCloudPerturb2(1:end-1,:) == pointCloudPerturb2(2:end,:);
% findSame = findSame(:,1) & findSame(:,2);
pointCloudPerturb2(findSame,:) = [];


figure,
plot(pointCloud1(:,1),pointCloud1(:,2),'o',...
     pointCloud2(:,1),pointCloud2(:,2),'o')
axis([0 1 0 1])

pointCloudPerturb = [pointCloudPerturb1;pointCloudPerturb2];
pointCloudSensitivity = (pointCloudPerturb - pointCloud) / dyt;

dlmwrite('point_cloud.txt',pointCloud);
dlmwrite('point_cloud_sensitivity.txt',pointCloudSensitivity);

