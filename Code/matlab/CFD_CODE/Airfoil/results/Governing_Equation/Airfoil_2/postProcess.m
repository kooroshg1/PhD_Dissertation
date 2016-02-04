clc;
% clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
fontsize = 52;
linewidth = 5.0;
markersize = 15.0;
skip = 10;
% ----------------------------------------------------------------------- %
%% Load data
loadData = false;
if loadData
clear all;
fontsize = 52;
linewidth = 5.0;
markersize = 15.0;
skip = 10;
xImm = dlmread('400X100_0/xImm.txt');
yImm = dlmread('400X100_0/yImm.txt');

% -------------------------------- 400X100_0 ------------------------------- %
U1 = dlmread('400X100_0/U.txt');
V1 = dlmread('400X100_0/V.txt');
P1 = dlmread('400X100_0/P.txt');
Xu1 = dlmread('400X100_0/Xu.txt');
Yu1 = dlmread('400X100_0/Yu.txt');
Xv1 = dlmread('400X100_0/Xv.txt');
Yv1 = dlmread('400X100_0/Yv.txt');
Xp1 = dlmread('400X100_0/Xp.txt');
Yp1 = dlmread('400X100_0/Yp.txt');
u1 = dlmread('400X100_0/uBoundary.txt');
v1 = dlmread('400X100_0/vBoundary.txt');
p1 = dlmread('400X100_0/pBoundary.txt');

[in, on] = inpolygon(Xu1, Yu1, xImm, yImm); U1(in) = nan;
[in, on] = inpolygon(Xv1, Yv1, xImm, yImm); V1(in) = nan;
[in, on] = inpolygon(Xp1, Yp1, xImm, yImm); P1(in) = nan;

% ------------------------------- 400X100_0001 ------------------------------- %
U2 = dlmread('400X100_0001/U.txt');
V2 = dlmread('400X100_0001/V.txt');
P2 = dlmread('400X100_0001/P.txt');
Xu2 = dlmread('400X100_0001/Xu.txt');
Yu2 = dlmread('400X100_0001/Yu.txt');
Xv2 = dlmread('400X100_0001/Xv.txt');
Yv2 = dlmread('400X100_0001/Yv.txt');
Xp2 = dlmread('400X100_0001/Xp.txt');
Yp2 = dlmread('400X100_0001/Yp.txt');
u2 = dlmread('400X100_0001/uBoundary.txt');
v2 = dlmread('400X100_0001/vBoundary.txt');
p2 = dlmread('400X100_0001/pBoundary.txt');

[in, on] = inpolygon(Xu2, Yu2, xImm, yImm); U2(in) = nan;
[in, on] = inpolygon(Xv2, Yv2, xImm, yImm); V2(in) = nan;
[in, on] = inpolygon(Xp2, Yp2, xImm, yImm); P2(in) = nan;

end
% ======================================================================= %
% ============================ U - Velocity ============================= %
% ======================================================================= %
%% Sensitivity Calculation
dy = 0.001;
dUdC1 = (U2 - U1) / dy; dudC1 = (u2 - u1) / dy;
dVdC1 = (V2 - V1) / dy; dvdC1 = (v2 - v1) / dy;
dPdC1 = (P2 - P1) / dy; dpdC1 = (p2 - p1) / dy;

figure,
contourf(Xu1, Yu1, U1, 50, 'linestyle', 'none')
colorbar()
hold on
plot(xImm, yImm, 'k', 'linewidth', linewidth)
axis equal;
set(gcf,'renderer','painters')
set(gca, 'fontsize', fontsize)

figure,
contourf(Xp1, Yp1, P1, 50, 'linestyle', 'none')
colorbar()
hold on
plot(xImm, yImm, 'k', 'linewidth', linewidth)
axis equal;
set(gcf,'renderer','painters')
set(gca, 'fontsize', fontsize)


% figure,
% contourf(Xu1, Yu1, dUdC1, 50, 'linestyle', 'none')
% colorbar()
% hold on
% plot(xImm, yImm, 'k', 'linewidth', linewidth)
% axis equal;
% set(gcf,'renderer','painters')
% set(gca, 'fontsize', fontsize)
% 
% figure,
% contourf(Xv1, Yv1, dVdC1, 40, 'linestyle', 'none')
% colorbar()
% hold on
% plot(xImm, yImm, 'k', 'linewidth', linewidth)
% axis equal;
% set(gcf,'renderer','painters')
% set(gca, 'fontsize', fontsize)
% 
% figure,
% contourf(Xp1, Yp1, dPdC1, 50, 'linestyle', 'none')
% colorbar()
% hold on
% plot(xImm, yImm, 'k', 'linewidth', linewidth)
% axis equal;
% set(gcf,'renderer','painters')
% set(gca, 'fontsize', fontsize)

% figure,
% plot(xImm, sgolayfilt(dpdC1, 9, 61), 'k', ...
%      xImm, sgolayfilt(dpdC1, 9, 61), 'r--', ...
%      'linewidth', linewidth)
% legend('CSA', 'CS')
% xlabel('Chord', 'fontsize', fontsize);
% ylabel('dP/dC', 'fontsize', fontsize);
% set(gcf,'renderer','painters')
% set(gca, 'fontsize', fontsize)

 
 
