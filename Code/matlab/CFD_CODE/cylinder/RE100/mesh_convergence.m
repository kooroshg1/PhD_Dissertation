clc;
% clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
fontsize = 52;
linewidth = 5.0;
% ----------------------------------------------------------------------- %
% xImm = dlmread('100X75/xImm');
% yImm = dlmread('100X75/yImm');
% 
% U1 = dlmread('100X75/U.txt');
% V1 = dlmread('100X75/V.txt');
% P1 = dlmread('100X75/P.txt');
% uB1 = dlmread('100X75/uBoundary.txt');
% vB1 = dlmread('100X75/vBoundary.txt');
% pB1 = dlmread('100X75/pBoundary.txt');
% Xu1 = dlmread('100X75/Xu.txt');
% Yu1 = dlmread('100X75/Yu.txt');
% Xv1 = dlmread('100X75/Xv.txt');
% Yv1 = dlmread('100X75/Yv.txt');
% Xp1 = dlmread('100X75/Xp.txt');
% Yp1 = dlmread('100X75/Yp.txt');
% uB1 = dlmread('100X75/uBoundary.txt');
% vB1 = dlmread('100X75/vBoundary.txt');
% pB1 = dlmread('100X75/pBoundary.txt');
% 
% [in, on] = inpolygon(Xu1, Yu1, xImm, yImm); U1(in) = nan;
% [in, on] = inpolygon(Xv1, Yv1, xImm, yImm); V1(in) = nan;
% [in, on] = inpolygon(Xp1, Yp1, xImm, yImm); P1(in) = nan;
% 
% U2 = dlmread('200X150/U.txt');
% V2 = dlmread('200X150/V.txt');
% P2 = dlmread('200X150/P.txt');
% uB2 = dlmread('200X150/uBoundary.txt');
% vB2 = dlmread('200X150/vBoundary.txt');
% pB2 = dlmread('200X150/pBoundary.txt');
% Xu2 = dlmread('200X150/Xu.txt');
% Yu2 = dlmread('200X150/Yu.txt');
% Xv2 = dlmread('200X150/Xv.txt');
% Yv2 = dlmread('200X150/Yv.txt');
% Xp2 = dlmread('200X150/Xp.txt');
% Yp2 = dlmread('200X150/Yp.txt');
% uB2 = dlmread('200X150/uBoundary.txt');
% vB2 = dlmread('200X150/vBoundary.txt');
% pB2 = dlmread('200X150/pBoundary.txt');
% 
% [in, on] = inpolygon(Xu2, Yu2, xImm, yImm); U2(in) = nan;
% [in, on] = inpolygon(Xv2, Yv2, xImm, yImm); V2(in) = nan;
% [in, on] = inpolygon(Xp2, Yp2, xImm, yImm); P2(in) = nan;
% 
% U3 = dlmread('300X225/U.txt');
% V3 = dlmread('300X225/V.txt');
% P3 = dlmread('300X225/P.txt');
% uB3 = dlmread('300X225/uBoundary.txt');
% vB3 = dlmread('300X225/vBoundary.txt');
% pB3 = dlmread('300X225/pBoundary.txt');
% Xu3 = dlmread('300X225/Xu.txt');
% Yu3 = dlmread('300X225/Yu.txt');
% Xv3 = dlmread('300X225/Xv.txt');
% Yv3 = dlmread('300X225/Yv.txt');
% Xp3 = dlmread('300X225/Xp.txt');
% Yp3 = dlmread('300X225/Yp.txt');
% uB3 = dlmread('300X225/uBoundary.txt');
% vB3 = dlmread('300X225/vBoundary.txt');
% pB3 = dlmread('300X225/pBoundary.txt');
% 
% [in, on] = inpolygon(Xu3, Yu3, xImm, yImm); U3(in) = nan;
% [in, on] = inpolygon(Xv3, Yv3, xImm, yImm); V3(in) = nan;
% [in, on] = inpolygon(Xp3, Yp3, xImm, yImm); P3(in) = nan;
% 
% U4 = dlmread('400X300/U.txt');
% V4 = dlmread('400X300/V.txt');
% P4 = dlmread('400X300/P.txt');
% uB4 = dlmread('400X300/uBoundary.txt');
% vB4 = dlmread('400X300/vBoundary.txt');
% pB4 = dlmread('400X300/pBoundary.txt');
% Xu4 = dlmread('400X300/Xu.txt');
% Yu4 = dlmread('400X300/Yu.txt');
% Xv4 = dlmread('400X300/Xv.txt');
% Yv4 = dlmread('400X300/Yv.txt');
% Xp4 = dlmread('400X300/Xp.txt');
% Yp4 = dlmread('400X300/Yp.txt');
% uB4 = dlmread('400X300/uBoundary.txt');
% vB4 = dlmread('400X300/vBoundary.txt');
% pB4 = dlmread('400X300/pBoundary.txt');
% 
% [in, on] = inpolygon(Xu4, Yu4, xImm, yImm); U4(in) = nan;
% [in, on] = inpolygon(Xv4, Yv4, xImm, yImm); V4(in) = nan;
% [in, on] = inpolygon(Xp4, Yp4, xImm, yImm); P4(in) = nan;

theta = [linspace(90, -90, length(xImm) / 2), ...
         linspace(-90, 90, length(xImm) / 2)];

% ======================================================================= %
% ============================ U - Velocity ============================= %
% ======================================================================= %
%% Velocity on the boundary
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(theta, uB4, 'k', ...
%      'linewidth', linewidth)
% xlabel('\theta', 'fontsize', fontsize)
% ylabel('u-velocity (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
% grid('on')
% grid minor
% set(gca,'LooseInset',get(gca,'TightInset'))

%% Contour plots
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% contourf(Xu4, Yu4, U4, 50, 'linestyle', 'none')
% colorbar('eastoutside')
% hold on
% plot(xImm, yImm, 'k', ...
%      'linewidth', linewidth)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence plot
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Xu1(:, 38), U1(:, 38), ...
%      Xu2(:, 75), U2(:, 75), ...
%      Xu3(:, 112), U3(:, 112), ...
%      Xu4(:, 150), U4(:, 150), ...
%      'linewidth', linewidth)
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('U (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% legend('100X75', '200X150', '300X225', '400X300', 'location', 'best')
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
x = 2;
[a, b] = find(Xu1(:, 38) < x); U1ind = U1(a(end) + 1, 38);
[a, b] = find(Xu2(:, 75) < x); U2ind = U2(a(end) + 1, 75);
[a, b] = find(Xu3(:, 112) < x); U3ind = U3(a(end) + 1, 112);
[a, b] = find(Xu4(:, 150) < x); U4ind = U4(a(end) + 1, 150);

err = abs([U1ind, U2ind, U3ind, U4ind] - U4ind) ./ U4ind;
meshSize = [100, 200, 300, 400];
f = fit(meshSize',err','power1')
figure,
set(gcf, 'Position', get(0,'Screensize'));
loglog(meshSize, err, 'k', ...
    meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
    'linewidth', linewidth)
xlabel('Number of nodes in X direction', 'fontsize', fontsize)
ylabel('Absolute error', 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'renderer','painters')
grid('on')
grid minor
set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

% ======================================================================= %
% ============================ V - Velocity ============================= %
% ======================================================================= %
%% Velocity on the boundary
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(theta, vB4, 'k', ...
%      'linewidth', linewidth)
% xlabel('\theta', 'fontsize', fontsize)
% ylabel('v-velocity (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
% grid('on')
% grid minor
% set(gca,'LooseInset',get(gca,'TightInset'))

%% Contour plots
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% contourf(Xv4, Yv4, V4, 50, 'linestyle', 'none')
% colorbar('eastoutside')
% hold on
% plot(xImm, yImm, 'k', ...
%      'linewidth', linewidth)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence plot
% x = 0.5;
% [a, b] = find(Xv1(:, 1) < 0.5); y1ind = a(end) + 1;
% [a, b] = find(Xv2(:, 1) < 0.5); y2ind = a(end) + 1;
% [a, b] = find(Xv3(:, 1) < 0.5); y3ind = a(end) + 1;
% [a, b] = find(Xv4(:, 1) < 0.5); y4ind = a(end) + 1;
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Yv1(y1ind, :), V1(y1ind, :), ...
%      Yv2(y2ind, :), V2(y2ind, :), ...
%      Yv3(y3ind, :), V3(y3ind, :), ...
%      Yv4(y4ind, :), V4(y4ind, :), ...
%      'linewidth', linewidth)
% xlabel('Y (m)', 'fontsize', fontsize)
% ylabel('V (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% legend('100X75', '200X150', '300X225', '400X300')
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
x = 0.5;
y = 0.75;

[a, b] = find(Xv1(:, 1) < x); x1ind = a(end) + 1;
[a, b] = find(abs(Yv1(1, :)) < y); y1ind = b(end) + 1;
[Xv1(x1ind, y1ind) Yv1(x1ind, y1ind)];
V1ind = V1(x1ind, y1ind);

[a, b] = find(Xv2(:, 1) < x); x2ind = a(end) + 1;
[a, b] = find(abs(Yv2(1, :)) < y); y2ind = b(end) + 1;
[Xv2(x2ind, y2ind) Yv2(x2ind, y2ind)];
V2ind = V2(x2ind, y2ind);

[a, b] = find(Xv3(:, 1) < x); x3ind = a(end) + 1;
[a, b] = find(abs(Yv3(1, :)) < y); y3ind = b(end) + 1;
[Xv3(x3ind, y3ind) Yv3(x3ind, y3ind)];
V3ind = V3(x3ind, y3ind);

[a, b] = find(Xv4(:, 1) < x); x4ind = a(end) + 1;
[a, b] = find(abs(Yv4(1, :)) < y); y4ind = b(end) + 1;
[Xv4(x4ind, y4ind) Yv4(x4ind, y4ind)];
V4ind = V4(x4ind, y4ind);

err = abs([V1ind, V2ind, V3ind, V4ind] - V4ind) ./ V4ind;
meshSize = [100, 200, 300, 400];
f = fit(meshSize',err','power1')
figure,
set(gcf, 'Position', get(0,'Screensize'));
loglog(meshSize, err, 'k', ...
    meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
    'linewidth', linewidth)
xlabel('Number of nodes in X direction', 'fontsize', fontsize)
ylabel('Absolute error', 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'renderer','painters')
grid('on')
grid minor
set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

% ======================================================================= %
% ============================== Pressure =============================== %
% ======================================================================= %
%% Velocity on the boundary
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(theta, pB4, 'k', ...
%      'linewidth', linewidth)
% xlabel('\theta', 'fontsize', fontsize)
% ylabel('Pressure (Pa)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
% grid('on')
% grid minor
% set(gca,'LooseInset',get(gca,'TightInset'))

%% Contour plots
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% contourf(Xp4, Yp4, P4, 50, 'linestyle', 'none')
% colorbar('eastoutside')
% hold on
% plot(xImm, yImm, 'k', ...
%      'linewidth', linewidth)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence plot
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Xp1(:, 38), P1(:, 38), ...
%      Xp2(:, 75), P2(:, 75), ...
%      Xp3(:, 112), P3(:, 112), ...
%      Xp4(:, 150), P4(:, 150), ...
%      'linewidth', linewidth)
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('P (Pa)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% legend('100X75', '200X150', '300X225', '400X300')
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
x = 2;
[a, b] = find(Xp1(:, 38) < x); P1ind = P1(a(end) + 1, 38);
[a, b] = find(Xp2(:, 75) < x); P2ind = P2(a(end) + 1, 75);
[a, b] = find(Xp3(:, 112) < x); P3ind = P3(a(end) + 1, 112);
[a, b] = find(Xp4(:, 150) < x); P4ind = P4(a(end) + 1, 150);

err = abs([P1ind, P2ind, P3ind, P4ind] - P4ind) ./ abs(P4ind);
meshSize = [100, 200, 300, 400];
f = fit(meshSize',err','power1')
figure,
set(gcf, 'Position', get(0,'Screensize'));
loglog(meshSize, err, 'k', ...
    meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
    'linewidth', linewidth)
xlabel('Number of nodes in X direction', 'fontsize', fontsize)
ylabel('Absolute error', 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'renderer','painters')
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
grid('on')
grid minor
set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.