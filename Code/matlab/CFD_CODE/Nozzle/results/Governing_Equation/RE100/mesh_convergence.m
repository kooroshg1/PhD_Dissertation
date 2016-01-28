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
xImm = dlmread('200X50/xImm.txt');
yImm = dlmread('200X50/yImm.txt');

% -------------------------------- 200X50 ------------------------------- %
U1 = dlmread('200X50/U.txt');
V1 = dlmread('200X50/V.txt');
P1 = dlmread('200X50/P.txt');
Xu1 = dlmread('200X50/Xu.txt');
Yu1 = dlmread('200X50/Yu.txt');
Xv1 = dlmread('200X50/Xv.txt');
Yv1 = dlmread('200X50/Yv.txt');
Xp1 = dlmread('200X50/Xp.txt');
Yp1 = dlmread('200X50/Yp.txt');
u1 = dlmread('200X50/uBoundary.txt');
v1 = dlmread('200X50/vBoundary.txt');
p1 = dlmread('200X50/pBoundary.txt');

[in, on] = inpolygon(Xu1, Yu1, xImm, yImm); U1(in) = nan;
[in, on] = inpolygon(Xv1, Yv1, xImm, yImm); V1(in) = nan;
[in, on] = inpolygon(Xp1, Yp1, xImm, yImm); P1(in) = nan;

% ------------------------------- 300X75 -------------------------------- %
U2 = dlmread('300X75/U.txt');
V2 = dlmread('300X75/V.txt');
P2 = dlmread('300X75/P.txt');
Xu2 = dlmread('300X75/Xu.txt');
Yu2 = dlmread('300X75/Yu.txt');
Xv2 = dlmread('300X75/Xv.txt');
Yv2 = dlmread('300X75/Yv.txt');
Xp2 = dlmread('300X75/Xp.txt');
Yp2 = dlmread('300X75/Yp.txt');
u2 = dlmread('300X75/uBoundary.txt');
v2 = dlmread('300X75/vBoundary.txt');
p2 = dlmread('300X75/pBoundary.txt');

[in, on] = inpolygon(Xu2, Yu2, xImm, yImm); U2(in) = nan;
[in, on] = inpolygon(Xv2, Yv2, xImm, yImm); V2(in) = nan;
[in, on] = inpolygon(Xp2, Yp2, xImm, yImm); P2(in) = nan;

% ------------------------------- 400X100 ------------------------------- %
U3 = dlmread('400X100/U.txt');
V3 = dlmread('400X100/V.txt');
P3 = dlmread('400X100/P.txt');
Xu3 = dlmread('400X100/Xu.txt');
Yu3 = dlmread('400X100/Yu.txt');
Xv3 = dlmread('400X100/Xv.txt');
Yv3 = dlmread('400X100/Yv.txt');
Xp3 = dlmread('400X100/Xp.txt');
Yp3 = dlmread('400X100/Yp.txt');
u3 = dlmread('400X100/uBoundary.txt');
v3 = dlmread('400X100/vBoundary.txt');
p3 = dlmread('400X100/pBoundary.txt');

[in, on] = inpolygon(Xu3, Yu3, xImm, yImm); U3(in) = nan;
[in, on] = inpolygon(Xv3, Yv3, xImm, yImm); V3(in) = nan;
[in, on] = inpolygon(Xp3, Yp3, xImm, yImm); P3(in) = nan;

% ------------------------------- 450X112 ------------------------------- %
U31 = dlmread('450X112/U.txt');
V31 = dlmread('450X112/V.txt');
P31 = dlmread('450X112/P.txt');
Xu31 = dlmread('450X112/Xu.txt');
Yu31 = dlmread('450X112/Yu.txt');
Xv31 = dlmread('450X112/Xv.txt');
Yv31 = dlmread('450X112/Yv.txt');
Xp31 = dlmread('450X112/Xp.txt');
Yp31 = dlmread('450X112/Yp.txt');
u31 = dlmread('450X112/uBoundary.txt');
v31 = dlmread('450X112/vBoundary.txt');
p31 = dlmread('450X112/pBoundary.txt');

[in, on] = inpolygon(Xu31, Yu31, xImm, yImm); U31(in) = nan;
[in, on] = inpolygon(Xv31, Yv31, xImm, yImm); V31(in) = nan;
[in, on] = inpolygon(Xp31, Yp31, xImm, yImm); P31(in) = nan;

% ------------------------------- 500X125 ------------------------------- %
U4 = dlmread('500X125/U.txt');
V4 = dlmread('500X125/V.txt');
P4 = dlmread('500X125/P.txt');
Xu4 = dlmread('500X125/Xu.txt');
Yu4 = dlmread('500X125/Yu.txt');
Xv4 = dlmread('500X125/Xv.txt');
Yv4 = dlmread('500X125/Yv.txt');
Xp4 = dlmread('500X125/Xp.txt');
Yp4 = dlmread('500X125/Yp.txt');
u4 = dlmread('500X125/uBoundary.txt');
v4 = dlmread('500X125/vBoundary.txt');
p4 = dlmread('500X125/pBoundary.txt');

[in, on] = inpolygon(Xu4, Yu4, xImm, yImm); U4(in) = nan;
[in, on] = inpolygon(Xv4, Yv4, xImm, yImm); V4(in) = nan;
[in, on] = inpolygon(Xp4, Yp4, xImm, yImm); P4(in) = nan;

% ------------------------------- 600X150 ------------------------------- %
U5 = dlmread('600X150/U.txt');
V5 = dlmread('600X150/V.txt');
P5 = dlmread('600X150/P.txt');
Xu5 = dlmread('600X150/Xu.txt');
Yu5 = dlmread('600X150/Yu.txt');
Xv5 = dlmread('600X150/Xv.txt');
Yv5 = dlmread('600X150/Yv.txt');
Xp5 = dlmread('600X150/Xp.txt');
Yp5 = dlmread('600X150/Yp.txt');
u5 = dlmread('600X150/uBoundary.txt');
v5 = dlmread('600X150/vBoundary.txt');
p5 = dlmread('600X150/pBoundary.txt');

[in, on] = inpolygon(Xu5, Yu5, xImm, yImm); U5(in) = nan;
[in, on] = inpolygon(Xv5, Yv5, xImm, yImm); V5(in) = nan;
[in, on] = inpolygon(Xp5, Yp5, xImm, yImm); P5(in) = nan;

end
% ======================================================================= %
% ============================ U - Velocity ============================= %
% ======================================================================= %
%% Complex step validation
% y = 0.0;
% [a, b] = min(abs(Yu5(1, :) - y));
% 
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Xu5(:, b), Up5(:, b), 'k',...
%      Xu5(:, b), dUdB5(:, b), 'r--',...
%      'linewidth', linewidth, 'markersize', markersize)
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('u-velocity (m/s)', 'fontsize', fontsize)
% legend('CSA', 'CS', 'fontsize', fontsize)
% title(['RMSE = ' num2str(calcRMSE(dUdB5(:, b), Up5(:, b)))], 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
% grid('on')
% grid minor
% set(gca,'LooseInset',get(gca,'TightInset'))

%% Contour plots
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% contourf(Xu8, Yu8, Up8, 50, 'linestyle', 'none')
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
% plot(Xu1(:, 25), U1(:, 25), ...
%      Xu2(:, 37), U2(:, 37), ...
%      Xu3(:, 50), U3(:, 50), ...
%      Xu31(:, 56), U31(:, 56), ...
%      Xu4(:, 62), U4(:, 62), ...
%      'linewidth', linewidth)
% legend('200X50', '300X75', '400X100', '500X125', '600X150', ...
%        'location', 'best')
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('U (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
x = 0.3;
[a, b] = find(Xu1(:, 25) < x); U1ind = U1(a(end) + 1, 25);
[a, b] = find(Xu2(:, 37) < x); U2ind = U2(a(end) + 1, 37);
[a, b] = find(Xu3(:, 50) < x); U3ind = U3(a(end) + 1, 50);
[a, b] = find(Xu31(:, 56) < x); U31ind = U31(a(end) + 1, 56);
[a, b] = find(Xu4(:, 62) < x); U4ind = U4(a(end) + 1, 62);

err = abs([U1ind, U2ind, U3ind, U31ind, U4ind] - U4ind) ./ abs(U4ind);

meshSize = [200*50, 300*75, 400*100, 500*125, 600*150];
f = fit(meshSize',err','power1')

figure,
set(gcf, 'Position', get(0,'Screensize'));
loglog(meshSize, err, 'k', ...
       meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
       'linewidth', linewidth)
xlabel('Number of nodes', 'fontsize', fontsize)
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
%% Complex step validation
% x = 0.75;
% [a, b] = min(abs(Xv5(:, 1) - x));
% 
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Yv5(b, :), Vp5(b, :), 'k',...
%      Yv5(b, :), dVdB5(b, :), 'r--',...
%      'linewidth', linewidth)
% xlabel('Y (m)', 'fontsize', fontsize)
% ylabel('dV/dr (1/s)', 'fontsize', fontsize)
% legend('CSA', 'CS')
% title(['RMSE = ' num2str(calcRMSE(dVdB5(b, :), Vp5(b, :)))], 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
% grid('on')
% grid minor
% set(gca,'LooseInset',get(gca,'TightInset'))

%% Contour plots
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% contourf(Xv8, Yv8, Vp8, 20, 'linestyle', 'none')
% colorbar('eastoutside')
% hold on
% plot(xImm, yImm, 'k', ...
%      'linewidth', linewidth)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence plot
% x = 0.3;
% [a, b] = find(Xv1(:, 1) < x); y1ind = a(end) + 1;
% [a, b] = find(Xv2(:, 1) < x); y2ind = a(end) + 1;
% [a, b] = find(Xv3(:, 1) < x); y3ind = a(end) + 1;
% [a, b] = find(Xv31(:, 1) < x); y31ind = a(end) + 1;
% [a, b] = find(Xv4(:, 1) < x); y4ind = a(end) + 1;
% 
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Yv1(y1ind, :), V1(y1ind, :), ...
%      Yv31(y31ind, :), V31(y31ind, :), ...
%      Yv3(y3ind, :), V3(y3ind, :), ...
%      Yv2(y2ind, :), V2(y2ind, :), ...
%      Yv4(y4ind, :), V4(y4ind, :), ...
%      'linewidth', linewidth)
% legend('200X50', '300X75', '400X100', '500X125', '600X150', ...
%        'location', 'best')
% xlabel('Y (m)', 'fontsize', fontsize)
% ylabel('V (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.


%% Convergence order
x = 0.3;
y = 0.55;

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

[a, b] = find(Xv31(:, 1) < x); x31ind = a(end) + 1;
[a, b] = find(abs(Yv31(1, :)) < y); y31ind = b(end) + 1;
[Xv31(x31ind, y31ind) Yv31(x31ind, y31ind)]
V31ind = V31(x31ind, y31ind);

[a, b] = find(Xv4(:, 1) < x); x4ind = a(end) + 1;
[a, b] = find(abs(Yv4(1, :)) < y); y4ind = b(end) + 1;
[Xv4(x4ind, y4ind) Yv4(x4ind, y4ind)];
V4ind = V4(x4ind, y4ind);

err = abs([V1ind, V2ind, V3ind, V31ind, V4ind] - V4ind) ./ abs(V4ind);
err = sort(err, 'descend');
meshSize = [200*50, 300*75, 400*100, 500*125, 600*150];
f = fit(meshSize',err','power1')

figure,
set(gcf, 'Position', get(0,'Screensize'));
loglog(meshSize, err, 'k', ...
    meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
    'linewidth', linewidth)
xlabel('Number of nodes', 'fontsize', fontsize)
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
%% Complex step validation
% y = 0.2;
% [a, b] = min(abs(Yu5(1, :) - y));
% 
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Xp5(:, b), Pp5(:, b), 'k',...
%      Xp5(:, b), dPdB5(:, b), 'r--',...
%      'linewidth', linewidth)
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('dP/dr (Pa/m)', 'fontsize', fontsize)
% legend('CSA', 'CS')
% title(['RMSE = ' num2str(calcRMSE(dPdB5(:, b), Pp5(:, b)))], 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
% grid('on')
% grid minor
% set(gca,'LooseInset',get(gca,'TightInset'))

%% Contour plots
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% contourf(Xp8, Yp8, Pp8 / 10, 35, 'linestyle', 'none')
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
% plot(Xp1(:, 25), P1(:, 25), ...
%      Xp2(:, 37), P2(:, 37), ...
%      Xp3(:, 50), P3(:, 50), ...
%      Xp31(:, 56), P31(:, 56), ...
%      Xp4(:, 62), P4(:, 62), ...
%      'linewidth', linewidth)
% legend('200X50', '300X75', '400X100', '500X125', '600X150', ...
%        'location', 'best')
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('P (Pa)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
x = 0.3;
[a, b] = find(Xp1(:, 25) < x); P1ind = P1(a(end) + 1, 25);
[a, b] = find(Xp2(:, 37) < x); P2ind = P2(a(end) + 1, 37);
[a, b] = find(Xp3(:, 50) < x); P3ind = P2(a(end) + 1, 50);
[a, b] = find(Xp31(:, 56) < x); P31ind = P31(a(end) + 1, 56);
[a, b] = find(Xp4(:, 62) < x); P4ind = P4(a(end) + 1, 62);


err = abs([P1ind, P2ind, P3ind, P31ind, P4ind] - P4ind) ./ abs(P4ind);
err = sort(err, 'descend');
meshSize = [200*50, 300*75, 400*100, 500*125, 600*150];
f = fit(meshSize',err','power1')
figure,
set(gcf, 'Position', get(0,'Screensize'));
loglog(meshSize, err, 'k', ...
    meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
    'linewidth', linewidth)
xlabel('Number of nodes', 'fontsize', fontsize)
ylabel('Absolute error', 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
set(gca,'LooseInset',get(gca,'TightInset'))
set(gcf,'renderer','painters')
set(gca,'YMinorTick','on')
set(gca,'XMinorTick','on')
grid('on')
grid minor
set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Pressure sensitivity on surface
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(theta, sgolayfilt(pp5, 1, 15), 'k', ...
%      theta, sgolayfilt(pp6, 1, 15) - 41, 'r', ...
%      'linewidth', linewidth)
% legend('CSA', 'CS')
% xlabel('\theta', 'fontsize', fontsize)
% ylabel('dP/dr (Pa/m)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% set(gca,'YMinorTick','on')
% set(gca,'XMinorTick','on')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
