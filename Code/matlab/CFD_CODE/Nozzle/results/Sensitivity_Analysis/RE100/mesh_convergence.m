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
xImm = dlmread('200X50/xImm.txt');
yImm = dlmread('200X50/yImm.txt');

% ------------------------------- 200X50 -------------------------------- %
dUdB1 = dlmread('200X50/dUdB.txt');
dVdB1 = dlmread('200X50/dVdB.txt');
dPdB1 = dlmread('200X50/dPdB.txt');
Xu1 = dlmread('200X50/Xu.txt');
Yu1 = dlmread('200X50/Yu.txt');
Xv1 = dlmread('200X50/Xv.txt');
Yv1 = dlmread('200X50/Yv.txt');
Xp1 = dlmread('200X50/Xp.txt');
Yp1 = dlmread('200X50/Yp.txt');
dudB1 = dlmread('200X50/dudB.txt');
dvdB1 = dlmread('200X50/dvdB.txt');
dpdB1 = dlmread('200X50/dpdB.txt');

Up1 = dlmread('200X50_CSA/Up.txt');
Vp1 = dlmread('200X50_CSA/Vp.txt');
Pp1 = dlmread('200X50_CSA/Pp.txt');
pp1 = dlmread('200X50_CSA/ppB.txt');

[in, on] = inpolygon(Xu1, Yu1, xImm, yImm); dUdB1(in) = nan; Up1(in) = nan;
[in, on] = inpolygon(Xv1, Yv1, xImm, yImm); dVdB1(in) = nan; Vp1(in) = nan;
[in, on] = inpolygon(Xp1, Yp1, xImm, yImm); dPdB1(in) = nan; Pp1(in) = nan;

% ------------------------------- 300X75 -------------------------------- %
dUdB2 = dlmread('300X75/dUdB.txt');
dVdB2 = dlmread('300X75/dVdB.txt');
dPdB2 = dlmread('300X75/dPdB.txt');
Xu2 = dlmread('300X75/Xu.txt');
Yu2 = dlmread('300X75/Yu.txt');
Xv2 = dlmread('300X75/Xv.txt');
Yv2 = dlmread('300X75/Yv.txt');
Xp2 = dlmread('300X75/Xp.txt');
Yp2 = dlmread('300X75/Yp.txt');
dudB2 = dlmread('300X75/dudB.txt');
dvdB2 = dlmread('300X75/dvdB.txt');
dpdB2 = dlmread('300X75/dpdB.txt');

Up2 = dlmread('300X75_CSA/Up.txt');
Vp2 = dlmread('300X75_CSA/Vp.txt');
Pp2 = dlmread('300X75_CSA/Pp.txt');
pp2 = dlmread('300X75_CSA/ppB.txt');

[in, on] = inpolygon(Xu2, Yu2, xImm, yImm); dUdB2(in) = nan; Up2(in) = nan;
[in, on] = inpolygon(Xv2, Yv2, xImm, yImm); dVdB2(in) = nan; Vp2(in) = nan;
[in, on] = inpolygon(Xp2, Yp2, xImm, yImm); dPdB2(in) = nan; Pp2(in) = nan;

% ------------------------------- 400X100 ------------------------------- %
dUdB3 = dlmread('400X100/dUdB.txt');
dVdB3 = dlmread('400X100/dVdB.txt');
dPdB3 = dlmread('400X100/dPdB.txt');
Xu3 = dlmread('400X100/Xu.txt');
Yu3 = dlmread('400X100/Yu.txt');
Xv3 = dlmread('400X100/Xv.txt');
Yv3 = dlmread('400X100/Yv.txt');
Xp3 = dlmread('400X100/Xp.txt');
Yp3 = dlmread('400X100/Yp.txt');
dudB3 = dlmread('400X100/dudB.txt');
dvdB3 = dlmread('400X100/dvdB.txt');
dpdB3 = dlmread('400X100/dpdB.txt');

Up3 = dlmread('400X100_CSA/Up.txt');
Vp3 = dlmread('400X100_CSA/Vp.txt');
Pp3 = dlmread('400X100_CSA/Pp.txt');
pp3 = dlmread('400X100_CSA/ppB.txt');

[in, on] = inpolygon(Xu3, Yu3, xImm, yImm); dUdB3(in) = nan; Up3(in) = nan;
[in, on] = inpolygon(Xv3, Yv3, xImm, yImm); dVdB3(in) = nan; Vp3(in) = nan;
[in, on] = inpolygon(Xp3, Yp3, xImm, yImm); dPdB3(in) = nan; Pp3(in) = nan;

% ------------------------------- 450X112 ------------------------------- %
dUdB4 = dlmread('450X112/dUdB.txt');
dVdB4 = dlmread('450X112/dVdB.txt');
dPdB4 = dlmread('450X112/dPdB.txt');
Xu4 = dlmread('450X112/Xu.txt');
Yu4 = dlmread('450X112/Yu.txt');
Xv4 = dlmread('450X112/Xv.txt');
Yv4 = dlmread('450X112/Yv.txt');
Xp4 = dlmread('450X112/Xp.txt');
Yp4 = dlmread('450X112/Yp.txt');
dudB4 = dlmread('450X112/dudB.txt');
dvdB4 = dlmread('450X112/dvdB.txt');
dpdB4 = dlmread('450X112/dpdB.txt');

Up4 = dlmread('450X112_CSA/Up.txt');
Vp4 = dlmread('450X112_CSA/Vp.txt');
Pp4 = dlmread('450X112_CSA/Pp.txt');
pp4 = dlmread('450X112_CSA/ppB.txt');

[in, on] = inpolygon(Xu4, Yu4, xImm, yImm); dUdB4(in) = nan; Up4(in) = nan;
[in, on] = inpolygon(Xv4, Yv4, xImm, yImm); dVdB4(in) = nan; Vp4(in) = nan;
[in, on] = inpolygon(Xp4, Yp4, xImm, yImm); dPdB4(in) = nan; Pp4(in) = nan;

end
% ======================================================================= %
% ============================ U - Velocity ============================= %
% ======================================================================= %
%% Complex step validation
y = 0.5;
[a, b] = min(abs(Yu3(1, :) - y));

figure,
set(gcf, 'Position', get(0,'Screensize'));
plot(Xu3(:, b), Up3(:, b), 'k',...
     Xu3(:, b), dUdB3(:, b), 'r--',...
     'linewidth', linewidth, 'markersize', markersize)
xlabel('X (m)', 'fontsize', fontsize)
ylabel('u-velocity (m/s)', 'fontsize', fontsize)
legend('CSA', 'CS', 'fontsize', fontsize)
title(['RMSE = ' num2str(calcRMSE(dUdB3(:, b), Up3(:, b)))], 'fontsize', fontsize)
set(gca, 'fontsize', fontsize)
set(gcf,'renderer','painters')
set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
grid('on')
grid minor
set(gca,'LooseInset',get(gca,'TightInset'))

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
% plot(Xu1(:, 25), Up1(:, 25), ...
%      Xu2(:, 37), Up2(:, 37), ...
%      Xu3(:, 50), Up3(:, 50), ...
%      Xu4(:, 56), Up4(:, 56), ...
%      'linewidth', linewidth)
% legend('300X225', '400X300', '500X375', '600X450', ...
%        'location', 'southeast')
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('dU/dr (1/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
% x = 1.5;
% [a, b] = find(Xu3(:, 112) < x); Up3ind = Up3(a(end) + 1, 112);
% [a, b] = find(Xu4(:, 150) < x); Up4ind = Up4(a(end) + 1, 150);
% [a, b] = find(Xu5(:, 187) < x); Up5ind = Up5(a(end) + 1, 187);
% [a, b] = find(Xu6(:, 225) < x); Up6ind = Up6(a(end) + 1, 225);
% 
% err = abs([Up3ind, Up4ind, Up5ind, Up6ind] - Up6ind) ./ abs(Up6ind);
% err = sort(err, 'descend');
% meshSize = [300*225, 400*300, 500*375, 600*450];
% f = fit(meshSize',err','power1')
% 
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% loglog(meshSize, err, 'k', ...
%     meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
%     'linewidth', linewidth)
% xlabel('Number of nodes', 'fontsize', fontsize)
% ylabel('Absolute error', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

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
% x = 0.75;
% [a, b] = find(Xv3(:, 1) < x); y3ind = a(end) + 1;
% [a, b] = find(Xv4(:, 1) < x); y4ind = a(end) + 1;
% [a, b] = find(Xv5(:, 1) < x); y5ind = a(end) + 1;
% [a, b] = find(Xv6(:, 1) < x); y6ind = a(end) + 1;
% [a, b] = find(Xv7(:, 1) < x); y7ind = a(end) + 1;
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(Yv3(y3ind, :), Vp3(y3ind, :), ...
%      Yv4(y4ind, :), Vp4(y4ind, :), ...
%      Yv5(y5ind, :), Vp5(y5ind, :), ...
%      Yv6(y6ind, :), Vp6(y6ind, :), ...
%      'linewidth', linewidth)
% legend('300X225', '400X300', '500X375', '600X450', ...
%        'location', 'southeast')
% xlabel('Y (m)', 'fontsize', fontsize)
% ylabel('dV/dr (1/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
% x = 0.5;
% y = 0.75;
% 
% [a, b] = find(Xv3(:, 1) < x); x3ind = a(end) + 1;
% [a, b] = find(abs(Yv3(1, :)) < y); y3ind = b(end) + 1;
% [Xv3(x3ind, y3ind) Yv3(x3ind, y3ind)];
% Vp3ind = Vp3(x3ind, y3ind);
% 
% [a, b] = find(Xv4(:, 1) < x); x4ind = a(end) + 1;
% [a, b] = find(abs(Yv4(1, :)) < y); y4ind = b(end) + 1;
% [Xv4(x4ind, y4ind) Yv4(x4ind, y4ind)];
% Vp4ind = Vp4(x4ind, y4ind);
% 
% [a, b] = find(Xv5(:, 1) < x); x5ind = a(end) + 1;
% [a, b] = find(abs(Yv5(1, :)) < y); y5ind = b(end) + 1;
% [Xv5(x5ind, y5ind) Yv5(x5ind, y5ind)];
% Vp5ind = Vp5(x5ind, y5ind);
% 
% [a, b] = find(Xv6(:, 1) < x); x6ind = a(end) + 1;
% [a, b] = find(abs(Yv6(1, :)) < y); y6ind = b(end) + 1;
% [Xv6(x6ind, y6ind) Yv6(x6ind, y6ind)];
% Vp6ind = Vp6(x6ind, y6ind);
% 
% err = abs([Vp3ind, Vp4ind, Vp5ind, Vp6ind] - Vp6ind) ./ Vp6ind;
% err = sort(err, 'descend')
% meshSize = [300*225, 400*300, 500*375, 600*450];
% f = fit(meshSize',err','power1')
% 
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% loglog(meshSize, err, 'k', ...
%     meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
%     'linewidth', linewidth)
% xlabel('Number of nodes', 'fontsize', fontsize)
% ylabel('Absolute error', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

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
% plot(Xp3(:, 112), Pp3(:, 112), ...
%      Xp4(:, 150), Pp4(:, 150), ...
%      Xp5(:, 187), Pp5(:, 187), ...
%      Xp6(:, 225), Pp6(:, 225), ...
%      'linewidth', linewidth)
% legend('300X225', '400X300', '500X375', '600X450', ...
%        'location', 'southeast')
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('dP/dr (Pa/m)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence order
% x = 1.0;
% [a, b] = find(Xp3(:, 112) < x); Pp3ind = Pp3(a(end) + 1, 112);
% [a, b] = find(Xp4(:, 150) < x); Pp4ind = Pp4(a(end) + 1, 150);
% [a, b] = find(Xp5(:, 187) < x); Pp5ind = Pp5(a(end) + 1, 187);
% [a, b] = find(Xp6(:, 225) < x); Pp6ind = Pp6(a(end) + 1, 225);
% 
% err = abs([Pp3ind, Pp4ind, Pp5ind, Pp6ind] - Pp6ind) ./ abs(Pp6ind);
% err = sort(err, 'descend')
% meshSize = [300*225, 400*300, 500*375, 600*450];
% f = fit(meshSize',err','power1')
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% loglog(meshSize, err, 'k', ...
%     meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
%     'linewidth', linewidth)
% xlabel('Number of nodes', 'fontsize', fontsize)
% ylabel('Absolute error', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% set(gca,'YMinorTick','on')
% set(gca,'XMinorTick','on')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

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
