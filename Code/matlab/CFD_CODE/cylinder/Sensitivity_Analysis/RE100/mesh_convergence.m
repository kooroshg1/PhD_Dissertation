clc;
% clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
fontsize = 52;
linewidth = 5.0;
% ----------------------------------------------------------------------- %
%% Load data
loadData = false;
if loadData
clear all;
xImm = dlmread('100X75/xImm');
yImm = dlmread('100X75/yImm');

% -------------------------------- 100X75 ------------------------------- %
dUdB1 = dlmread('100X75/dUdB.txt');
dVdB1 = dlmread('100X75/dVdB.txt');
dPdB1 = dlmread('100X75/dPdB.txt');
Xu1 = dlmread('100X75/Xu.txt');
Yu1 = dlmread('100X75/Yu.txt');
Xv1 = dlmread('100X75/Xv.txt');
Yv1 = dlmread('100X75/Yv.txt');
Xp1 = dlmread('100X75/Xp.txt');
Yp1 = dlmread('100X75/Yp.txt');
dudB1 = dlmread('100X75/dudB.txt');
dvdB1 = dlmread('100X75/dvdB.txt');
dpdB1 = dlmread('100X75/dpdB.txt');

[in, on] = inpolygon(Xu1, Yu1, xImm, yImm); dUdB1(in) = nan;
[in, on] = inpolygon(Xv1, Yv1, xImm, yImm); dVdB1(in) = nan;
[in, on] = inpolygon(Xp1, Yp1, xImm, yImm); dPdB1(in) = nan;

% ------------------------------- 200X150 ------------------------------- %
dUdB2 = dlmread('200X150/dUdB.txt');
dVdB2 = dlmread('200X150/dVdB.txt');
dPdB2 = dlmread('200X150/dPdB.txt');
Xu2 = dlmread('200X150/Xu.txt');
Yu2 = dlmread('200X150/Yu.txt');
Xv2 = dlmread('200X150/Xv.txt');
Yv2 = dlmread('200X150/Yv.txt');
Xp2 = dlmread('200X150/Xp.txt');
Yp2 = dlmread('200X150/Yp.txt');
dudB2 = dlmread('200X150/dudB.txt');
dvdB2 = dlmread('200X150/dvdB.txt');
dpdB2 = dlmread('200X150/dpdB.txt');

Up2 = dlmread('200X150_CSA/Up.txt');
Vp2 = dlmread('200X150_CSA/Vp.txt');
Pp2 = dlmread('200X150_CSA/Pp.txt');
pp2 = dlmread('200X150_CSA/ppB.txt');

[in, on] = inpolygon(Xu2, Yu2, xImm, yImm); dUdB2(in) = nan; Up2(in) = nan;
[in, on] = inpolygon(Xv2, Yv2, xImm, yImm); dVdB2(in) = nan; Vp2(in) = nan;
[in, on] = inpolygon(Xp2, Yp2, xImm, yImm); dPdB2(in) = nan; Pp2(in) = nan;

% ------------------------------- 300X225 ------------------------------- %
dUdB3 = dlmread('300X225/dUdB.txt');
dVdB3 = dlmread('300X225/dVdB.txt');
dPdB3 = dlmread('300X225/dPdB.txt');
Xu3 = dlmread('300X225/Xu.txt');
Yu3 = dlmread('300X225/Yu.txt');
Xv3 = dlmread('300X225/Xv.txt');
Yv3 = dlmread('300X225/Yv.txt');
Xp3 = dlmread('300X225/Xp.txt');
Yp3 = dlmread('300X225/Yp.txt');
dudB3 = dlmread('300X225/dudB.txt');
dvdB3 = dlmread('300X225/dvdB.txt');
dpdB3 = dlmread('300X225/dpdB.txt');

Up3 = dlmread('300X225_CSA/Up.txt');
Vp3 = dlmread('300X225_CSA/Vp.txt');
Pp3 = dlmread('300X225_CSA/Pp.txt');
pp3 = dlmread('300X225_CSA/ppB.txt');

[in, on] = inpolygon(Xu3, Yu3, xImm, yImm); dUdB3(in) = nan; Up3(in) = nan;
[in, on] = inpolygon(Xv3, Yv3, xImm, yImm); dVdB3(in) = nan; Vp3(in) = nan;
[in, on] = inpolygon(Xp3, Yp3, xImm, yImm); dPdB3(in) = nan; Pp3(in) = nan;

% ------------------------------- 400X300 ------------------------------- %
dUdB4 = dlmread('400X300/dUdB.txt');
dVdB4 = dlmread('400X300/dVdB.txt');
dPdB4 = dlmread('400X300/dPdB.txt');
Xu4 = dlmread('400X300/Xu.txt');
Yu4 = dlmread('400X300/Yu.txt');
Xv4 = dlmread('400X300/Xv.txt');
Yv4 = dlmread('400X300/Yv.txt');
Xp4 = dlmread('400X300/Xp.txt');
Yp4 = dlmread('400X300/Yp.txt');
dudB4 = dlmread('400X300/dudB.txt');
dvdB4 = dlmread('400X300/dvdB.txt');
dpdB4 = dlmread('400X300/dpdB.txt');

Up4 = dlmread('400X300_CSA/Up.txt');
Vp4 = dlmread('400X300_CSA/Vp.txt');
Pp4 = dlmread('400X300_CSA/Pp.txt');
pp4 = dlmread('400X300_CSA/ppB.txt');

[in, on] = inpolygon(Xu4, Yu4, xImm, yImm); dUdB4(in) = nan; Up4(in) = nan;
[in, on] = inpolygon(Xv4, Yv4, xImm, yImm); dVdB4(in) = nan; Vp4(in) = nan;
[in, on] = inpolygon(Xp4, Yp4, xImm, yImm); dPdB4(in) = nan; Pp4(in) = nan;

% ------------------------------- 500x375 ------------------------------- %
dUdB5 = dlmread('500X375/dUdB.txt');
dVdB5 = dlmread('500X375/dVdB.txt');
dPdB5 = dlmread('500X375/dPdB.txt');
Xu5 = dlmread('500X375/Xu.txt');
Yu5 = dlmread('500X375/Yu.txt');
Xv5 = dlmread('500X375/Xv.txt');
Yv5 = dlmread('500X375/Yv.txt');
Xp5 = dlmread('500X375/Xp.txt');
Yp5 = dlmread('500X375/Yp.txt');
dudB5 = dlmread('500X375/dudB.txt');
dvdB5 = dlmread('500X375/dvdB.txt');
dpdB5 = dlmread('500X375/dpdB.txt');

Up5 = dlmread('500X375_CSA/Up.txt');
Vp5 = dlmread('500X375_CSA/Vp.txt');
Pp5 = dlmread('500X375_CSA/Pp.txt');
pp5 = dlmread('500X375_CSA/ppB.txt');

[in, on] = inpolygon(Xu5, Yu5, xImm, yImm); dUdB5(in) = nan; Up5(in) = nan;
[in, on] = inpolygon(Xv5, Yv5, xImm, yImm); dVdB5(in) = nan; Vp5(in) = nan;
[in, on] = inpolygon(Xp5, Yp5, xImm, yImm); dPdB5(in) = nan; Pp5(in) = nan;

% ------------------------------- 600X450 ------------------------------- %
dUdB6 = dlmread('600X450/dUdB.txt');
dVdB6 = dlmread('600X450/dVdB.txt');
dPdB6 = dlmread('600X450/dPdB.txt');
Xu6 = dlmread('600X450/Xu.txt');
Yu6 = dlmread('600X450/Yu.txt');
Xv6 = dlmread('600X450/Xv.txt');
Yv6 = dlmread('600X450/Yv.txt');
Xp6 = dlmread('600X450/Xp.txt');
Yp6 = dlmread('600X450/Yp.txt');
dudB6 = dlmread('600X450/dudB.txt');
dvdB6 = dlmread('600X450/dvdB.txt');
dpdB6 = dlmread('600X450/dpdB.txt');

Up6 = dlmread('600X450_CSA/Up.txt');
Vp6 = dlmread('600X450_CSA/Vp.txt');
Pp6 = dlmread('600X450_CSA/Pp.txt');
pp6 = dlmread('600X450_CSA/ppB.txt');

[in, on] = inpolygon(Xu6, Yu6, xImm, yImm); dUdB6(in) = nan; Up6(in) = nan;
[in, on] = inpolygon(Xv6, Yv6, xImm, yImm); dVdB6(in) = nan; Vp6(in) = nan;
[in, on] = inpolygon(Xp6, Yp6, xImm, yImm); dPdB6(in) = nan; Pp6(in) = nan;

% ------------------------------- 700X525 ------------------------------- %
dUdB7 = dlmread('700X525/dUdB.txt');
dVdB7 = dlmread('700X525/dVdB.txt');
dPdB7 = dlmread('700X525/dPdB.txt');
Xu7 = dlmread('700X525/Xu.txt');
Yu7 = dlmread('700X525/Yu.txt');
Xv7 = dlmread('700X525/Xv.txt');
Yv7 = dlmread('700X525/Yv.txt');
Xp7 = dlmread('700X525/Xp.txt');
Yp7 = dlmread('700X525/Yp.txt');
dudB7 = dlmread('700X525/dudB.txt');
dvdB7 = dlmread('700X525/dvdB.txt');
dpdB7 = dlmread('700X525/dpdB.txt');

Up7 = dlmread('700X525_CSA/Up.txt');
Vp7 = dlmread('700X525_CSA/Vp.txt');
Pp7 = dlmread('700X525_CSA/Pp.txt');
pp7 = dlmread('700X525_CSA/ppB.txt');

[in, on] = inpolygon(Xu7, Yu7, xImm, yImm); dUdB7(in) = nan; Up7(in) = nan;
[in, on] = inpolygon(Xv7, Yv7, xImm, yImm); dVdB7(in) = nan; Vp7(in) = nan;
[in, on] = inpolygon(Xp7, Yp7, xImm, yImm); dPdB7(in) = nan; Pp7(in) = nan;

% ------------------------------- 800X600 ------------------------------- %
dUdB8 = dlmread('800X600/dUdB.txt');
dVdB8 = dlmread('800X600/dVdB.txt');
dPdB8 = dlmread('800X600/dPdB.txt');
Xu8 = dlmread('800X600/Xu.txt');
Yu8 = dlmread('800X600/Yu.txt');
Xv8 = dlmread('800X600/Xv.txt');
Yv8 = dlmread('800X600/Yv.txt');
Xp8 = dlmread('800X600/Xp.txt');
Yp8 = dlmread('800X600/Yp.txt');
dudB8 = dlmread('800X600/dudB.txt');
dvdB8 = dlmread('800X600/dvdB.txt');
dpdB8 = dlmread('800X600/dpdB.txt');

Up8 = dlmread('800X600_CSA/Up.txt');
Vp8 = dlmread('800X600_CSA/Vp.txt');
Pp8 = dlmread('800X600_CSA/Pp.txt');
pp8 = dlmread('800X600_CSA/ppB.txt');

[in, on] = inpolygon(Xu8, Yu8, xImm, yImm); dUdB8(in) = nan; Up8(in) = nan;
[in, on] = inpolygon(Xv8, Yv8, xImm, yImm); dVdB8(in) = nan; Vp8(in) = nan;
[in, on] = inpolygon(Xp8, Yp8, xImm, yImm); dPdB8(in) = nan; Pp8(in) = nan;

theta = [linspace(90, -90, length(xImm) / 2), ...
         linspace(-90, 90, length(xImm) / 2)];
end
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
% contourf(Xu3, Yu3, abs(dUdB3 - Up3) ./ Up3, 50, 'linestyle', 'none')
% colorbar('eastoutside')
% 
% hold on
% plot(xImm, yImm, 'k', ...
%      'linewidth', linewidth)
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

%% Convergence plot
fontsize = 12;
linewidth = 3;
figure,
set(gcf, 'Position', get(0,'Screensize'));
plot(Xu2(:, 75), dUdB2(:, 75), 'g--',...
     Xu3(:, 112), dUdB3(:, 112), 'y--',...
     Xu4(:, 150), dUdB4(:, 150), 'b--',...
     Xu5(:, 187), dUdB5(:, 187), 'c--',...
     Xu6(:, 225), dUdB6(:, 225), 'k--',...
     Xu7(:, 262), dUdB7(:, 262), 'm--',...
     Xu8(:, 300), dUdB8(:, 300), 'r--',...
     Xu2(:, 75), Up2(:, 75), 'g',...
     Xu3(:, 112), Up3(:, 112), 'y',...
     Xu4(:, 150), Up4(:, 150), 'b',...
     Xu5(:, 187), Up5(:, 187), 'c',...
     Xu6(:, 225), Up6(:, 225), 'k',...
     Xu7(:, 262), Up7(:, 262), 'm',...
     Xu8(:, 300), Up8(:, 300), 'r',...
     'linewidth', linewidth)
legend('CS2', 'CS3', 'CS4', 'CS5', 'CS6', 'CS7', 'CS8', ...
       'CSA2', 'CSA3', 'CSA4', 'CSA5', 'CSA6', 'CSA7', 'CSA8')

% figure,
% plot(Xu1(:, 38), dUdB1(:, 38), ...
%      Xu2(:, 75), dUdB2(:, 75), ...
%      Xu3(:, 112), dUdB3(:, 112), ...
%      Xu4(:, 150), dUdB4(:, 150), ...
%      Xu5(:, 187), dUdB5(:, 187), ...
%      Xu6(:, 225), dUdB6(:, 225), ...
%      Xu7(:, 262), dUdB7(:, 262), ...
%      Xu8(:, 300), dUdB8(:, 300), ...
%      'linewidth', linewidth)
% legend('CS1', 'CS2', 'CS3', 'CS4', 'CS5', 'CS6', 'CS7', 'CS8')
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('U (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% % legend('100X75', '200X150', '300X225', '400X300', 'CSA', 'location', 'best')
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.

% figure,
% plot(Xu2(:, 75), Up2(:, 75), ...
%      Xu3(:, 112), Up3(:, 112), ...
%      Xu4(:, 150), Up4(:, 150), ...
%      Xu5(:, 187), Up5(:, 187), ...
%      Xu6(:, 225), Up6(:, 225), ...
%      Xu7(:, 262), Up7(:, 262), ...
%      Xu8(:, 300), Up8(:, 300), ...
%      'linewidth', linewidth)
% legend('CSA2', 'CSA3', 'CSA4', 'CSA5', 'CSA6', 'CSA7', 'CSA8')
% xlabel('X (m)', 'fontsize', fontsize)
% ylabel('U (m/s)', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% % legend('100X75', '200X150', '300X225', '400X300', 'CSA', 'location', 'best')
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
break
%% Convergence order
% x = 2;
% [a, b] = find(Xu1(:, 38) < x); U1ind = U1(a(end) + 1, 38);
% [a, b] = find(Xu2(:, 75) < x); U2ind = U2(a(end) + 1, 75);
% [a, b] = find(Xu3(:, 112) < x); U3ind = U3(a(end) + 1, 112);
% [a, b] = find(Xu4(:, 150) < x); U4ind = U4(a(end) + 1, 150);
% 
% err = abs([U1ind, U2ind, U3ind, U4ind] - U4ind) ./ U4ind;
% meshSize = [100, 200, 300, 400];
% f = fit(meshSize',err','power1')
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% loglog(meshSize, err, 'k', ...
%     meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
%     'linewidth', linewidth)
% xlabel('Number of nodes in X direction', 'fontsize', fontsize)
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
% x = 0.5;
% y = 0.75;
% 
% [a, b] = find(Xv1(:, 1) < x); x1ind = a(end) + 1;
% [a, b] = find(abs(Yv1(1, :)) < y); y1ind = b(end) + 1;
% [Xv1(x1ind, y1ind) Yv1(x1ind, y1ind)];
% V1ind = V1(x1ind, y1ind);
% 
% [a, b] = find(Xv2(:, 1) < x); x2ind = a(end) + 1;
% [a, b] = find(abs(Yv2(1, :)) < y); y2ind = b(end) + 1;
% [Xv2(x2ind, y2ind) Yv2(x2ind, y2ind)];
% V2ind = V2(x2ind, y2ind);
% 
% [a, b] = find(Xv3(:, 1) < x); x3ind = a(end) + 1;
% [a, b] = find(abs(Yv3(1, :)) < y); y3ind = b(end) + 1;
% [Xv3(x3ind, y3ind) Yv3(x3ind, y3ind)];
% V3ind = V3(x3ind, y3ind);
% 
% [a, b] = find(Xv4(:, 1) < x); x4ind = a(end) + 1;
% [a, b] = find(abs(Yv4(1, :)) < y); y4ind = b(end) + 1;
% [Xv4(x4ind, y4ind) Yv4(x4ind, y4ind)];
% V4ind = V4(x4ind, y4ind);
% 
% err = abs([V1ind, V2ind, V3ind, V4ind] - V4ind) ./ V4ind;
% meshSize = [100, 200, 300, 400];
% f = fit(meshSize',err','power1')
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% loglog(meshSize, err, 'k', ...
%     meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
%     'linewidth', linewidth)
% xlabel('Number of nodes in X direction', 'fontsize', fontsize)
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
% x = 2;
% [a, b] = find(Xp1(:, 38) < x); P1ind = P1(a(end) + 1, 38);
% [a, b] = find(Xp2(:, 75) < x); P2ind = P2(a(end) + 1, 75);
% [a, b] = find(Xp3(:, 112) < x); P3ind = P3(a(end) + 1, 112);
% [a, b] = find(Xp4(:, 150) < x); P4ind = P4(a(end) + 1, 150);
% 
% err = abs([P1ind, P2ind, P3ind, P4ind] - P4ind) ./ abs(P4ind);
% meshSize = [100, 200, 300, 400];
% f = fit(meshSize',err','power1')
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% loglog(meshSize, err, 'k', ...
%     meshSize(1:end-1), f.a .* meshSize(1:end-1) .^ f.b, 'k--', ...
%     'linewidth', linewidth)
% xlabel('Number of nodes in X direction', 'fontsize', fontsize)
% ylabel('Absolute error', 'fontsize', fontsize)
% set(gca, 'fontsize', fontsize)
% set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'renderer','painters')
% set(gca,'YMinorTick','on')
% set(gca,'XMinorTick','on')
% grid('on')
% grid minor
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.