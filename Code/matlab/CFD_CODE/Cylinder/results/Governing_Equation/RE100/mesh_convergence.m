clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
fontsize = 52;
linewidth = 5.0;
markersize = 15.0;
skip = 10;
% ----------------------------------------------------------------------- %
uB1 = dlmread('100X75/uBoundary.txt');
vB1 = dlmread('100X75/vBoundary.txt');
pB1 = dlmread('100X75/pBoundary.txt');

uB2 = dlmread('200X150/uBoundary.txt');
vB2 = dlmread('200X150/vBoundary.txt');
pB2 = dlmread('200X150/pBoundary.txt');

uB3 = dlmread('300X225/uBoundary.txt');
vB3 = dlmread('300X225/vBoundary.txt');
pB3 = dlmread('300X225/pBoundary.txt');

uB4 = dlmread('400X300/uBoundary.txt');
vB4 = dlmread('400X300/vBoundary.txt');
pB4 = dlmread('400X300/pBoundary.txt');

theta = [linspace(-90, 0, length(pB1) / 2), ...
         linspace(0, 90, length(pB1) / 2)];
     
% figure,
% set(gcf, 'Position', get(0,'Screensize'));
% plot(theta, uB4, 'k', ...
%      'linewidth', linewidth)
% xlabel('\theta', 'fontsize', fontsize)
% ylabel('u-velocity (m/s)', 'fontsize', fontsize)
% ylim([-1, 1])
% set(gca, 'fontsize', fontsize)
% set(gcf,'renderer','painters')
% set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
% grid('on')
% grid minor
% set(gca,'LooseInset',get(gca,'TightInset'))

figure,
set(gcf, 'Position', get(0,'Screensize'));
plot(theta, vB4, 'k', ...
     'linewidth', linewidth)
xlabel('\theta', 'fontsize', fontsize)
ylabel('v-velocity (m/s)', 'fontsize', fontsize)
ylim([-1, 1])
set(gca, 'fontsize', fontsize)
set(gcf,'renderer','painters')
set(gcf, 'PaperPosition', [0.25 2.5 8 6]); % last 2 are width/height.
grid('on')
grid minor
set(gca,'LooseInset',get(gca,'TightInset'))