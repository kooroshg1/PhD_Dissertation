clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
U = dlmread('U.txt');
Upert = dlmread('Upert.txt');
dUdB = (Upert - U) / 0.01;

dUdB_CSA = dlmread('dUdB_CSA.txt');

xStart = -1.0;
xEnd = 2.0;
yStart = 0.0;
yEnd = 1.0;

pointCloud = dlmread('point_cloud.txt');
xImm = pointCloud(:,1);
yImm = pointCloud(:,2);

xu = linspace(xStart,xEnd,size(U,1));
yu = linspace(yStart,yEnd,size(U,2));
[Yu,Xu] = meshgrid(yu,xu);
figure,
contourf(Xu,Yu,dUdB,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal

figure,
contourf(Xu,Yu,dUdB_CSA,50,'linestyle','none')
hold all
plot(xImm,yImm,'w.','linewidth',2)
set(gca,'XTickLabel','','YTickLabel','')
axis equal

figure,
plot(xu,dUdB_CSA(:,50),'k',...
     xu,dUdB(:,50),'r')
legend('CSA','FD')