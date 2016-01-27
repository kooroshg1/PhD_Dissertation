clc;
close all;
fontsize = 30;
linewidth = 5;

figure,
plot(linspace(0,tf,length(fxHist)),fxHist([1,35,50],:))
legend('\theta = 0','\theta = 135','\theta = 180')
xlabel('Time','fontsize',fontsize)
ylabel('f_x','fontsize',fontsize)
title('force value in x direction at different location on the cylinder','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
plot(linspace(0,tf,length(fyHist)),fyHist([1,35,50],:))
legend('\theta = 0','\theta = 135','\theta = 180')
xlabel('Time','fontsize',fontsize)
ylabel('f_y','fontsize',fontsize)
title('force value in y direction at different location on the cylinder','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
plot(linspace(0,360,length(uBoundary)),uBoundary,'k','linewidth',linewidth)
xlabel('\theta on cylinder','fontsize',fontsize)
ylabel('Flow velocity','fontsize',fontsize)
title('Absolute velocity of fluid on cylinder surface','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
plot(linspace(0,tf,nt),centerLocation,'k','linewidth',linewidth)
xlabel('Time','fontsize',fontsize)
ylabel('Cylinder center location','fontsize',fontsize)
title('Time history of location of center of cylinder','fontsize',fontsize)
set(gca,'fontsize',fontsize)

figure,
plot(linspace(0,tf,nt-1),diff(centerLocation)/dt,'k','linewidth',linewidth)
xlabel('Time','fontsize',fontsize)
ylabel('Cylinder center velocity','fontsize',fontsize)
title('Time history of velocity of center of cylinder','fontsize',fontsize)
set(gca,'fontsize',fontsize)

% figure,
% plot(linspace(0,tf,nt-1-1),diff(diff(centerLocation)/dt)/dt,'k','linewidth',linewidth)
% xlabel('Time','fontsize',fontsize)
% ylabel('Cylinder center velocity','fontsize',fontsize)
% title('Time history of velocity of center of cylinde','fontsize',fontsize)
% set(gca,'fontsize',fontsize)