clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
% data = dlmread('NACA2408.txt');
% data = dlmread('NACA4418.txt');
% data = dlmread('NACA1410.txt');
data = dlmread('NACA001264.txt');

x = data(:,1);
y = data(:,2);

[ind,vla] = find(x == 0);
xTop = x(1:ind); yTop = y(1:ind);
xTop = flipud(xTop); yTop = flipud(yTop);
xBottom = x(ind + 1:end); yBottom = y(ind + 1:end);

area = 0;
xArea = 0;
yArea = 0;
for i = 1:length(xTop) - 1
    area = area + (xTop(i+1) - xTop(i)) * (yTop(i) + yTop(i+1)) * 0.5;
    xArea = xArea + (xTop(i+1) + xTop(i)) * 0.5 * ...
        (xTop(i+1) - xTop(i)) * (yTop(i) + yTop(i+1)) * 0.5;
    yArea = yArea + (yTop(i+1) + yTop(i)) * 0.5 * ...
        (xTop(i+1) - xTop(i)) * (yTop(i) + yTop(i+1)) * 0.5;
end
for i = 1:length(xBottom) - 1
    area = area + (xBottom(i+1) - xBottom(i)) * abs(yBottom(i) + yBottom(i+1)) * 0.5;
    xArea = xArea + (xBottom(i+1) + xBottom(i)) * 0.5 * ...
        (xBottom(i+1) - xBottom(i)) * abs(yBottom(i) + yBottom(i+1)) * 0.5;
    yArea = yArea + (yBottom(i+1) + yBottom(i)) * 0.5 * ...
        (xBottom(i+1) - xBottom(i)) * abs(yBottom(i) + yBottom(i+1)) * 0.5;
end
xBar = xArea / area;
yBar = yArea / area;

theta = -10;
rotMat = [cosd(theta),-sind(theta);sind(theta),cosd(theta)];
xTop = xTop - xBar;
xBottom = xBottom - xBar;
yTop = yTop - yBar;
yBottom = yBottom - yBar;

for i=1:length(xTop)
    temp = rotMat * [xTop(i);yTop(i)];
    xTop(i) = temp(1);
    yTop(i) = temp(2);
end

for i=1:length(xBottom)
    temp = rotMat * [xBottom(i);yBottom(i)];
    xBottom(i) = temp(1);
    yBottom(i) = temp(2);
end

xTop = xTop + xBar;
xBottom = xBottom + xBar;
yTop = yTop + yBar;
yBottom = yBottom + yBar;

x = [xBottom(1);xTop(1:2)];
y = [yBottom(1);yTop(1:2)];

coef = [y(1)^2,y(1),1;y(2)^2,y(2),1;y(3)^2,y(3),1] \ [x(1);x(2);x(3)];

X = @(Y) coef(1) * Y.^2 + coef(2) * Y + coef(3);

% figure(1),
% plot(xTop,yTop,'ro',...
%     xBottom,yBottom,'ko',...
%     linspace(x(1),x(end),10),Y(linspace(x(1),x(end),10)),'g+')
% axis equal
% pause(3)

c = [];
airfoilInter = [X(linspace(y(1),y(end),10)'),linspace(y(1),y(end),10)'];

for iTop=2:length(xTop)-1
    x = [xTop(iTop-1:iTop+1)];
    y = [yTop(iTop-1:iTop+1)];
    coef = [x(1)^2,x(1),1;x(2)^2,x(2),1;x(3)^2,x(3),1] \ [y(1);y(2);y(3)];
    Y = @(X) coef(1) * X.^2 + coef(2) * X + coef(3);
    
    airfoilInter = [airfoilInter;linspace(x(2),x(3),10)',Y(linspace(x(2),x(3),10)')];
    
%     figure(1),
%     plot(xTop,yTop,'ro',...
%         xBottom,yBottom,'ko',...
%      linspace(x(1),x(end),10),Y(linspace(x(1),x(end),10)),'g+')
%     axis equal
%     pause(3)
end

for iBottom=2:length(xBottom)-1
    x = [xBottom(iBottom-1:iBottom+1)];
    y = [yBottom(iBottom-1:iBottom+1)];
    coef = [x(1)^2,x(1),1;x(2)^2,x(2),1;x(3)^2,x(3),1] \ [y(1);y(2);y(3)];
    Y = @(X) coef(1) * X.^2 + coef(2) * X + coef(3);
    
    airfoilInter = [airfoilInter;linspace(x(2),x(3),10)',Y(linspace(x(2),x(3),10)')];
%     figure(1),
%     plot(xTop,yTop,'ro',...
%         xBottom,yBottom,'ko',...
%         linspace(x(1),x(end),10),Y(linspace(x(1),x(end),10)),'g+')
%     axis equal
%     pause(3)
end

figure(1),
plot(data(:,1),data(:,2),'r*',...
    xTop,yTop,'ro',...
    xBottom,yBottom,'ko',...
    airfoilInter(:,1),airfoilInter(:,2),'g+')
axis equal

pointCloud = airfoilInter;
pointCloudSensitivity = airfoilInter;

dlmwrite('point_cloud.txt',pointCloud);
dlmwrite('point_cloud_sensitivity.txt',pointCloudSensitivity);