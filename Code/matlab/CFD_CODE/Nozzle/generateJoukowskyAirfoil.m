clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
% function flag = openFOAM_joukowsky_airfoil_generator(x0,y0)
fontsize = 40;
linewidth = 5;
markersize = 20;
% ----------------------------------------------------------------------- %
% Joukowski showed that the image of a circle passing through (1.0,0) and
% containing the point (-1,0) is mapped onto a curve shaped like the
% cross section of an airplane wing.

N = 150; % Number of points on airfoil surface
% x0 = -0.05; % Defines maximum thickness. [-1 0]
x0 = -0.05; % Defines maximum thickness. [-1 0]
dx0 = 0.00; % Perturbation in thickness
xNew = x0 + dx0;
x0 = x0 + dx0;

% y0 = 0.3 + 0.0001; % Defines camber. Positive values for up camber [0 1]
y0 = 0.3 + 0.0000; % Defines camber. Positive values for up camber [0 1]
dy0 = 0.000; % Perturbation in camber
yNew = y0 + dy0;
% y0 = y0 - dy0;

b = sqrt(0.9);

R = sqrt((1 - x0)^2 + y0^2); % Calculate R based on (1.0,0) on circle
RNew = sqrt((1 - xNew)^2 + yNew^2); % Calculate R based on (1.0,0) on circle
% Checking the (-1,0) is inside
if (1 + x0)^2 + y0^2 - R^2 >= 0
    disp(['Cannot generate airfoil shape'])
%     return;
end

if (1 + xNew)^2 + yNew^2 - RNew^2 >= 0
    disp(['Cannot generate airfoil shape'])
%     return;
end

theta = linspace(0,2*pi,N);

x = x0 + R .* cos(theta);
y = y0 + R .* sin(theta);
z = x + i * y;

xNew = xNew + RNew .* cos(theta);
yNew = yNew + RNew .* sin(theta);
zNew = xNew + i * yNew;

z = z + b ./ z; % Joukowski transform
X = real(z);
Y = imag(z);

zNew = zNew + b ./ zNew; % Joukowski transform
XNew = real(zNew);
YNew = imag(zNew);

Y = Y / (max(X) - min(X));
X = X / (max(X) - min(X));
X = X + 0.5 - max(X);
[a,b] = max(X);
Y = Y - Y(b);

YNew = YNew / (max(XNew) - min(XNew));
XNew = XNew / (max(XNew) - min(XNew));
XNew = XNew + 0.5 - max(XNew);
[a,b] = max(XNew);
YNew = YNew - Y(b);

pointCloud = [X',Y'];
if dy0 == 0
    pointCloudSensitivity = ([XNew',YNew'] - [X',Y']) / dx0;
elseif dx0 == 0
	pointCloudSensitivity = ([XNew',YNew'] - [X',Y']) / dy0;
end

distX = abs(pointCloud(2:end,1) - pointCloud(1:end-1,1));
meanDistX = mean(distX);
j = 1;
for i=j:length(pointCloud)-1
    for j=i+1:length(pointCloud)-1
        if j > length(pointCloud)
            break
        end
        if abs(pointCloud(i,1) - pointCloud(j,1)) < meanDistX
            pointCloud(j,:) = [];
        else
            break
        end
    end
end

% figure,
% plot(X,Y,'ko',...
%      pointCloud(:,1),pointCloud(:,2),'r+')
% axis([-1 1 -1 1])
% length(pointCloud)

figure,
plot(X,Y,'k',...
    'linewidth',linewidth)
xlabel('X','fontsize',fontsize)
ylabel('Y','fontsize',fontsize)
axis([-0.6 0.6 -0.6 0.6])
set(gca,'fontsize',fontsize,'linewidth',linewidth)
axis('equal')


dlmwrite('point_cloud.txt',pointCloud);
dlmwrite('point_cloud_sensitivity.txt',pointCloudSensitivity);
% copyfile('point_cloud.txt','current_simulation_results/')

