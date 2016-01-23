clc;
% clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
xc = 0.5;
yc = 0.5;
dR = 0.0000;
% dR = 0.0001i;
R = 0.10 - dR;
N = 100;

x = xc + R * cos(linspace(0,2*pi,N)'); %x(end) = [];
y = yc + R * sin(linspace(0,2*pi,N)'); %y(end) = [];

pointCloud = [x y];

dxdR = cos(linspace(0,2*pi,N)'); %dxdR(end) = [];
dydR = sin(linspace(0,2*pi,N)'); %dydR(end) = [];

pointCloudSensitivity = [dxdR dydR];

dlmwrite('point_cloud.txt',pointCloud);
dlmwrite('point_cloud_sensitivity.txt',pointCloudSensitivity);