clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
xFront = 0.475;
xBack = 0.525;
yBottom = 0.0;
yTop = 0.3;

nx = 10 * 3;
ny = 30 * 3;

pointCloud = [xFront*ones(ny,1) linspace(yBottom,yTop,ny)';...
              linspace(xFront,xBack,nx)' yTop*ones(nx,1);...
              xBack*ones(ny,1) linspace(yTop,yBottom,ny)']
findSame = pointCloud(1:end-1,:) == pointCloud(2:end,:);
findSame = findSame(:,1) == findSame(:,2);
pointCloud(findSame,:) = [];

xFront = 0.475;
xBack = 0.525;
yBottom = 0.0;
yTop = 0.3 + 0.01;

% nx = 10;
% ny = 30;

pointCloudPerturb = [xFront*ones(ny,1) linspace(yBottom,yTop,ny)';...
              linspace(xFront,xBack,nx)' yTop*ones(nx,1);...
              xBack*ones(ny,1) linspace(yTop,yBottom,ny)'];
findSame = pointCloudPerturb(1:end-1,:) == pointCloudPerturb(2:end,:);
findSame = findSame(:,1) == findSame(:,2);
pointCloudPerturb(findSame,:) = [];

pointCloudSensitivity = (pointCloudPerturb - pointCloud) / 0.01;

dlmwrite('point_cloud.txt',pointCloud);
dlmwrite('point_cloud_sensitivity.txt',pointCloudSensitivity);