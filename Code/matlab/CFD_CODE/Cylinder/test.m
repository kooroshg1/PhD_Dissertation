clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
pointCloud = dlmread('point_cloud.txt');
nodeEffectiveArea = zeros(size(pointCloud,1),1);
nodeUnitNormalVector = zeros(size(pointCloud,1),2);
pBoundary = dlmread('surface_pressure.txt');
for i=1:size(pointCloud,1)
    if i == 1
        nodeEffectiveArea(i) = 0.5 * sqrt(...
        (pointCloud(i+1,1) - pointCloud(end,1))^2 + ...
        (pointCloud(i+1,2) - pointCloud(end,2))^2);

        alpha = pointCloud(i+1,1) - pointCloud(end,1);
        beta = pointCloud(i+1,2) - pointCloud(end,2);
        nodeUnitNormalVector(i,1) = beta / sqrt(alpha^2 + beta^2);
        nodeUnitNormalVector(i,2) = -alpha / sqrt(alpha^2 + beta^2);
    elseif i == size(pointCloud,1)
        nodeEffectiveArea(i) = 0.5 * sqrt(...
        (pointCloud(1,1) - pointCloud(end-1,1))^2 + ...
        (pointCloud(1,2) - pointCloud(end-1,2))^2);

        alpha = pointCloud(1,1) - pointCloud(end-1,1);
        beta = pointCloud(1,2) - pointCloud(end-1,2);
        nodeUnitNormalVector(i,1) = beta / sqrt(alpha^2 + beta^2);
        nodeUnitNormalVector(i,2) = -alpha / sqrt(alpha^2 + beta^2);
    else
        nodeEffectiveArea(i) = 0.5 * sqrt(...
            (pointCloud(i+1,1) - pointCloud(i-1,1))^2 + ...
            (pointCloud(i+1,2) - pointCloud(i-1,2))^2);

        alpha = pointCloud(i+1,1) - pointCloud(i-1,1);
        beta = pointCloud(i+1,2) - pointCloud(i-1,2);
        nodeUnitNormalVector(i,1) = beta / sqrt(alpha^2 + beta^2);
        nodeUnitNormalVector(i,2) = -alpha / sqrt(alpha^2 + beta^2);
    end
end

[pBoundary .* nodeEffectiveArea .* -nodeUnitNormalVector(:,1),...
    pBoundary .* nodeEffectiveArea,...
    nodeUnitNormalVector(:,1)]

sum(pBoundary .* nodeEffectiveArea .* -nodeUnitNormalVector(:,1))

figure,
quiver(pointCloud(:,1),pointCloud(:,2),nodeUnitNormalVector(:,1),nodeUnitNormalVector(:,2))