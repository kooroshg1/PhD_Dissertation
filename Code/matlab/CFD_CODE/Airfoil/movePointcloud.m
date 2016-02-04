function movePointcloud(pBoundary,dt,nodeEffectiveArea,nodeUnitNormalVector)

pointCloud = dlmread('point_cloud.txt');
pointCloudVelocity = dlmread('point_cloud_velocity.txt');
% pBoundary = dlmread('surface_pressure.txt');

m = 0.01; % Cylinder mass

% forceX = trapz(linspace(0,2*pi*0.1,length(pBoundary)),pBoundary);
% forceY = zeros(size(forceX,1),size(forceX,2));

forceX = sum(pBoundary .* nodeEffectiveArea .* -nodeUnitNormalVector(:,1));
forceY = sum(pBoundary .* nodeEffectiveArea .* -nodeUnitNormalVector(:,2));

ax = forceX / m;
ay = forceY / m;

pointCloudVelocity(:,1) = pointCloudVelocity(:,1) + ax * dt;
pointCloudVelocity(:,2) = pointCloudVelocity(:,2) + ay * dt;

pointCloud(:,1) = pointCloud(:,1) + pointCloudVelocity(:,1) * dt;
pointCloud(:,2) = pointCloud(:,2) + pointCloudVelocity(:,2) * dt;

dlmwrite('point_cloud.txt',pointCloud);
dlmwrite('point_cloud_velocity.txt',pointCloudVelocity);