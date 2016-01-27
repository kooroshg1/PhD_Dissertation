% clc;
% clear all;
% close all;
% format short g;
% ----------------------------------------------------------------------- %
function [deltaMat,deltaMatSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xStar,yStar,dxStar_db,dyStar_db,meshType)
% N = 7; % Number of lagrange points
% xStart = 0.0;
% xEnd = 1.0;
% yStart = 0.0;
% yEnd = 1.0;
% nx = 20;      % number of x-gridpoints
% ny = 20;      % number of y-gridpoints
% xStar = 0.5 + 0.1 * cos(linspace(0,2*pi,N));
% yStar = 0.5 + 0.1 * sin(linspace(0,2*pi,N));
% meshType = 'U';
% dxStar_db = cos(linspace(0,2*pi,N));
% dyStar_db = sin(linspace(0,2*pi,N));

dx = (xEnd - xStart) / nx;
dy = (yEnd - yStart) / ny;
if strcmp(meshType,'U')
    x = xStart+dx:dx:xEnd-dx;
    y = yStart+dy/2:dy:yEnd-dy/2;
    nx = nx - 1;
elseif strcmp(meshType,'V')
    x = xStart+dx/2:dx:xEnd-dx/2;
    y = yStart+dy:dy:yEnd-dy;
    ny = ny - 1;
elseif strcmp(meshType,'P')
    x = (xStart + dx / 2):dx:(xEnd - dx / 2);
    y = (yStart + dy / 2):dy:(yEnd - dy / 2);
end
index = 1:nx*ny;
[Y,X] = meshgrid(y,x);

xVec = reshape(X,[],1);
yVec = reshape(Y,[],1);
% nearPointsIndex = sparse(length(xStar),8);
% nearPointsIndexSensitivity = sparse(length(xStar),8);
% nearPointsCoordinate = sparse(4*length(xStar),2);
% 
% for i=1:length(xStar)
%     distVec = sqrt((xVec - xStar(i)).^2 + (yVec - yStar(i)).^2);
%     ddistVec_db = ( 2 * (xVec - xStar(i)) * dxStar_db(i) + ...
%                     2 * (yVec - yStar(i)) * dyStar_db(i) ) ./ ...
%                   ( 2 * distVec );
%     closePoints = sortrows([distVec,xVec,yVec,index',ddistVec_db],1);
%     
%     nearPointsCoordinate((i-1)*4+1:i*4,:) = closePoints(1:4,[2,3]);
%     
%     nearPointsIndex(i,:) = [closePoints(1:4,4)' closePoints(1:4,1)'];
%     nearPointsIndexSensitivity(i,:) = [closePoints(1:4,4)' closePoints(1:4,5)'];
% end
% 
% nearPointsIndexSensitivity = [nearPointsIndexSensitivity(:,1:4) ...
%     (nearPointsIndexSensitivity(:,5:8) .* (sum(nearPointsIndex(:,5:8)')' * ones(1,4)) + ...
%     (sum(nearPointsIndexSensitivity(:,5:8)')' * ones(1,4)) .* nearPointsIndex(:,5:8)) ./ ...
%     (sum(nearPointsIndex(:,5:8)')' * ones(1,4)).^2];
% 
% nearPointsIndex = [nearPointsIndex(:,1:4) (nearPointsIndex(:,5:8)) ./ ...
%     (sum(nearPointsIndex(:,[5:8])')' * ones(1,4))];
% 
% deltaMat = zeros(size(nearPointsIndex,1),nx*ny);
% deltaMatSensitivity = zeros(size(nearPointsIndexSensitivity,1),nx*ny);
% 
% for i=1:size(nearPointsIndex,1)
%     deltaMat(i,nearPointsIndex(i,1:4)) = nearPointsIndex(i,5:8);
%     deltaMatSensitivity(i,nearPointsIndexSensitivity(i,1:4)) = nearPointsIndexSensitivity(i,5:8);
% end

nearIndex = sparse(length(xStar),4);
M = sparse(length(xStar),nx*ny);
dMdB = sparse(length(xStar),nx*ny);

leftBoundaryIndex = 1:nx:nx*ny;
rightBoundaryIndex = nx:nx:nx*ny;
topBoundaryIndex = leftBoundaryIndex(end):rightBoundaryIndex(end);
bottomBoundaryIndex = leftBoundaryIndex(1):rightBoundaryIndex(1);
for iImm = 1:length(xStar)
    L2Edist = sortrows([sqrt((xVec - xStar(iImm)).^2 + (yVec - yStar(iImm)).^2), (1:nx*ny)'],1);
    temp = L2Edist(1:2,2);
    temp = sort(temp);
    if (temp(end) - temp(1)) == nx && (xVec(temp(end)) - xStar(iImm)) <= 0
        temp = [temp;temp + 1];
    elseif (temp(end) - temp(1)) == nx && (xVec(temp(end)) - xStar(iImm)) > 0
        temp = [temp;temp - 1];
    elseif (temp(end) - temp(1)) == 1 && (yVec(temp(end)) - yStar(iImm)) > 0
        temp = [temp;temp - nx];
    elseif (temp(end) - temp(1)) == 1 && (yVec(temp(end)) - yStar(iImm)) <= 0
        temp = [temp;temp + nx];
    end
    temp = sort(temp);
    if mod(temp(1),nx) == 0
        temp(2) = temp(1) - 1;
        temp(4) = temp(3) - 1;
        temp = sort(temp);
    end
    if max(temp(1) == topBoundaryIndex) == 1
        temp(3) = temp(1) - nx;
        temp(4) = temp(2) - nx;
        temp = sort(temp);
    end
    if max(temp(3) == bottomBoundaryIndex) == 1
        temp(1) = temp(3) + nx;
        temp(2) = temp(4) + nx;
        temp = sort(temp);
    end
%     [xStar(iImm),yStar(iImm)]
%     figure(1),
%     plot(xStar(iImm),yStar(iImm),'r*',...
%          xVec,yVec,'bo',...
%          xVec(temp),yVec(temp),'r+')
% 	pause(0.1)
%     temp
    nearIndex(iImm,:) = temp';
    
    % Generate interpolation matrix
    u = (xStar(iImm) - xVec(temp(1))) / (xVec(temp(2)) - xVec(temp(1)));
    v = (yStar(iImm) - yVec(temp(1))) / (yVec(temp(3)) - yVec(temp(1)));
%     if (v < 0)
%         disp('HA!!')
%     end
    
    dudb = dxStar_db(iImm) / (xVec(temp(2)) - xVec(temp(1)));
    dvdb = dyStar_db(iImm) / (yVec(temp(3)) - yVec(temp(1)));
    
%     Minter = [(1-u)*(1-v) u*(1-v) (1-u)*v u*v];
    M(iImm,temp(1)) = (1-u)*(1-v);
    M(iImm,temp(2)) = u*(1-v);
    M(iImm,temp(3)) = (1-u)*v;
    M(iImm,temp(4)) = u*v;
    
    dMdB(iImm,temp(1)) = (-dudb)*(1-v) + (1-u)*(-dvdb);
    dMdB(iImm,temp(2)) = dudb*(1-v) + u*(-dvdb);
    dMdB(iImm,temp(3)) = (-dudb)*v + (1-u)*dvdb;
    dMdB(iImm,temp(4)) = dudb*v + u*dvdb;
end

deltaMat = M;
deltaMatSensitivity = dMdB;

% figure,
% plot(xVec,yVec,'ko',...
%     xStar,yStar,'r+',...
%     nearPointsCoordinate(:,1),nearPointsCoordinate(:,2),'g+')
% axis equal
