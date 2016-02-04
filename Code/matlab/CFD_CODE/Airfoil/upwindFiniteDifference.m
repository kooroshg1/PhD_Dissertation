% clc;
% clear all;
% close all;
% format short g;
% % ----------------------------------------------------------------------- %
% originalDirection = 'y';
% neededDirection = 'y';
% nx = 6;
% ny = 6;
% dt = 0.1;
% upwindMatrix = floor(10 * rand(nx+1,ny-1) - 5);
% flowDirection = floor(10 * rand(nx+1,ny-1) - 5);
% gamma = 1;
function upwindMatrix = upwindFiniteDifference(upwindMatrix,flowDirection,originalDirection,neededDirection,gamma)

if strcmp(neededDirection,'x') && strcmp(originalDirection,'x')
    flowDirection = upwindMatrix;
    flowDirectionAtMidPoint = (flowDirection(1:end-1,:) + ...
                               flowDirection(2:end,:)) / 2;
                           
    flowDirectionAtMidPoint = 2 * ceil(heaviside(flowDirectionAtMidPoint)) - 1;
    
    upwindMatrix = (upwindMatrix(1:end-1,:) + upwindMatrix(2:end,:)) / 2 - gamma .* flowDirectionAtMidPoint .* ...
    (upwindMatrix(2:end,:) - upwindMatrix(1:end-1,:)) / 2;
elseif strcmp(neededDirection,'y') && strcmp(originalDirection,'x')
    flowDirectionAtMidPoint = (flowDirection(1:end-1,:) + ...
                               flowDirection(2:end,:)) / 2;
	
	flowDirectionAtMidPoint = 2 * ceil(heaviside(flowDirectionAtMidPoint)) - 1;
    
    upwindMatrix = (upwindMatrix(:,2:end) + upwindMatrix(:,1:end-1)) / 2 - ...
        gamma .* flowDirectionAtMidPoint .* (upwindMatrix(:,2:end) - upwindMatrix(:,1:end-1)) / 2;
elseif strcmp(neededDirection,'x') && strcmp(originalDirection,'y')
    flowDirectionAtMidPoint = (flowDirection(:,2:end) + ...
                               flowDirection(:,1:end-1)) / 2;
                           
	flowDirectionAtMidPoint = 2 * ceil(heaviside(flowDirectionAtMidPoint)) - 1;
    
    upwindMatrix = (upwindMatrix(2:end,:) + upwindMatrix(1:end-1,:)) / 2 - ...
        gamma .* flowDirectionAtMidPoint .* (upwindMatrix(2:end,:) - upwindMatrix(1:end-1,:)) / 2;
elseif strcmp(neededDirection,'y') && strcmp(originalDirection,'y')
    flowDirection = upwindMatrix;
    flowDirectionAtMidPoint = (flowDirection(:,2:end) + ...
                               flowDirection(:,1:end-1)) / 2;
                           
	flowDirectionAtMidPoint = 2 * ceil(heaviside(flowDirectionAtMidPoint)) - 1;
    
    upwindMatrix = (upwindMatrix(:,2:end) + upwindMatrix(:,1:end-1)) / 2 - ...
        gamma  .* flowDirectionAtMidPoint .*(upwindMatrix(:,2:end) - upwindMatrix(:,1:end-1)) / 2;
end
