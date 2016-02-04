clc;
clear all;
close all;
format short g;
% ----------------------------------------------------------------------- %
xStart = 0.0;
xEnd = 1.0;
yStart = 0.0;
yEnd = 1.0;
nx = 200;      % number of x-gridpoints
ny = 200;      % number of y-gridpoints

xc = 0.5;
yc = 0.5;
R = 0.1;
dR = 0.000001;
N = 300;

xImm = xc + R * cos(linspace(0,2*pi,N)');
yImm = yc + R * sin(linspace(0,2*pi,N)');

dxImm_db = cos(linspace(0,2*pi,N)');
dyImm_db = sin(linspace(0,2*pi,N)');

[deltaMatU1,deltaMatUSensitivity] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'U');

Rnew = R + dR;

xImm = xc + Rnew * cos(linspace(0,2*pi,N)');
yImm = yc + Rnew * sin(linspace(0,2*pi,N)');

dxImm_db = cos(linspace(0,2*pi,N)');
dyImm_db = sin(linspace(0,2*pi,N)');

[deltaMatU2,deltaMatUSensitivity2] = findBoundary(xStart,xEnd,yStart,yEnd,nx,ny,xImm,yImm,dxImm_db,dyImm_db,'U');

deltaMatUSensitivity_FD = (deltaMatU2 - deltaMatU1) / (Rnew - R);
max(max(abs(full(deltaMatUSensitivity_FD - deltaMatUSensitivity))))