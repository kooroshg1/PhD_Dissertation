Re = 1e2;     % Reynolds number
dt = 1e-2;    % time step
tf = 5;    % final time
xStart = 0.0;
xEnd = 4.0;
yStart = -1.5;
yEnd = 1.5;
nx = 400;      % number of x-gridpoints
ny = 300;      % number of y-gridpoints
nsteps = 10;  % number of steps with graphic output

alpha = -1e4;
beta = -5e1;

xc = 0.5;
yc = 0.0;
dR = 0.0000;
% dR = 0.0001i;
R = 0.10 - dR;
N = 100;
