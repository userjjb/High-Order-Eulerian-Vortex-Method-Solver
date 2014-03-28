close all
clear all
clc
% Driver script for solving the 1D advection equations

%Global vars and RK4 coeffs
Globals1D;

% Order of polymomials used for approximation 
N = 3;
alpha = 1; %control central flux - upwind flux

% Generate simple 1D mesh
[Nv, VX, K, EToV] = MeshGen1D(0.0,2.0,10);

% Initialize solver and construct grid and metric
StartUp1D;

% Set initial conditions
u = sin(x);



% Solve Problem
FinalTime = 10; %Length of runtime

% compute time step size
xmin = min(abs(x(1,:)-x(2,:)));
CFL=0.75; dt   = CFL/(2*pi)*xmin; dt = .5*dt;
Nsteps = ceil(FinalTime/dt); dt = FinalTime/Nsteps; 

LogTime=0.02;   %Time period to log
TimeTol=dt*1.05;   %How close mod(time,LogTime) should be to 0
[u,saved] = Advec1D(u,FinalTime,LogTime,TimeTol,alpha);
beep;

for i=2:nT
    plot(x,saved(:,:,i),'o-')
    axis([0 2 -1.2 1.2])
%     hold on
%     plot([0,2],[1.2,-1.2],'w.')
%     hold off
    text(2,1,num2str((i/nT)*FinalTime));
    pause(0.001)
end