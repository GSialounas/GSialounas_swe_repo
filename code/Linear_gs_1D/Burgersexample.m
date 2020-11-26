% This matlab codes implement the 5th order RBF-WENO scheme
% method for 1D Burgers' equation
%
% Implementation : 1D WENO Finite Difference Lax-Friedrich splitting
% Time-integration: 3rd order Runge-Kutta method
%
% For reference, see the following paper:
% Radial basis function ENO and WENO finite difference
% methods based on the optimization of shape parameters, in revision
%
% Jingyang Guo and Jae-Hun Jung    May, 2016
%
% Test problem:
% Burger's Equation
% u_t - (1/2 u^2)_x = 0, -1 < x < 1
% u(x,0)=-sin(pi*x)
% u(0,t)=u(1,t)=0
clc
clear
close all

clear all;

global N; N = 400; 
C = 0.4;
T = 4; 

x = linspace(-1,1,N+1);
global dx; dx = 2/N;
dt = C*dx;
time_step = T/dt;

ub = -sin(pi*x);

% ghost cells and impose B.C.
vb_plus(4:N+4) = 1/2*(1/2*ub.*ub + ub);
vb_plus(1:3) = vb_plus(4);
vb_plus(N+5:N+7) = vb_plus(N+4);

vb_minus(4:N+4) = 1/2*(1/2*ub.*ub - ub);
vb_minus(1:3) = vb_minus(4);
vb_minus(N+5:N+7) = vb_minus(N+4);

time = 0;
for j=1:round(time_step)   
       
time = time + dt;
y0 = ub;

% RK1 =====================================================================
[vminus,~] = weno_RBF_k3(vb_plus); [~,vplus] = weno_RBF_k3(vb_minus);

fminus = vminus(2:N+2)+vplus(3:N+3);
fplus  = vminus(1:N+1)+vplus(2:N+2);

flux = -(fminus-fplus)/dx;
y1 = y0 + dt*flux;
ub = y1;
vb_plus(4:N+4) = 1/2*(1/2*ub.*ub + ub);
vb_minus(4:N+4) = 1/2*(1/2*ub.*ub - ub);

% RK2 =====================================================================
[vminus,~] = weno_RBF_k3(vb_plus); [~,vplus] = weno_RBF_k3(vb_minus);

fminus = vminus(2:N+2)+vplus(3:N+3);
fplus  = vminus(1:N+1)+vplus(2:N+2);

flux = -(fminus-fplus)/dx;
y1 = 3/4*y0+1/4*(ub + dt*flux);
ub = y1;
vb_plus(4:N+4) = 1/2*(1/2*ub.*ub + ub);
vb_minus(4:N+4) = 1/2*(1/2*ub.*ub - ub);

% RK3 =====================================================================
[vminus,~] = weno_RBF_k3(vb_plus); [~,vplus] = weno_RBF_k3(vb_minus);

fminus = vminus(2:N+2)+vplus(3:N+3);
fplus  = vminus(1:N+1)+vplus(2:N+2);

flux = -(fminus-fplus)/dx;
y1 = 1/3*y0+2/3*(ub + dt*flux);
ub = y1;
vb_plus(4:N+4) = 1/2*(1/2*ub.*ub + ub);
vb_minus(4:N+4) = 1/2*(1/2*ub.*ub - ub);
% =========================================================================
% for animation
clf();
plot(x,ub,'o-');
axis([-1,1,-1.2,1.2]);
pause(0.001)
end
