clc;
clear;
close all;

global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();
% In this file we will run linear advection using SSP3 for temporal
% discretization and WENO5 for the following conditions:
% Initial conditions: smooth (sinusoidal), hat, step (discontinuous)
% Boudnary conditions: periodic: 

for i=1:1
    if i ==1 
        init_conds = 'sinIC';
   
    else
    end
    
    if init_conds == "sinIC"
        N_pts_burger= 100;
        fn_ex = @(x,t) burger_sol(x,t,N_pts_burger);
    else 
    end

    var_out = fn_SSP3WENO_burger(init_conds,fn_ex);
    var_out = fn_plot_WENO_burger(init_conds);   
    
end
function [] = gauss
% For n point gauss quadrature, return evaluation points and weights for
% gauss quadrature over [-1,1]
global nq xq wq
if nq == 1
    xq = 0;
    wq = 1;
elseif nq == 2
    xq = [-0.57735026918962576451, 0.57735026918962576451];
    wq = [0.5, 0.5];
elseif nq == 3
    xq = [-0.7745966692414834, 0, 0.7745966692414834];
    wq = 0.5*[0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
elseif nq > 3
    fprintf(1,'No Gauss quadrature implemented of this degree\n')
    return
end
end

