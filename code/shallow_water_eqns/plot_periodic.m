clc;
clear;
close all;
T=40;
u_init=1;
v_init =1;
which_rec = "normal_one"; 
pow=1;

xlim_p = T;
ylim_p1 = 1e-10;
ylim_p2 = 5e-7;
EI_plot_y_lim = 10;
R_plot_lims = [1e-11, 2e-7];
scheme_name = 'SHW_RK3_WENO3_rec3_fixed';
% var_out = fn_plot_shw_periodic(1,1,ylim_p1, ylim_p2, xlim_p , EI_plot_y_lim, R_plot_lims);
var_out = fn_plot_shw_periodic(1,1,ylim_p1, ylim_p2, xlim_p , EI_plot_y_lim, R_plot_lims,pow,scheme_name);
