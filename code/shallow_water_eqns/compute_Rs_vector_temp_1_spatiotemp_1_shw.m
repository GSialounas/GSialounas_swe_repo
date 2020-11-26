function [L2Rt,L2Rt_h,L2Rt_hv,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv,IUt_h, IUt_hv] = compute_Rs_vector_temp_1_spatiotemp_1_shw(x,dist_x_pl,dist_x_min,h_old,h_new,evalt,tj,dt,f_h_old,f_h_new,spatial_disc,flux_fn)
global nq xq wq

g=9.81;
dx= dist_x_pl(1);
L2Rt = 0;
L2Rt_h = 0;

L2Rt_hv = 0;

L2Rt_arr = zeros(size(x));
c_0_coeff_arr_new=0; c_0_coeff_arr_old=0;


c_0_t = h_old;
c_1_t = 1/dt*(h_new-h_old);



diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
ht = c_0_t + c_1_t * diff_t ;
ht_t = c_1_t ;

c_0_ts = ht;
c_1_ts =  1/dx * (circshift(ht,-1,2)-ht);


c_0_ts_t = ht_t;
c_1_ts_t = 1/dx * (circshift(ht_t,-1,2)-ht_t);



%Now we separate the two reconstructions
c_0_ts_h = c_0_ts(1,:);
c_0_ts_hv = c_0_ts(2,:);

c_1_ts_h = c_1_ts(1,:);
c_1_ts_hv = c_1_ts(2,:);


c_0_ts_t_h = c_0_ts_t(1,:);
c_0_ts_t_hv = c_0_ts_t(2,:);

c_1_ts_t_h = c_1_ts_t(1,:);
c_1_ts_t_hv = c_1_ts_t(2,:);


for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU_h   = c_0_ts_h + c_1_ts_h .* diff_x ;
    IU_hv  = c_0_ts_hv + c_1_ts_hv .* diff_x ;
    IUx_h  =  c_1_ts_h ;
    IUx_hv =  c_1_ts_hv ;
    
    IUt_h =  c_0_ts_t_h + c_1_ts_t_h .* diff_x ;
    IUt_hv =  c_0_ts_t_hv + c_1_ts_t_hv .* diff_x ;
    
    R_h = IUt_h +IUx_hv;
    R_hv = IUt_hv + ((IU_h).^(-2)).*(2*IU_hv.*IUx_hv.*IU_h - (IU_hv.^2).*IUx_h) + g*IU_h.*IUx_h;
    
    integral_txq =  (wq(iq)*dist_x_pl.*(R_h).^2) +(wq(iq)*dist_x_pl.*(R_hv).^2);
    integral_txq_h =  (wq(iq)*dist_x_pl.*(R_h).^2);
    integral_txq_hv =  (wq(iq)*dist_x_pl.*(R_hv).^2);

   
    
    L2Rt = L2Rt + sum(integral_txq);
    L2Rt_h = L2Rt_h + sum(integral_txq_h);
    L2Rt_hv = L2Rt_hv + sum(integral_txq_hv);
    L2Rt_arr = L2Rt_arr + (integral_txq);

    
    
end

end