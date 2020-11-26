function [L2Rt, L2Rt_h, L2Rt_hv,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, IUt_h, IUx_hv, IUt_hv, IUx_hv_flux ] = compute_Rs_vector_temp_3_spatiotemp_3_shw_weno(x,dist_x_pl,dist_x_min,h_old,h_new,evalt,tj,dt,f_h_old,f_h_new,spatial_disc,flux_fn,flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2)
global nq xq wq

g=9.81;
dx= dist_x_pl(1);
L2Rt = 0;
L2Rt_h = 0;
L2Rt_hv = 0;
L2Rt_arr = zeros(size(x));
c_0_coeff_arr_new=0; c_0_coeff_arr_old=0;


fluxfun='linear'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'   % Scalar Advection, CFL_max: 0.65
        c=1; flux = @(w) c*w; 
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40  
        flux = @(w) w.^2/2; 
        dflux = @(w) w; 
    case 'buckley' % Buckley-Leverett, CFL_max: 0.20 & tEnd: 0.40
        flux = @(w) 4*w.^2./(4*w.^2+(1-w).^2);
        dflux = @(w) 8*w.*(1-w)./(5*w.^2-2*w+1).^2;
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        S = @(w) 0.1*w.^2;
    case 'dont'
        S = @(w) zeros(size(w));
end


c_0_t = h_old;
c_1_t = -f_h_old;
c_2_t = (1/dt)*(f_h_new-f_h_old) + (3/(dt^2))*(h_new - h_old + dt*f_h_old);
c_3_t =  -(1/(dt^2))*(f_h_new-f_h_old) - (2/(dt^3))*(h_new - h_old + dt*f_h_old);


diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
ht = c_0_t + c_1_t * diff_t +c_2_t *diff_t^2 + c_3_t*diff_t^3;
ht_t = c_1_t + 2*c_2_t*diff_t +3 * c_3_t*diff_t^2;

c_0_ts = ht;
% c_1_ts = fn_central_diff(dist_x_pl,dt, ht);
c_1_ts =  WENO5resAdv1d_fdm_gs_system_just_for_rec([ht(1,:);ht(2,:)],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2);%fn_central_diff(dist_x_pl,dt, ht);
c_1_ts_weno =  WENO5resAdv1d_fdm_gs_system_just_for_rec([ht(1,:);ht(2,:)],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2);%fn_central_diff(dist_x_pl,dt, ht);

c_2_ts = -(1/(dx)) * (circshift(c_1_ts,-1,2)-c_1_ts) + (3/(dx^2))*(circshift(ht,-1,2) - ht - dx*c_1_ts);
c_3_ts = (1/(dx^2)) * (circshift(c_1_ts,-1,2)-c_1_ts) - (2/(dx^3))*(circshift(ht,-1,2) - ht - dx*c_1_ts);


c_0_ts_t = ht_t;
% c_1_ts_t =  fn_central_diff(dist_x_pl,dt, ht_t);
c_1_ts_t = WENO5resAdv1d_fdm_gs_system_just_for_rec([ht_t(1,:);ht_t(2,:)],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2);% fn_central_diff(dist_x_pl,dt, ht_t);
c_1_ts_t_weno = WENO5resAdv1d_fdm_gs_system_just_for_rec([ht_t(1,:);ht_t(2,:)],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2);% fn_central_diff(dist_x_pl,dt, ht_t);

c_2_ts_t = -(1/(dx)) * (circshift(c_1_ts_t,-1,2)-c_1_ts_t) + (3/(dx^2))*(circshift(ht_t,-1,2) - ht_t - dx*c_1_ts_t);
c_3_ts_t = (1/(dx^2)) * (circshift(c_1_ts_t,-1,2)-c_1_ts_t) - (2/(dx^3))*(circshift(ht_t,-1,2) - ht_t - dx*c_1_ts_t);




%Now we separate the two reconstructions
c_0_ts_h = c_0_ts(1,:);
c_0_ts_hv = c_0_ts(2,:);

c_1_ts_h = c_1_ts(1,:);
c_1_ts_hv = c_1_ts(2,:);

c_2_ts_h = c_2_ts(1,:);
c_2_ts_hv = c_2_ts(2,:);

c_3_ts_h = c_3_ts(1,:);
c_3_ts_hv = c_3_ts(2,:);

c_0_ts_t_h = c_0_ts_t(1,:);
c_0_ts_t_hv = c_0_ts_t(2,:);

c_1_ts_t_h = c_1_ts_t(1,:);
c_1_ts_t_hv = c_1_ts_t(2,:);

c_2_ts_t_h = c_2_ts_t(1,:);
c_2_ts_t_hv = c_2_ts_t(2,:);

c_3_ts_t_h = c_3_ts_t(1,:);
c_3_ts_t_hv = c_3_ts_t(2,:);


for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU_h = c_0_ts_h + c_1_ts_h .* diff_x + c_2_ts_h .* diff_x .^2 +c_3_ts_h.*diff_x.^3;
    IU_hv = c_0_ts_hv + c_1_ts_hv .* diff_x + c_2_ts_hv .* diff_x .^2 +c_3_ts_hv.*diff_x.^3;
    
    IUx_h = ( c_1_ts_h + 2* c_2_ts_h .*diff_x +  3* c_2_ts_h .*diff_x.^2);
    IUx_hv = ( c_1_ts_hv + 2* c_2_ts_hv .*diff_x +  3* c_2_ts_hv .*diff_x.^2);
    
    IUt_h  =  c_0_ts_t_h + c_1_ts_t_h .* diff_x + c_2_ts_t_h .* diff_x .^2 +c_3_ts_t_h.*diff_x.^3;
    IUt_hv =  c_0_ts_t_hv + c_1_ts_t_hv .* diff_x + c_2_ts_t_hv .* diff_x .^2 +c_3_ts_t_hv.*diff_x.^3;
    
    R_h  = IUt_h +IUx_hv;
    IUx_hv_flux = ((IU_h).^(-2)).*(2*IU_hv.*IUx_hv.*IU_h - (IU_hv.^2).*IUx_h) + g*IU_h.*IUx_h;
    R_hv = IUt_hv + ((IU_h).^(-2)).*(2*IU_hv.*IUx_hv.*IU_h - (IU_hv.^2).*IUx_h) + g*IU_h.*IUx_h;
    
    integral_txq = (wq(iq)*dist_x_pl.*(R_h.^2));% + (wq(iq)*dist_x_pl.*(R_hv.^2));
    integral_txq_hv = (wq(iq)*dist_x_pl.*(R_hv.^2)) ;%+ (wq(iq)*dist_x_pl.*(R_hv).^2);
    
    integral_txq_h = (wq(iq)*dist_x_pl.*(R_h.^2));
    
    
    
    
    L2Rt = L2Rt + sum(integral_txq);
    L2Rt_h = L2Rt_h + sum(integral_txq_h);
    L2Rt_hv = L2Rt_hv + sum(integral_txq_hv);
    
    L2Rt_arr = L2Rt_arr + (integral_txq);
    
    
    
end

end