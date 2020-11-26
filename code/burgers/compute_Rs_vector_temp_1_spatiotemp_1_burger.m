
function [L2Rt,L2Rt_arr,c_0_coeff_arr_new, c_0_coeff_arr_old,L2Rt_arr_t,L2Rt_arr_x] = compute_Rs_vector_temp_1_spatiotemp_1_burger(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc,flux_fn)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
L2Rt_arr = zeros(size(x));
L2Rt_arr_t = zeros(size(x));
L2Rt_arr_x = zeros(size(x));
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients

% reconstruction with derivative matching at t_n
c_0_t = uold; % this is u at t_n
c_1_t = (1/dt)*(u-uold);

% reconstruction with derivative matching at t_{n+1}
% c_0_t = uold; % this is u at t_n
% c_1_t = (1/dt) * (2 * u - 2 * uold + dt * u_x);
% c_2_t = -1/(dt^2) * ( u - uold + dt * u_x );


diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t * diff_t ;
Ut_t = c_1_t;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
% h=dist_x_pl(1);
dx_fine = x(2)-x(1);
if spatial_disc == 'CS'
    Ut_x = fn_central_nonu(dist_x_pl, dist_x_min, Ut);
    Ut_xt = fn_central_nonu(dist_x_pl, dist_x_min, Ut_t);
    %     Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
    %     Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
elseif spatial_disc == 'BS'
    Ut_x  = flux_fn(dist_x_pl,dt,Ut);%1./(dist_x_pl).*(flux_lxf(dist_x_pl, Ut, dt)-flux_lxf(dist_x_min,circshift(Ut,1),dt));% fn_backward_burger(dist_x_pl, dist_x_min, Ut);
    Ut_xt = flux_fn(dist_x_pl,dt,Ut_t);%1./dist_x_pl.*(flux_lxf(dist_x_pl, Ut_t, dt)-flux_lxf(dist_x_min,circshift(Ut_t,1),dt));  %fn_backward_burger(dist_x_pl, dist_x_min, Ut_t);
elseif spatial_disc == '2S'
    Ut_x  = fn_2S(dist_x_pl,Ut);
    Ut_xt = fn_2S(dist_x_pl,Ut_t);
    %     Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
    %     Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
elseif spatial_disc == '3S'
    Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2));
    Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
elseif spatial_disc == '4CS'
    Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));
    Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
else
end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
% Reconstruction where we match the spatial derivative at x_j
c_0_ts = Ut; % this is u at t_n
c_1_ts = 1./dist_x_pl.*(circshift(Ut,-1)-circshift(Ut,0));%Ut_x; % - f_h(U^n)

% c_2_ts = (1./(dist_x_pl.^2)) .* (circshift(Ut,-1) - (Ut + dist_x_pl .* Ut_x));

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = 1./dist_x_pl.*(circshift(Ut_t,-1)-circshift(Ut_t,0));%Ut_xt; % - f_h(U^n)
% c_2_ts_t = (1./(dist_x_pl.^2)) .* (circshift(Ut_t,-1) - (Ut_t + dist_x_pl .* Ut_xt));

% reconstruction where we match teh spatial derivative at x_{j+1}
% c_0_ts = Ut; % this is u at t_n
% c_1_ts = (1./dist_x_pl) .* (2*circshift(Ut,-1) - 2 * Ut - dist_x_pl.*circshift(Ut_x,-1)); % - f_h(U^n)
% c_2_ts = -(1./(dist_x_pl.^2)) .* (circshift(Ut,-1) - Ut - dist_x_pl .* circshift(Ut_x,-1));
%
% c_0_ts_t = Ut_t; % this is u at t_n
% c_1_ts_t = (1./dist_x_pl) .* (2*circshift(Ut_t,-1) - 2*Ut_t - dist_x_pl.*circshift(Ut_xt,-1)); % - f_h(U^n)
% c_2_ts_t = -(1./(dist_x_pl.^2)) .* (circshift(Ut_t,-1) - Ut_t - dist_x_pl .* circshift(Ut_xt,-1));


% c_old and c_new with matching derivative at x_j
c_0_old = uold;
c_1_old = 1./dist_x_pl.*(circshift(uold,-1)-uold);%uold_x;
% c_2_old = (1./(dist_x_pl.^2)) .* (circshift(uold,-1) - (uold + dist_x_pl.*uold_x));
% c_2_old = (1./(dist_x_pl.^2)) .* (circshift(uold,-1) - (uold + dist_x_pl.*fn_central_nonu(dist_x_pl,dist_x_min,uold)));



c_0_new = u;
c_1_new = 1./dist_x_pl.*(circshift(u,-1)-u);
% c_2_new = (1./(dist_x_pl.^2)) .* (circshift(u,-1) - (u + dist_x_pl.*u_x));
% c_2_new = (1./(dist_x_pl.^2)) .* (circshift(u,-1) - (u + dist_x_pl.*fn_central_nonu(dist_x_pl,dist_x_min,u)));


% % c_old and c_new with matching derivative at x_{j+1}
% c_0_old = uold;
% c_1_old = (1./dist_x_pl).* (2 * circshift(uold,-1) - 2 * uold - dist_x_pl.*circshift(uold_x,-1) );
% c_2_old = -(1./(dist_x_pl.^2)) .* (circshift(uold,-1) - uold - dist_x_pl .* circshift(uold_x,-1));
%
%
% c_0_new = u;
% c_1_new = (1./dist_x_pl).* (2 * circshift(u,-1) - 2 * u - dist_x_pl.*circshift(u_x,-1) );
% c_2_new = -(1./(dist_x_pl.^2)) .* (circshift(u,-1) - u - dist_x_pl .* circshift(u_x,-1));
%


c_0_coeff_arr_old  =  [c_0_old; c_1_old];
c_0_coeff_arr_new  =  [c_0_new; c_1_new];



% not_include = ones(size(dist_x_pl));
% n_i = find(dist_x_pl~=dist_x_pl(1),1);
% not_include(n_i-1:n_i+1)=0;

%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU  = c_0_ts + c_1_ts .* diff_x ;
    IUx =  c_1_ts ;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x;
    
    integral_txq = (wq(iq)*dist_x_pl.*( (IUt+IU.*IUx)).^2);
    integral_txq_t = (wq(iq)*dist_x_pl.*( IUt).^2);
    integral_txq_x = (wq(iq)*dist_x_pl.*( IUx).^2);
    
    L2Rt = L2Rt + sum(integral_txq);
    L2Rt_arr = L2Rt_arr + (integral_txq);
    L2Rt_arr_t =  c_1_ts;%c_2_ts_t.*diff_x.^2 ;% L2Rt_arr_t + (integral_txq_t);
    L2Rt_arr_x =  c_1_ts_t; %c_2_ts_t.*diff_x.^2 ;%L2Rt_arr_x + (integral_txq_x);
    
    
end
% end

end