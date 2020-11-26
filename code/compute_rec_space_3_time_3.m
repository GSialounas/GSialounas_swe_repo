function [L2Rt,c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_rec_space_3_time_3(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc,c)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element

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


L2Rt = 0;
dx =x(2)- x(1);
%
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients
c_0_t = uold; % this is u at t_n
c_1_t = -uold_x; % - f_h(U^n)
c_2_t = 3/(dt^2) * (u - uold) + (1/dt) * (u_x + 2 * uold_x);
c_3_t = -1/(dt^2) .* (u_x + uold_x) - 2/(dt^3) .* (u - uold);

diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t *diff_t +c_2_t *diff_t^2 + c_3_t*diff_t^3;
Ut_t = c_1_t + 2*c_2_t*diff_t +3*c_3_t*diff_t.^2;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
h=dist_x_pl(1);
if spatial_disc == "CS"
    Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
    Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
elseif spatial_disc == "BS"
    Ut_x = 1./(dist_x_min).*(-circshift(Ut,1) + circshift(Ut,0));
    Ut_xt = 1./(dist_x_min).*(-circshift(Ut_t,1) + circshift(Ut_t,0));
elseif spatial_disc == "2S"
    Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
    Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
elseif spatial_disc == "3S"
        Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2));
        Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
    
%     Ut_x = WENO5resAdv1d_fdm_gs(u,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);
%     Ut_xt = WENO5resAdv1d_fdm_gs(u,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t);
elseif spatial_disc == "4CS"
    Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));
    Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
elseif spatial_disc == "WENO"
    Ut_x =WENO5resAdv1d_fdm_gs(Ut,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
    Ut_xt = WENO5resAdv1d_fdm_gs(Ut_t,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
elseif spatial_disc == "WENO3"
    Ut_x =WENO3resAdv1d(Ut,flux,dflux,S,dx);%resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);%
    Ut_xt = WENO3resAdv1d(Ut_t,flux,dflux,S,dx);%resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t); %
    
else
end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
c_0_ts = Ut; % this is u at t_n
c_1_ts = Ut_x; % - f_h(U^n)
c_2_ts = 3./(dist_x_pl.^2) .* (circshift(Ut,-1) - Ut) - (1./dist_x_pl) .* (circshift(Ut_x,-1) + 2 * Ut_x);
c_3_ts = 1./(dist_x_pl.^2) .* (circshift(Ut_x,-1) + Ut_x) - 2./(dist_x_pl.^3) .* (circshift(Ut,-1) - Ut);

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = Ut_xt; % - f_h(U^n)
c_2_ts_t = 3./(dist_x_pl.^2) .* (circshift(Ut_t,-1) - Ut_t) - (1./dist_x_pl) .* (circshift(Ut_xt,-1) + 2 * Ut_xt);
c_3_ts_t = 1./(dist_x_pl.^2) .* (circshift(Ut_xt,-1) + Ut_xt) - 2./(dist_x_pl.^3) .* (circshift(Ut_t,-1) - Ut_t);



c_0_old = uold;
c_1_old = uold_x;
c_2_old = 3./(dist_x_pl.^2) .* (circshift(uold,-1) - uold) - (1./dist_x_pl).*(circshift(uold_x,-1) + 2 * uold_x);
c_3_old = 1./(dist_x_pl.^2) .* (circshift(uold_x,-1) + uold_x) - 2./(dist_x_pl.^3) .* (circshift(uold,-1) - uold);


c_0_new = u;
c_1_new = u_x;
c_2_new = 3./(dist_x_pl.^2) .* (circshift(u,-1)-u) - (1./dist_x_pl).*(circshift(u_x,-1) + 2 * u_x);
c_3_new = 1./(dist_x_pl.^2) .* (circshift(u_x,-1) + u_x) - (2./dist_x_pl.^3) .* (circshift(u,-1) - u);

c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old; c_3_old];
c_0_coeff_arr_new=  [c_0_new; c_1_new; c_2_new; c_3_new];



%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU = c_0_ts + c_1_ts .* diff_x + c_2_ts .* diff_x .^2 + c_3_ts .* diff_x.^3;
    IUx = c_1_ts + 2* c_2_ts .*diff_x + 3* c_3_ts .*diff_x.^2;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x +  c_2_ts_t.*diff_x.^2 + c_3_ts_t.*diff_x.^3;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUt + IUx).^2);
end
% end

end
