function f_h_t = flux_richtmayer_partial_temp(dist_x_pl,dt,Ut,Ut_t)
% Temporal derivative of the time discretization for SECOND ORDER
% RECONSTRUCTION
% Ut = c_0 + c_1*(t-t_n) + c_2*(t-t_n)^2; % temporal rec.
% Ut_t = c_1 + 2*c_2*(t-t_n); % temp. rec. time derivative
u_pl  = .5 * (circshift(Ut,-1) + Ut) - (dt./(2*dist_x_pl)) .* (.5 * (circshift(Ut,-1).^2) - .5 * Ut.^2);
u_min = .5 * (circshift(Ut,1) + Ut) - (dt./(2*dist_x_pl)) .* (.5 * Ut.^2 - .5 * circshift(Ut,1).^2);

u_pl_t  =  .5*(circshift(Ut_t,-1) + Ut_t) - (dt./(2*dist_x_pl)) .* ( circshift(Ut.*Ut_t,-1) -Ut.*Ut_t);
u_min_t =  .5*(Ut_t + circshift(Ut_t,1)) - (dt./(2.*circshift(dist_x_pl,1))) .* (Ut.*Ut_t - circshift(Ut.*Ut_t,1));
f_h_t   = (1./dist_x_pl).*(u_pl.*u_pl_t - u_min.*u_min_t);

end