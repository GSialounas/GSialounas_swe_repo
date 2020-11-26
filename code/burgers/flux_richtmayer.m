function f_h = flux_richtmayer(dist_x_pl,dt,u)
u_pl  = .5 * (circshift(u,-1) + u) - (dt./(2*dist_x_pl)) .* (.5 * (circshift(u,-1).^2) - .5 * u.^2);
u_min = .5 * (circshift(u,1) + u) - (dt./(2*dist_x_pl)) .* (.5 * u.^2 - .5 * circshift(u,1).^2);
f_h = (1./dist_x_pl) .* .5 .* (u_pl.^2  - u_min.^2);
end