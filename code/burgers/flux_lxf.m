function F =flux_lxf(dist_x_pl,dt,uold)
f_pl = ( dist_x_pl./(2*dt).*(uold-circshift(uold,-1)) +.5*.5*(uold.^2 +circshift(uold,-1).^2  ));
f_min = ( circshift(dist_x_pl,1)./(2*dt).*(circshift(uold,1)-circshift(uold,0)) +.5*.5*(circshift(uold,1).^2 +circshift(uold,0).^2  ));
F = (1./dist_x_pl) .*(f_pl-f_min);
end