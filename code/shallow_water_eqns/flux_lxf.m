function F_out =flux_lxf(dist_x_pl,dt,uold,F)
f_pl = ( dist_x_pl./(2*dt)).*(uold-circshift(uold,-1)) +.5*(F +circshift(F,-1));
f_min = ( circshift(dist_x_pl,1)./(2*dt)).*(circshift(uold,1)-uold) +.5*(circshift(F,1) + F);

F_out = (1./dist_x_pl) .*(f_pl-f_min);
end
