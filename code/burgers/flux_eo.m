function f_h = flux_eo(dist_x_pl,dt,u)

f_pl = w_eo(u,circshift(u,-1),dist_x_pl);
f_min = w_eo(circshift(u,1),u,dist_x_pl);
f_h = (1./dist_x_pl).*(f_pl-f_min);
end
