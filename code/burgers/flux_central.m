function f_h = flux_central(dist_x_pl,dt,u)

f_h = 1./(2*dist_x_pl).*((circshift(u,-1)-circshift(u,1)));
end