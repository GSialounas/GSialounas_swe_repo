function f_h = fn_central_diff(dist_x_pl,dt,u)

f_h = 1./(2*dist_x_pl).*(circshift(u,-1,2)-circshift(u,1,2));
end

