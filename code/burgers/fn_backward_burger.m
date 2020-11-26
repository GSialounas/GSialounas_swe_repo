function f = fn_backward_burger(dist_x_pl, dist_x_min, u)
f = .5./dist_x_min .* (circshift(u.^2,0)-circshift(u.^2,1));
end
