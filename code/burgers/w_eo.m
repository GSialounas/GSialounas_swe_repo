function F = w_eo(a,b,dist_x_pl)
F =(sqrt( .5*(a.^2).*(1+sign(a)) + .5*(b.^2).*(1-sign(b)) )) ;
end
