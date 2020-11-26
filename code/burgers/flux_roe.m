function f_h  = flux_roe(a,b,dx_fine)
f_h = (sqrt( .5*a.^2.*(1+sign(a+b)) + .5*b.^2.*(1-sign(a+b)) )) ;
end
