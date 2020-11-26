function u = burger_sol(x,t,N)
u  =zeros(size(x));
for k  = 1 :N
    u  = u + ( besselj(k,k*t)/(k*t)*sin(k*x));
end
u  = -2*u ;
end