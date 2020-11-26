function output = fn_exact_solution_shw(x,t,d,a,omega,k,g)


h =  d+a*sin(omega*t-k*x);
v = a*omega*(cosh(k*h)/sinh(k*d)).*sin(omega*t-k*x);
output=  [h; h.*v];
end