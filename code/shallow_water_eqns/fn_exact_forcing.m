function forcing = fn_exact_forcing(x,t,d,a,omega,k,g)

h =  d+a*sin(omega*t-k*x);
h_t = a*omega*cos(omega*t-k*x);
h_x = -a*k*cos(omega*t-k*x);

v = a*omega*(cosh(k*h)/sinh(k*d)).*sin(omega*t-k*x);
v_t = a*omega*(a*omega*k*cos(omega*t-k*x)).*sinh(k*h)/sinh(k*d).*sin(omega*t-k*x)+...
      a*omega*(cosh(k*h)/sinh(k*d)).*(omega*cos(omega*t-k*x));
  
v_x  = a*omega*(-a*(k^2)*cos(omega*t-k*x).*sinh(k*h)/sinh(k*d)).*sin(omega*t-k*x)+...
     a*omega*(cosh(k*h)/sinh(k*d)).*(-k*cos(omega*t-k*x));
 
 forcing_h = h_t + h_x.*v +v_x.*h;
 forcing_hv = h_t.*v + v_t.*h + h.*(2*v_x.*v)+(v.^2.*h_x)+ g*h.*h_x;
 forcing = [forcing_h; forcing_hv];

end