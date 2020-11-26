function hat_fun =fn_hat_exact(x,cntr,t,Lx,h_kmin,h_kpl)

% hat_fun =zeros(size(x));
x_kmin = cntr - h_kmin;
x_kpl  = cntr + h_kpl;
% 
% hat_fun(((x<=cntr) & (x>x_kmin))) = (x((x<=cntr) & (x>x_kmin))-x_kmin)/h_kmin;
% hat_fun(((x>cntr) & (x<x_kpl))) = (x_kpl-x((x>cntr) & (x<x_kpl)))/h_kpl;

hat_fun_minLx = zeros(size(x));
hat_fun_plLx = zeros(size(x));

x_minLx = mod(x-t,-Lx);
x_plLx = mod(x-t,Lx);

hat_fun_minLx(((x_minLx<=cntr) & (x_minLx>x_kmin))) = (x_minLx((x_minLx<=cntr) & (x_minLx>x_kmin))-x_kmin)/h_kmin;
hat_fun_minLx(((x_minLx>cntr) & (x_minLx<x_kpl))) = (x_kpl-x_minLx((x_minLx>cntr) & (x_minLx<x_kpl)))/h_kpl;

hat_fun_plLx(((x_plLx<=cntr) & (x_plLx>x_kmin))) = (x_plLx((x_plLx<=cntr) & (x_plLx>x_kmin))-x_kmin)/h_kmin;
hat_fun_plLx(((x_plLx>cntr) & (x_plLx<x_kpl))) = (x_kpl-x_plLx((x_plLx>cntr) & (x_plLx<x_kpl)))/h_kpl;

hat_fun = max(hat_fun_minLx,hat_fun_plLx);


end