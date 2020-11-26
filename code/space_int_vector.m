function int = space_int_vector(x,dist_x_pl,dist_x_min,u,evalt,ex,c_0_coeff_arr) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
 
int = 0;
l_coeffs = length(c_0_coeff_arr);
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
    diff_x = xiq-x;
    IU = sum(c_0_coeff_arr.*[ones(1,l_coeffs); diff_x; diff_x.^2; diff_x.^3],1);    
    int = int + sum(wq(iq).* dist_x_pl .* (IU - ex(xiq,evalt)).^2);
end

end
