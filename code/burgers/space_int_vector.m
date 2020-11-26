function int = space_int_vector(x,dist_x_pl,dist_x_min,u,evalt,ex,c_0_coeff_arr) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq

int = 0;
% for j = 1 : length(x)-1
l_coeffs = length(c_0_coeff_arr);



for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
    diff_x = xiq-x;
    for i = 1:size(c_0_coeff_arr,1)
        diff_x_arr(i,:) = [diff_x.^(i-1)];
    end
    IU = sum(c_0_coeff_arr.*diff_x_arr,1);
    int = int + sum(wq(iq).* dist_x_pl .* (IU - ex(xiq,evalt)).^2);
end


end