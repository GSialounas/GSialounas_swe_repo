function x = create_grid(Lx, ratio_cTf, dx_fine, part_fine)


% In the future add position of fine to coarse relative to fine i.e. fine
% coarse or coarse then fine?
dx_coarse = ratio_cTf * dx_fine;
x_fine = 0 : dx_fine : (part_fine * Lx);
if part_fine ~=1
    x_coarse =  x_fine(end) + dx_coarse : dx_coarse : Lx;
else
    x_coarse= [];
end

x = [x_fine x_coarse];