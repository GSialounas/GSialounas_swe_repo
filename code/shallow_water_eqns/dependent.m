function[un1, un2, F1, F2] = dependent(h,v,g)
un1 = h;                 % Function to calculate initial
un2 = h.*v;              % dependent variables.
% un2= v;
F1 = h.*v;
F2 = h.*v.^2 + (g.*h.^2)/2;
% F2 = .5*v.^2 + g*h;

end
