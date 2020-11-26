function f = exact_hat_burger(x)
f = zeros(size(x));
f(abs(x-2)<=1)= 1-abs(x(abs(x-2)<=1)-2);

end