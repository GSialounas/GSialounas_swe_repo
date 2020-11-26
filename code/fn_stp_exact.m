function step_exact = fn_stp_exact(x,cntr,t,Lx,radius)

u0_exact_stppl = zeros(size(x));
u0_exact_stpmin = zeros(size(x));
u0_exact_stppl(abs(mod((x-cntr-t),Lx))<radius)=1;
u0_exact_stpmin(abs(mod((x-cntr-t),-Lx))<radius)=1;

step_exact =  max(u0_exact_stppl,u0_exact_stpmin);
end