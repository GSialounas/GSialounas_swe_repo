function u_x = fn_3BS(u,c,dx)
u_x = 1./(6*dx)*( 2*circshift(u,-1)  + 3*circshift(u,0) - 6* circshift(u,1) +circshift(u,2));
end