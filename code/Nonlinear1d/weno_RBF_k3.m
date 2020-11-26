function [vm,vp] = weno_RBF_k3(vb)
% RBF-WENO-JS scheme with k=3
% 5th order
% Jingyang Guo and Jae-Hun Jung    May, 2016
global N;
d0 = 3/10;     d1 = 3/5;      d2 = 1/10;
d0tilda = d2;  d1tilda = d1;  d2tilda = d0;

for i = 1:N+3  
    ep2dx2_0 = (vb(i+1)-3*vb(i+2)+3*vb(i+3)-vb(i+4))/(vb(i+1)-15*vb(i+2)+15*vb(i+3)-vb(i+4) +1.e-14);
    ep2dx2_1 = (vb(i)  -3*vb(i+1)+3*vb(i+2)-vb(i+3))/(vb(i)  -15*vb(i+1)+15*vb(i+2)-vb(i+3) +1.e-14);

    root2 = (-2*vb(i+1)+3*vb(i+2)-vb(i+3))/(-vb(i+1)+2*vb(i+2)-vb(i+3) + 1.e-14);

    if and(root2<=3,root2>=0)
      ep2dx2_0 = 0;
      ep2dx2_1 = 0;
    end

    vm0 = (1/3+5/6*ep2dx2_0)*vb(i+2)+(5/6-2/3*ep2dx2_0)*vb(i+3)+(-1/6-1/6*ep2dx2_0)*vb(i+4);
    vp0 = (11/6-9/2*ep2dx2_1)*vb(i+2)+(-7/6+6*ep2dx2_1)*vb(i+3)+(1/3-3/2*ep2dx2_1)*vb(i+4);  
    vm1 = (-1/6-1/6*ep2dx2_0)*vb(i+1)+(5/6-2/3*ep2dx2_0)*vb(i+2)+(1/3+5/6*ep2dx2_0)*vb(i+3);
    vp1 = (1/3+5/6*ep2dx2_1)*vb(i+1)+(5/6-2/3*ep2dx2_1)*vb(i+2)+(-1/6-1/6*ep2dx2_1)*vb(i+3);     
    vm2 = (1/3-3/2*ep2dx2_0)*vb(i)+(-7/6+6*ep2dx2_0)*vb(i+1)+(11/6-9/2*ep2dx2_0)*vb(i+2);
    vp2 = (-1/6-1/6*ep2dx2_1)*vb(i)+(5/6-2/3*ep2dx2_1)*vb(i+1)+(1/3+5/6*ep2dx2_1)*vb(i+2);

    beta0 = 13/12*(vb(i+2)-2*vb(i+3)+vb(i+4))^2 + 1/4*(3*vb(i+2)-4*vb(i+3)+vb(i+4))^2;  
    beta1 = 13/12*(vb(i+1)-2*vb(i+2)+vb(i+3))^2 + 1/4*(vb(i+1)-vb(i+3))^2;
    beta2 = 13/12*(vb(i)-2*vb(i+1)+vb(i+2))^2 + 1/4*(vb(i)-4*vb(i+1)+3*vb(i+2))^2;    

    epsilon = 1.e-6;

    alpha0 = d0/(epsilon+beta0)^2; 
    alpha1 = d1/(epsilon+beta1)^2; 
    alpha2 = d2/(epsilon+beta2)^2;

    alpha0tilda = d0tilda/(epsilon+beta0)^2; 
    alpha1tilda = d1tilda/(epsilon+beta1)^2; 
    alpha2tilda = d2tilda/(epsilon+beta2)^2;

    w0 = alpha0/(alpha0+alpha1+alpha2); 
    w1 = alpha1/(alpha0+alpha1+alpha2); 
    w2 = alpha2/(alpha0+alpha1+alpha2);
    w0tilda = alpha0tilda/(alpha0tilda+alpha1tilda+ alpha2tilda); 
    w1tilda = alpha1tilda/(alpha0tilda+alpha1tilda+ alpha2tilda);
    w2tilda = alpha2tilda/(alpha0tilda+alpha1tilda+ alpha2tilda);

    vm(i) = w0     *vm0+w1     *vm1+w2     *vm2;
    vp(i) = w0tilda*vp0+w1tilda*vp1+w2tilda*vp2;
end
end