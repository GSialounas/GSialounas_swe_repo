clc;
clear;
close all;

g= 9.81;
depth= 2;
lambda = 32*pi;
k = 2*pi/lambda;
omega = sqrt(g*k*tanh(k*depth));
T = 2*pi/omega;
amp = .1;
x=linspace(0,lambda,1001);
exact_eta=@(x,t) 2+amp*sin(omega*t-k*x);
exact_u =@(x,t) (amp*omega*cosh(k*(depth+amp*sin(omega*t-k*x))))/sinh(k*depth).*sin(omega*t-k*x);
dt = T/1000;

figure
for it = 1:10000
    t = it*dt;
    if mod(it,10)==0
        subplot(1,2,1)
        plot(x, exact_eta(x,t),'bo')
        axis([0 lambda 0 3])
        subplot(1,2,2)
        plot(x,exact_u(x,t),'ro')
        axis([0 lambda -5 5])

        pause(0.01)
    end
end