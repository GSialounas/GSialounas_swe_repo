clc;
clear;
close all;
% initiate quadrature rule
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();

show_plot = 0; % if show_plot = 0 plots are not shown
Max_refine = 4; % max refinement level in powers of 2
T_end = 2;
% fn_ex = @(x,t) fn_stp_exact(x,.25, t,1, .1);  % periodic step function
fn_ex = @(x,t) sin(2*pi*(x-t));  % periodic step function

init_conds = 'stepIC';
scheme = 'WENO';
cell_cell_arr = {};
for n_ref = 1 : Max_refine
%     dx = 2^(-(n_ref+5));
    c= 1; % advection speed
    x=linspace(0,1,2^(n_ref-1)*100+1);
    x=x(1:end-1);
    dx = x(2)-x(1);
    dt = .1 * dx; % Here the CFL condition is set to .1

    dist_x_pl = dx;
    dist_x_min = dx;
    
    t= 0 ;
    it = 1 ;
    disp('new_cycle')
    u = fn_ex(x,0);
    figure
    plot(x, fn_ex(x,0),'b')
    
    % Create arrays for: error, bound, EOC_error, EOC_bound and Effectivity
    % (EI) index
    time_arr = [];
    time_arr(1) = 0;
    error_arr = [];
    bound_arr = [];
    EOC_error = [];
    EOC_bound = [];
    EI_index  = [];
    
    % inititiate variable for L2 accumulation of bound : 
    L2L2R = 0 ;
    spat_disc = "WENO5";
    
    while t <= (T_end-dt/2)
        % First stage
        u0=u;
        
%         dF = resWENO5(u,c, dx);
        dF = fn_3BS(u,c, dx);

        dFold = dF;
        uold= u ;
        u = u0-dt*dF;
        
        % 2nd Stage
%         dF = resWENO5(u,c,dx);
        dF = fn_3BS(u,c, dx);

        u = 0.75*u0+0.25*(u-dt*dF);
        
        % 3rd stage
%         dF = resWENO5(u,c,dx);
        dF = fn_3BS(u,c, dx);

        u = (u0+2*(u-dt*dF))/3;
        dFnew = dF;
        unew = u;
        
        
        % Compute reconstruction    
        for iq = 1 : nq
            tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
            [RL2iq, c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_spatiotemp_3(x,dist_x_pl,dist_x_min,uold,unew,tq(iq),t,dt, dFold, dFnew,spat_disc); %compute R at gauss point          
            L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
            
            if (L2L2R<0)
                disp('bound accumulation cannot be negative!')
            end
        end
        if it ==1
            c_0_coeff_arr_0 = c_0_coeff_arr_old;
        end
        
        
        if (mod(it-1,20)==0 && show_plot == 1)
            l_coef= length(c_0_coeff_arr_new);
            IU = sum(c_0_coeff_arr_new.*(ones(1,l_coef).*[1;dist_x_pl;dist_x_pl.^2;dist_x_pl.^3]),1);
            plot(x,u,'b',x,IU,'ro')
            pause(0.1)
        end
        
        t = t+dt;
        bound_arr(it) = sqrt(L2L2R*exp(it*dt) + space_int_vector(x,dist_x_pl,dist_x_min,fn_ex(x,0),0,fn_ex,c_0_coeff_arr_0)); % the second part of the expression is the initial error(e_0)        
        error_arr(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,unew,it*dt,fn_ex,c_0_coeff_arr_new));        
        EI_index(it) = bound_arr(it)./error_arr(it);
        
        time_arr(it) = t-dt;
        it= it +1;
        cell_cell_arr{n_ref}=[time_arr; bound_arr; error_arr; EI_index];

        
    end
        
end
save([scheme,'_', init_conds,'_cell_arr_file_step_IC.mat'],'cell_cell_arr')

function [] = gauss
% For n point gauss quadrature, return evaluation points and weights for
% gauss quadrature over [-1,1]
global nq xq wq
if nq == 1
    xq = 0;
    wq = 1;
elseif nq == 2
    xq = [-0.57735026918962576451, 0.57735026918962576451];
    wq = [0.5, 0.5];
elseif nq == 3
    xq = [-0.7745966692414834, 0, 0.7745966692414834];
    wq = 0.5*[0.5555555555555556, 0.8888888888888888, 0.5555555555555556];
elseif nq > 3
    fprintf(1,'No Gauss quadrature implemented of this degree\n')
    return
end
end