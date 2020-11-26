clc;
clear;
close all;
global nq xq wq
nq = 2;
gauss();

% T=4*pi;
showplots = 0;
save_results_to_cellarr=1;
L_domain= 32*pi;
time_stepping = 'RK3';
rec ='rec_3';
spat_disc = 'LXF';
dt_dx_coupling = 'pow_one';
scheme_string= 'SHW_RK3_WENO3_rec3_fixed';
scheme_arr = {scheme_string};


g=9.81;
d= 2;
lambda = 32*pi;
k = 2*pi/lambda;
omega = sqrt(g*k*tanh(k*d));
a=.1;

fn_exact_soln = @(x,t) fn_exact_solution_shw(x,t,d,a,omega,k,g);
fluxfun='shallow_water'; % select flux function
% Define our Flux function
switch fluxfun
    case 'linear'   % Scalar Advection, CFL_max: 0.65
        c=1; flux = @(w) c*w;
        dflux = @(w) c*ones(size(w));
    case 'burgers' % Burgers, CFL_max: 0.40
        flux = @(w) w.^2/2;
        dflux = @(w) w;
    case 'buckley' % Buckley-Leverett, CFL_max: 0.20 & tEnd: 0.40
        flux = @(w) 4*w.^2./(4*w.^2+(1-w).^2);
        dflux = @(w) 8*w.*(1-w)./(5*w.^2-2*w+1).^2;
    case 'shallow_water'
        % w = [h; hv]
        flux = @(w) [w(2,:); (w(2,:).^2./w(1,:))+.5*g*w(1,:).^2];
        right_vecs = @(h,u) [1,1;u+sqrt(g*h), u-sqrt(g*h)];
        left_vecs = @(h,u) [(-u+sqrt(g*h))./(2*sqrt(g*h)) ,  1./(2*sqrt(g*h)); (u+sqrt(g*h))./(2*sqrt(g*h)) , -1./(2*sqrt(g*h))];
        
        eig_Jf_1 = @(w) w(2,:)./w(1,:) + sqrt(g*w(1,:));
        eig_Jf_2 = @(w) w(2,:)./w(1,:) - sqrt(g*w(1,:));
        dflux = @(w) [0 , 1; -(w(2,:).^2./(w(1,:).^2)) + g*w(1,:), 2*w(2)./w(1)];
end


sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        %         S = @(w) 0.1*w.^2;
        
        %          S = @(w) .1*w.^2;
        %         S =@(x,t) fn_exact_forcing(x,t,d,a,omega,k,g);
    case 'dont'
        S = @(w) zeros(size(w));
end

% scheme_arr = {'SHW_RK1_LXF_rec1_pow_two'};

x = linspace(0,L_domain,1001);
dist_x_pl = x(2:end)- x(1:end-1);
% dt = .05 *(x(2)-x(1));
x = x(1:end-1);


f_size = 14;
N_subplots = 3;
p=0; q=L_domain;


i_exponent = 1;
i_interpolant = 1;
i_scheme =1;
cell_cell_arr_shw = {};

flux_fn  = @(dist_x_pl,dt,Ut) ones(size(Ut));

maxit_arr = [4:7]; % max refinement iterations
axis_p = [0 q -3e-3 3e-3];
dt_arr = L_domain*(2.^(-(maxit_arr+7))); % mesh size
dx_arr = sqrt(10*dt_arr);

for m = 1:length(maxit_arr)
    t= 0;
    
    it =1;
    
    dx = L_domain*(2^(-(maxit_arr(m)+3))); % mesh size
    dx_coarsest = L_domain*(2^(-(maxit_arr(1)+3)));
    dx_finest   = L_domain*(2^(-(maxit_arr(end)+3)));
    dt_coarsest = .1* dx_coarsest;
    T= 1000*dt_coarsest;
    
    % Choice of dt;
    if dt_dx_coupling == "fixed"
        dt=.1*dx_finest;
    elseif dt_dx_coupling=="pow_one"
        dt=.1*dx;
    elseif dt_dx_coupling=="pow_two"
        dt=.1*dx_finest^2;
    else
    end
    x=0:dx:L_domain;
    dist_x_pl = x(2:end)-x(1:end-1);
    dist_x_min = dist_x_pl;
    x = x(1:end-1);
    %     h= 5+cos(2*pi*x);
    %     v=sin(cos(2*pi*x))./h;
    %initial conditions
    h = d+a*sin(k*x);
    %     v= zeros(size(h));
    v = a*omega*(cosh(k*h)/sinh(k*d)).*sin(omega*t-k*x);
    L2L2R = 0; %L2L2 accumulation of ||R||^2
    L2L2R_h = 0; %L2L2 accumulation of ||R_h||^2
    L2L2R_hv = 0; %L2L2 accumulation of ||R_hv||^2
    L2L2R_arr = zeros(size(x));
    
    
    bound_arr_h_eoc  = [0];%zeros(1,ceil(T/dt));
    bound_arr_hv_eoc = [0];%zeros(1,ceil(T/dt));
    bound_arr_eoc = [0];%zeros(1,ceil(T/dt));
    
    error_arr_eoc = [0];% zeros(1,ceil(T/dt));
    time_arr = [0];% zeros(1,ceil(T/dt));
    dofs_arr = [0];% zeros(1,ceil(T/dt));
    EI_index = [0];
    spat_disc="WENO3";
%     figgy_parasite = figure();

    while  t < (T-dt/2)
        %un1 = h and un2 = hv
        % F1= hv and F2 = hv^1+.5*gh^2
        if spat_disc=="WENO3"
            [un1, un2, F1, F2] = dependent(h,v,g);
            h_old = [un1;un2];
            hq_old = h_old;
            f_h_old_weno = WENO3resAdv1d_FD_system([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2);
            
            
            hq_stage1 = hq_old - dt* f_h_old_weno;
            h_stage1 = hq_stage1(1,:);
            v_stage1 = hq_stage1(2,:)./h_stage1;
            [un1_stage1, un2_stage1, F1_stage1, F2_stage1] = dependent(h_stage1,v_stage1,g);
            
            f_h_stage1_weno = WENO3resAdv1d_FD_system([un1_stage1;un2_stage1],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage1,F2_stage1);
            
            % Stage 2
            hq_stage2 = .75* hq_old +.25*hq_stage1 -.25*dt*f_h_stage1_weno;
            h_stage2 = hq_stage2(1,:);
            v_stage2 = hq_stage2(2,:)./h_stage2;
            
            [un1_stage2, un2_stage2, F1_stage2, F2_stage2] = dependent(h_stage2,v_stage2,g);
            f_h_stage2_weno = WENO3resAdv1d_FD_system([un1_stage2;un2_stage2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage2,F2_stage2);
            
            % Stage 3
            hq_stage3 = (1/3)* hq_old +(2/3)*hq_stage2 -(2/3)*dt*f_h_stage2_weno;
            h_stage3 = hq_stage3(1,:);
            v_stage3 = hq_stage3(2,:)./h_stage3;
            
            h_new = [hq_stage3(1,:);hq_stage3(2,:)];
            h = h_new(1,:);
            v = h_new(2,:)./h;
            
            [un1, un2, F1, F2] = dependent(h,v,g);
            
            f_h_new_weno = WENO3resAdv1d_FD_system([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2);
            
            
            
        elseif spat_disc=="WENO5"
            [un1, un2, F1, F2] = dependent(h,v,g);
            h_old = [un1;un2];
            hq_old = h_old;
            f_h_old_weno = WENO5resAdv1d_fdm_gs_system([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);
            
            
            hq_stage1 = hq_old - dt* f_h_old_weno;
            h_stage1 = hq_stage1(1,:);
            v_stage1 = hq_stage1(2,:)./h_stage1;
            [un1_stage1, un2_stage1, F1_stage1, F2_stage1] = dependent(h_stage1,v_stage1,g);
            
            f_h_stage1_weno = WENO5resAdv1d_fdm_gs_system([un1_stage1;un2_stage1],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage1,F2_stage1,x,t);
            
            % Stage 2
            hq_stage2 = .75* hq_old +.25*hq_stage1 -.25*dt*f_h_stage1_weno;
            h_stage2 = hq_stage2(1,:);
            v_stage2 = hq_stage2(2,:)./h_stage2;
            
            [un1_stage2, un2_stage2, F1_stage2, F2_stage2] = dependent(h_stage2,v_stage2,g);
            f_h_stage2_weno = WENO5resAdv1d_fdm_gs_system([un1_stage2;un2_stage2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1_stage2,F2_stage2,x,t);
            
            % Stage 3
            hq_stage3 = (1/3)* hq_old +(2/3)*hq_stage2 -(2/3)*dt*f_h_stage2_weno;
            h_stage3 = hq_stage3(1,:);
            v_stage3 = hq_stage3(2,:)./h_stage3;
            
            h_new = [hq_stage3(1,:);hq_stage3(2,:)];
            h = h_new(1,:);
            v = h_new(2,:)./h;
            
            [un1, un2, F1, F2] = dependent(h,v,g);
            
            f_h_new_weno = WENO5resAdv1d_fdm_gs_system([un1;un2],flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2,x,t);
            
            
            
        else
            
        end
        
        
        
        L2L2R_arr = zeros(size(x));
        for iq = 1 : nq
            tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
            if rec =="rec_1"
                [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv, IUt_h, IUt_hv] = compute_Rs_vector_temp_1_spatiotemp_1_shw(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old,f_h_new,'CS',flux_fn); %compute R at gauss point
            elseif rec =="rec_3"
                % this gives second order until the shock appears but has a
                % central flux
                %                 [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv,IUt_h, IUt_hv] = compute_Rs_vector_temp_3_spatiotemp_3_shw(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn); %compute R at gauss point
                [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, IUt_h, IUx_hv,IUt_hv, IUx_hv_flux ] = compute_Rs_vector_temp_3_spatiotemp_3_shw_weno3(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old_weno,f_h_new_weno,'CS',flux_fn,flux,dflux,S,dx,eig_Jf_1, eig_Jf_2,F1,F2); %compute R at gauss point
                
            else
            end
            L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
            L2L2R_h = L2L2R_h + wq(iq)*dt*(RL2iq_h); %quadrature formula
            L2L2R_hv = L2L2R_hv + wq(iq)*dt*(RL2iq_hv); %quadrature formula
            L2L2R_arr = L2L2R_arr + wq(iq)*dt*(RL2iq_arr); %quadrature formula
            
            if (L2L2R<0)
                disp('L2L2R<0')
            end
        end
        
        
        time_arr(it) = (it-1)*dt;
        bound_arr_eoc(it) = sqrt(L2L2R);
        bound_arr_h_eoc(it) = sqrt(L2L2R_h);
        bound_arr_hv_eoc(it) = sqrt(L2L2R_hv);
        
        % error_arr_eoc(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex));
        error_arr_eoc(it) = 1;%sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,it*dt,ex,c_0_coeff_arr_new));
        dofs_arr(it) = length(x);
        EI_index(it) = bound_arr_eoc(it)./error_arr_eoc(it);
        it = it + 1;
        t = t + dt  ;
        set(gca,'FontSize',16)
        if (mod(it-1,10)==0 &&showplots==1)
                      
            subplot(1,3,1)
            plot(x,h,'b')
            axis([0 L_domain 0 3])
            ylabel('$Height\left(x,t\right)$','interpreter','latex','FontSize',16)
            
            subplot(1,3,2)
            plot(x,h.*v,'b')    
            axis([0 L_domain -.4 .4])

            ylabel('$Momentum\left(x,t\right)$','interpreter','latex','FontSize',16)
            
            
            subplot(1,3,3)
            plot(x,RL2iq_arr,'b')     
            axis([0 L_domain -3e-2 3e-2])
            ylabel('$Local \, Residual$','interpreter','latex','FontSize',16)
%             title(it)
            
            pause(0.1)
            set(gcf, 'Position',  [100, 100, 1200, 400])
            
%             if ((m==1)&&(it==21|| it==621))
%                 saveas(figgy_parasite,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_',scheme_arr{i_scheme},'_',num2str(it),'_dtdx_',num2str(m)],'png')
%             end
        end
        
    end
    if (save_results_to_cellarr==1)
        cell_cell_arr_shw{m}=[time_arr;bound_arr_eoc;error_arr_eoc;EI_index;dofs_arr;bound_arr_h_eoc;bound_arr_hv_eoc];
    end
end
save([scheme_arr{i_scheme},'_cell_arr_file_shw_periodic.mat'],'cell_cell_arr_shw')




%
% function f_h = flux_richtmayer(dist_x_pl,dt,u,F)
% u_pl  = .5 * (circshift(u,-1) + u) - (dt./(2*dist_x_pl)) .* (circshift(F,-1) - F);
% u_min = .5 * (circshift(u,1) + u) - (dt./(2*dist_x_pl)) .* (F-  circshift(F,1));
% f_h = (1./dist_x_pl) .* .5 .* (u_pl.^2  - u_min.^2);
% end
%
%
function f_out = flux_richtmayer_h(dist_x_pl,dt,h, hv,F_h, F_hv,g)
h_pl  = .5 * (circshift(h,-1) + h) - (dt./(2*dist_x_pl)) .* (circshift(F_h,-1) - F_h);
h_min = .5 * (circshift(h,1) + h) - (dt./(2*dist_x_pl)) .* (F_h-  circshift(F_h,1));

hv_pl  = .5 * (circshift(hv,-1) + hv) - (dt./(2*dist_x_pl)) .* (circshift(F_hv,-1) - F_hv);
hv_min = .5 * (circshift(hv,1) + hv) - (dt./(2*dist_x_pl)) .* (F_hv-  circshift(F_hv,1));

f_h_pl = hv_pl;
f_h_min = hv_min;

f_hv_pl = (hv_pl.^2)./(h_pl) +  .5* g * h_pl.^2;
f_hv_min = (hv_min.^2)./(h_min) +  .5* g * h_min.^2;


f_h  = (1./dist_x_pl) .* (f_h_pl-f_h_min);
f_hv = (1./dist_x_pl) .* (f_hv_pl-f_hv_min);

f_out= [f_h; f_hv];
end



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
