clc;
clear;
close all;
global nq xq wq
nq = 2;
gauss();

% T=4*pi;
showplots = 1;
save_results_to_cellarr=1;
L_domain= 32*pi;
time_stepping = 'RK3';
rec ='rec_3';
spat_disc = 'LXF';
dt_dx_coupling = 'fixed';
scheme_string= 'SHW_RK1_LXF_rec1_fixed';
scheme_arr = {scheme_string};

% scheme_arr = {'SHW_RK1_LXF_rec1_pow_two'};

x = linspace(0,L_domain,1001);
dist_x_pl = x(2:end)- x(1:end-1);
% dt = .05 *(x(2)-x(1));
x = x(1:end-1);

h = .2+.1*sin(1/16*x);
v= zeros(size(h));

g= 9.81;
T = 4;
f_size = 14;
N_subplots = 3;
p=0; q=L_domain;


i_exponent = 1;
i_interpolant = 1;
i_scheme =1;
cell_cell_arr_shw = {};

flux_fn  = @(dist_x_pl,dt,Ut) flux_central(dist_x_pl,dt,Ut);

maxit_arr = [4:7]; % max refinement iterations
axis_p = [0 q -3e-3 3e-3];
dt_arr = L_domain*(2.^(-(maxit_arr+7))); % mesh size
dx_arr = sqrt(10*dt_arr);

for m = 1:length(maxit_arr)
    t= 0;
    
    it =2;
    
    dx = L_domain*(2^(-(maxit_arr(m)+3))); % mesh size
    dx_coarsest = L_domain*(2^(-(maxit_arr(1)+3))); 
    dx_finest   = L_domain*(2^(-(maxit_arr(end)+3))); 
    dt_coarsest = .1* dx_coarsest;
    T= 200*dt_coarsest;
    
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
    h = .2+.1*sin(1/16*x);
    v= zeros(size(h));
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
    
    while  t < (T-dt/2)
        %un1 = h and un2 = hv
        % F1= hv and F2 = hv^1+.5*gh^2
        [un1, un2, F1, F2] = dependent(h,v,g);
        h_old = [un1;un2];
        hq_old = h_old;
        f_h_old = [flux_lxf(dist_x_pl,dt,un1,F1);flux_lxf(dist_x_pl,dt,un2,F2)];
        %         f_h_old = flux_richtmayer_h(dist_x_pl,dt,un1, un2,F1, F2,g);
        
        %         un = [un1;un2] - dt* f_h_old;
        % stage 1:
        if time_stepping == "RK1"
            hq_stage1 = hq_old - dt* f_h_old;
            
            h_new = [hq_stage1(1,:);hq_stage1(2,:)];
            h = h_new(1,:);
            v = h_new(2,:)./h;
            [un1, un2, F1, F2] = dependent(h,v,g);
            f_h_new=[flux_lxf(dist_x_pl,dt,un1,F1);flux_lxf(dist_x_pl,dt,un2,F2)];
            %             f_h_new = flux_richtmayer_h(dist_x_pl,dt,un1, un2,F1, F2,g);
            
        elseif time_stepping =="RK2"
            hq_stage1 = hq_old - dt* f_h_old;
            h_stage1 = hq_stage1(1,:);
            v_stage1 = hq_stage1(2,:)./h_stage1;
            [un1_stage1, un2_stage1, F1_stage1, F2_stage1] = dependent(h_stage1,v_stage1,g);
            f_h_stage1 = [flux_lxf(dist_x_pl,dt,un1_stage1,F1_stage1);flux_lxf(dist_x_pl,dt,un2_stage1,F2_stage1)];
            %             f_h_stage1  =flux_richtmayer_h(dist_x_pl,dt,un1_stage1, un2_stage1, F1_stage1, F2_stage1,g);
            % Stage 2
            hq_stage2 = .5* hq_old +.5*hq_stage1 -.5*dt*f_h_stage1;
            h = hq_stage2(1,:);
            v = hq_stage2(2,:)./h;
            h_new =[hq_stage2(1,:);hq_stage2(2,:)];
            [un1, un2, F1, F2] = dependent(h,v,g);
            
            f_h_new=[flux_lxf(dist_x_pl,dt,un1,F1);flux_lxf(dist_x_pl,dt,un2,F2)];
            %              f_h_new = flux_richtmayer_h(dist_x_pl,dt,un1, un2,F1, F2,g);
        elseif time_stepping == "RK3"
            hq_stage1 = hq_old - dt* f_h_old;
            h_stage1 = hq_stage1(1,:);
            v_stage1 = hq_stage1(2,:)./h_stage1;
            [un1_stage1, un2_stage1, F1_stage1, F2_stage1] = dependent(h_stage1,v_stage1,g);
            f_h_stage1 = [flux_lxf(dist_x_pl,dt,un1_stage1,F1_stage1);flux_lxf(dist_x_pl,dt,un2_stage1,F2_stage1)];
            %              f_h_stage1  =flux_richtmayer_h(dist_x_pl,dt,un1_stage1, un2_stage1, F1_stage1, F2_stage1,g);
            % Stage 2
            hq_stage2 = .75* hq_old +.25*hq_stage1 -.25*dt*f_h_stage1;
            h_stage2 = hq_stage2(1,:);
            v_stage2 = hq_stage2(2,:)./h_stage2;
            
            [un1_stage2, un2_stage2, F1_stage2, F2_stage2] = dependent(h_stage2,v_stage2,g);
            f_h_stage2= [flux_lxf(dist_x_pl,dt,un1_stage2,F1_stage2);flux_lxf(dist_x_pl,dt,un2_stage2,F2_stage2)];
            %             f_h_stage2  =flux_richtmayer_h(dist_x_pl,dt,un1_stage2, un2_stage2, F1_stage2, F2_stage2,g);
            % Stage 3
            hq_stage3 = (1/3)* hq_old +(2/3)*hq_stage2 -(2/3)*dt*f_h_stage2;
            h_stage3 = hq_stage3(1,:);
            v_stage3 = hq_stage3(2,:)./h_stage3;
            
            h_new = [hq_stage3(1,:);hq_stage3(2,:)];
            h = h_new(1,:);
            v = h_new(2,:)./h;
            
            [un1, un2, F1, F2] = dependent(h,v,g);
            
            f_h_new=[flux_lxf(dist_x_pl,dt,un1,F1);flux_lxf(dist_x_pl,dt,un2,F2)];
            %             f_h_new = flux_richtmayer_h(dist_x_pl,dt,un1, un2,F1, F2,g);
            
        else
            
            
        end
        
        
        
        
        for iq = 1 : nq
            tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
            if rec =="rec_1"
                [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv, IUt_h, IUt_hv] = compute_Rs_vector_temp_1_spatiotemp_1_shw(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old,f_h_new,'CS',flux_fn); %compute R at gauss point
            elseif rec =="rec_3"
                [RL2iq,RL2iq_h,RL2iq_hv,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,IU_h, IU_hv, R_h, R_hv,IUx_h, IUx_hv,IUt_h, IUt_hv] = compute_Rs_vector_temp_3_spatiotemp_3_shw(x,dist_x_pl,dist_x_min,h_old,h_new,tq(iq),t,dt,f_h_old,f_h_new,'CS',flux_fn); %compute R at gauss point
                
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
        
        if (mod(it-1,100)==0 &&showplots==1)
            subplot(2,2,1)
            plot(x,h,'b',x, IU_h,'ro')
            axis([0 L_domain 0 .4])
            title('h')
            
            subplot(2,2,2)
            plot(x,IUx_h,'b',x,(1/(2*dx))*(circshift(IU_h,-1)-circshift(IU_h,1)),'ro')
            %             axis([])
            %             title('IU_h')
            title('R_h')
            
            
            subplot(2,2,3)
            plot(x,h.*v,'b', x, IU_hv,'ro')
            axis([0 L_domain -.3 .3])
            title('hv')
            subplot(2,2,4)
            
            plot(x,IUx_hv,'b',x,(1/(2*dx))*(circshift(IU_hv,-1)-circshift(IU_hv,1)),'ro')
            %             axis(axis_p)
            %             title('IU_hv')
            title('IUx_{hv}')
            
            pause(0.01)
            set(gcf, 'Position',  [100, 100, 800, 400])
%             if ((m==length(maxit_arr))&&(reps==20|| reps==1700))
%                 saveas(figgy_parasite,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_',scheme_arr{i_scheme},'_',num2str(reps),'_dtdx_',dt_dx_coupling,'_',num2str(m)],'png')               
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
