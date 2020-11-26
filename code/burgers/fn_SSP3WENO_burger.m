function var_out = fn_SSP3WENO_burger(init_conds,fn_ex)
global nq xq wq % gauss quadrature

showplot = 1;




N_pts_burger = 100;
ex = @ (x,t) burger_sol(x,t,N_pts_burger);
maxit_arr = [4:8]; % max refinement iterations



fps = 20;

U = 1;
K = 0;
T = .5;
scheme_arr = {'SSP3WENO3'};


fluxfun='burgers'; % select flux function
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
end

sourcefun='dont'; % add source term
% Source term
switch sourcefun
    case 'add'
        S = @(w) 0.1*w.^2;
    case 'dont'
        S = @(w) zeros(size(w));
end



for i_scheme = 1:length(scheme_arr)
    
    cell_cell_arr_burger = {};
    
    
    
    scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
    IC_string = 'sin';
    
    for m = 1:length(maxit_arr)
        Lx= 1; ratio_cTf = 1; part_fine = 1;
        h = 2^(-(maxit_arr(m)+3)); % mesh size
        x = create_grid(Lx, ratio_cTf, h, part_fine);
        x2 = -pi +2*pi*(x-x(1))/(x(end)-x(1));
        x=x2(1:end-1);
        dx_fine = x(2)-x(1);
        
        dt = .1*dx_fine^1;
        %                 x=x(1:end-1);
        Lx= pi;
        dist_x_pl = circshift(x,-1)-x;
        dist_x_pl(end)= Lx-x(end);
        
        dist_x_min = x- circshift(x,1);
        dist_x_min(1)= Lx-x(end);
        
        uold = -sin(x);%ex(x,0); % set initial condition
        u= -sin(x);%ex(x,0); % set initial condition
        t = 0; %initialise time
        L2L2R = 0; %L2L2 accumulation of ||R||^2
        exponential_factor = 0;
        
        % arrays used for EOC for bound and error
        bound_arr = [];%zeros(1,ceil(T/dt));
        error_arr = [];% zeros(1,ceil(T/dt));
        time_arr = [];% zeros(1,ceil(T/dt));
        dofs_arr = [];% zeros(1,ceil(T/dt));
        EI_index = [];
        it = 1;
        dx= x(2)-x(1);
        while t<=T-dt/2
            
            %             spat_disc = 'BS';
            %             f_h_old = flux_central(dist_x_pl,dt,uold);%flux_fn(dist_x_pl,dt,uold);
            %             u_stage_1 = uold - dt*uold.*f_h_old;
            %             u = .5 * uold+ .5 * u_stage_1  - .5 * dt * u_stage_1.* flux_central(dist_x_pl,dt,u_stage_1);
            %             f_h_new = flux_central(dist_x_pl,dt,u);%;flux_fn(dist_x_pl,dt,u);%1/dx_fine*(flux_eo(u,circshift(u,-1),dx_fine)- flux_eo(circshift(u,1),circshift(u,0),dx_fine));
            uold= u;
            spat_disc = "WENO3";
            c=1;
            if spat_disc =="WENO5"
                f_h_old  = WENO5resAdv1d_fdm_gs(uold,flux,dflux,S,dx);%resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu); %WENO5resAdv1d_fdm_gs(uold_nu,flux,dflux,S,dx);%
                
                % Three stage
                ustage_0 = uold;
                
                ustage_1 = uold - dt*f_h_old;
                f_h_stage_1 =  WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);%resWENO5(ustage_1,c,dx);% fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_1); %WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);%
                
                ustage_2 = (.75)*uold +.25*ustage_1 - .25 *dt *(f_h_stage_1);
                f_h_stage_2 = WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);%resWENO5(ustage_2,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_2); %WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);%
                
                ustage_3 = (1/3)*uold +(2/3)*ustage_2 - (2/3)*dt*(f_h_stage_2);
                
                u = ustage_3;
                f_h_new = WENO5resAdv1d_fdm_gs(u,flux,dflux,S,dx);%resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3); %WENO5resAdv1d_fdm_gs(u_nu,flux,dflux,S,dx);%
            elseif spat_disc== "WENO3"
                f_h_old  = WENO3resAdv1d(uold,flux,dflux,S,dx);%resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu); %WENO5resAdv1d_fdm_gs(uold_nu,flux,dflux,S,dx);%
                
                % Three stage
                ustage_0 = uold;
                
                ustage_1 = uold - dt*f_h_old;
                f_h_stage_1 =  WENO3resAdv1d(ustage_1,flux,dflux,S,dx);%resWENO5(ustage_1,c,dx);% fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_1); %WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);%
                
                ustage_2 = (.75)*uold +.25*ustage_1 - .25 *dt *(f_h_stage_1);
                f_h_stage_2 = WENO3resAdv1d(ustage_2,flux,dflux,S,dx);%resWENO5(ustage_2,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_2); %WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);%
                
                ustage_3 = (1/3)*uold +(2/3)*ustage_2 - (2/3)*dt*(f_h_stage_2);
                
                u = ustage_3;
                f_h_new = WENO3resAdv1d(u,flux,dflux,S,dx);%resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3); %WENO5resAdv1d_fdm_gs(u_nu,flux,dflux,S,dx);%
                
                
            end
            
            
            % Initial indicator calculation at time step t = it*dt
            L2L2R_arr = zeros(size(x));
            L2L2R_arr_t = zeros(size(x));
            L2L2R_arr_x = zeros(size(x));
            for iq = 1 : nq
                tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                [RL2iq,RL2iq_arr, c_0_coeff_arr_new, c_0_coeff_arr_old,R_h_burger,max_IUx] = compute_burger_rec_space_3_time_3(x,dist_x_pl,dist_x_min,uold,u,tq(iq),t,dt, f_h_old, f_h_new,spat_disc,c); %compute R at gauss point
                L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                L2L2R_arr = L2L2R_arr +wq(iq)*dt*(RL2iq_arr);
                % We will now also make a calculation of the factor that
                % goes into the exponential
                exponential_factor = exponential_factor + wq(iq)*dt*(max_IUx +1); %quadrature formula

                
                if (L2L2R<0)
                    disp('L2L2R<0')
                end
            end
            
            if it ==1
                c_0_coeff_arr_0 = c_0_coeff_arr_old;
                x_0 = x;
                dist_x_pl_0 = dist_x_pl;
                dist_x_min_0 = dist_x_min;
            end
            
            if (showplot && mod(it-1,fps)==0)
                l_coef= length(c_0_coeff_arr_new);
                IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl;dist_x_pl.^2;dist_x_pl.^3],1);
                subplot(1,2,1)
                plot(x,u,'r*',x,ex(x,t),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                axis([x(1) x(end) -1 1])
                subplot(1,2,2)
                plot(x,RL2iq_arr,'b')
                pause(0.01)
            end
            
            t = t +dt;
            uold = u;
            max(L2L2R_arr_x);
            time_arr(it) = it*dt;
            IUx = sum(c_0_coeff_arr_new(2:end, :).*[ones(size(x));2*dist_x_pl; 3*dist_x_pl.^2],1);
            
%             bound_arr(it) = sqrt(exp(it*dt*(1+max(abs(IUx))))*(L2L2R + space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)));
            bound_arr(it) = sqrt(exp(exponential_factor)*(L2L2R + space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex(x_0,0),10^-15,ex,c_0_coeff_arr_0)));

            % error_arr_eoc(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex));
            error_arr(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,it*dt,ex,c_0_coeff_arr_new));%sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,it*dt,ex,c_0_coeff_arr_new));
            dofs_arr(it) = length(x);
            EI_index(it) = bound_arr(it)./error_arr(it);
            it = it + 1;
        end
        finalL2err(m) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u,(it-1)*dt,ex,c_0_coeff_arr_new)); %compute final time L2 error
%         R(m) = sqrt(exp((it-1)*dt*(1+max(abs(IUx))))*(space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex,10^-10,ex,c_0_coeff_arr_0) + L2L2R)); %bound
        R(m) = sqrt(exp(exponential_factor)*(space_int_vector(x_0,dist_x_pl_0,dist_x_min_0,ex,10^-10,ex,c_0_coeff_arr_0) + L2L2R)); %bound

        EOCe(1) = 0;
        EOCR(1) = 0;
        if m > 1
            EOCe(m) = log(finalL2err(m-1)/finalL2err(m))/log(2);
            EOCR(m) = log(R(m-1)/R(m))/log(2);
            EOC_error_arr(m) = EOCe(m);
            EOC_bound_arr(m) = EOCR(m);
        end
        EI(m) = R(m)/finalL2err(m);
        fprintf(1,'||(u - IU)(T)||_L2 = %.5f  EOC = %1.2f\n',finalL2err(m),EOCe(m))
        fprintf(1,'      ||R||_L2(L2) = %.5f  EOC = %1.2f\n',R(m),EOCR(m))
        fprintf(1,'                EI = %.5f  \n',EI(m))
        
        % Arrays for plotting
        cell_cell_arr_burger{m}=[time_arr;bound_arr;error_arr;EI_index;dofs_arr];
    end
    
    
    save([scheme_arr{i_scheme},'_cell_arr_file_sinIC_burger.mat'],'cell_cell_arr_burger')
    
end




var_out =  1;
end