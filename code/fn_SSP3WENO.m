function var_out = fn_SSP3WENO(init_conds, fn_ex)
global nq xq wq % gauss quadrature


showplot = 0; %boolean for graphics
fps= 20;

scheme_arr = {'SSP3WENO3'};
intrpl_arr = {'spl'};
% interpolant_arr
exponent_arr = [1];
% global nq xq wq % gauss quadrature
% nq = 2; %number of quad points in 1d per cell
% gauss();


fluxfun='linear'; % select flux function
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
    
    
    
    scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
    
    T = 2; % final time
    maxit =4; % max refinement iterations
    
    cntr = .25;
    Lx=1;
    pow_exp=100;
    ratio_cTf=2;
    part_fine=.5;
    ex =  fn_ex;%@(x,t) sin(2*pi*(x-t));
    
    ratio_cTf = 2;
    
    for m = 1:maxit
        
        % x = create_grid(Lx, ratio_cTf, h, part_fine); % h is the step in the fine portion of the grid
        x = linspace(0,1,2^(m-1)*100+1);
        x=x(1:end-1);
        h= x(2)-x(1);
        dx=h;
        CFL= .1;
        dt = CFL*h;
        
        dist_x_pl = circshift(x,-1)-x;
        dist_x_pl(end)= Lx-x(end);
        
        dist_x_min = x - circshift(x,1);
        dist_x_min(1)= Lx - x(end);
        
        uold = ex(x,0); % set initial condition
        uold_nu = ex(x,0); % IC for non-uniform mesh
        
        t = 0; %initialise time
        it = 1;
        
        
        % arrays used for EOC for bound and error
        bound_arr = [];%zeros(1,ceil(T/dt));
        error_arr = [];% zeros(1,ceil(T/dt));
        time_arr = [];% zeros(1,ceil(T/dt));
        EI_index = [];
        L2L2R = 0; %L2L2 accumulation of ||R||^2
        
        while t <= (T-dt/2)
            
            
            spat_disc = "WENO3";
            c=1;
            if spat_disc  == "WENO5"
                uold_x_nu  = WENO5resAdv1d_fdm_gs(uold_nu,flux,dflux,S,dx);
                % Three stage
                ustage_0 = uold_nu;
                ustage_1 = uold_nu - dt*uold_x_nu;
                uold_x_stage_1 =  WENO5resAdv1d_fdm_gs(ustage_1,flux,dflux,S,dx);
                ustage_2 = (.75)*uold_nu +.25*ustage_1 - .25 *dt *(uold_x_stage_1);
                uold_x_stage_2 = WENO5resAdv1d_fdm_gs(ustage_2,flux,dflux,S,dx);
                ustage_3 = (1/3)*uold_nu +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
                u_nu = ustage_3;
                
                u_x_nu = WENO5resAdv1d_fdm_gs(u_nu,flux,dflux,S,dx);%resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3); %WENO5resAdv1d_fdm_gs(u_nu,flux,dflux,S,dx);%
                
            elseif spat_disc =="WENO3"
                uold_x_nu  = WENO3resAdv1d(uold_nu,flux,dflux,S,dx);
                % Three stage
                ustage_0 = uold_nu;
                ustage_1 = uold_nu - dt*uold_x_nu;
                uold_x_stage_1 =  WENO3resAdv1d(ustage_1,flux,dflux,S,dx);
                ustage_2 = (.75)*uold_nu +.25*ustage_1 - .25 *dt *(uold_x_stage_1);
                uold_x_stage_2 = WENO3resAdv1d(ustage_2,flux,dflux,S,dx);
                ustage_3 = (1/3)*uold_nu +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
                u_nu = ustage_3;
                u_x_nu = WENO3resAdv1d(u_nu,flux,dflux,S,dx);
                
            elseif spat_disc == "3S"
                uold_x_nu  = FD_backward_3(uold_nu,dx);
                % Three stage
                ustage_0 = uold_nu;
                ustage_1 = uold_nu - dt*uold_x_nu;
                uold_x_stage_1 =  FD_backward_3(ustage_1,dx);
                ustage_2 = (.75)*uold_nu +.25*ustage_1 - .25 *dt *(uold_x_stage_1);
                uold_x_stage_2 = FD_backward_3(ustage_2, dx);
                ustage_3 = (1/3)*uold_nu +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
                u_nu = ustage_3;
                u_x_nu = FD_backward_3(u_nu, dx);
            else
            end
            
            
            
            %Loop over quad points in time
            for iq = 1 : nq
                tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                [RL2iq, c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_rec_space_3_time_3(x,dist_x_pl,dist_x_min,uold_nu,u_nu,tq(iq),t,dt, uold_x_nu, u_x_nu,spat_disc,c); %compute R at gauss point
                L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                if (L2L2R<0)
                    disp('L2L2R cannot be negative')
                end
            end
            if it ==1
                c_0_coeff_arr_0 = c_0_coeff_arr_old;
            end
            
            if (showplot && mod(it-1,fps)==0)
                l_coef= length(c_0_coeff_arr_new);
                IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl;dist_x_pl.^2;dist_x_pl.^3],1);
                plot(x,u_nu,'r*',x,ex(x,t),'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                axis([x(1) x(end) -1 1])
                pause(0.01)
            end
            
            t = t + dt; %move in time
            uold_nu = u_nu;
            time_arr(it) = it*dt;
            bound_arr(it) = sqrt(exp(it*dt) *(L2L2R+ space_int_vector(x,dist_x_pl,dist_x_min,ex(x,0),0,ex,c_0_coeff_arr_0))); % the second part of the expression is the initial error(e_0)
            error_arr(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex,c_0_coeff_arr_new));
            EI_index(it) = bound_arr(it)./error_arr(it);
            it = it+1;
            % Need the error and the bound at every time step for
            % expected order of convergence (EOC) plots
            % error :
            
        end
        
        finalL2err(m) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,T,ex,c_0_coeff_arr_new)); %compute final time L2 error
        R(m) = sqrt(exp(T)*(space_int_vector(x,dist_x_pl,dist_x_min,ex(x,0),0,ex,c_0_coeff_arr_0) + L2L2R)); %bound
        
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
        cell_cell_arr{m}=[time_arr;bound_arr;error_arr;EI_index];
    end
    %     set(gcf, 'Position',  [100, 100, 1000, 500]) % position position width height
    
    save([scheme_arr{i_scheme},'_cell_arr_file_',init_conds,'.mat'],'cell_cell_arr')
end




    function uold_x_nu  = fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu)
        uold_x_nu = ( dist_alpha_coeffs(1,:).*circshift(uold_nu,-1)./dist_alpha_arr(1,:)...
            + dist_alpha_coeffs(2,:).*circshift(uold_nu,0)./dist_alpha_arr(2,:)...
            + dist_alpha_coeffs(3,:).* circshift(uold_nu,1)./dist_alpha_arr(3,:)...
            + dist_alpha_coeffs(4,:).*circshift(uold_nu,2)./dist_alpha_arr(4,:));
    end
var_out =1;
end
