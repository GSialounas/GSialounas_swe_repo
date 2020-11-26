clear all
close all
showplot = 0; %boolean for graphics
fps= 20;

scheme_arr = {'SSP3WENO'};
intrpl_arr = {'spl'};
% interpolant_arr
exponent_arr = [1];
global nq xq wq % gauss quadrature
nq = 2; %number of quad points in 1d per cell
gauss();

for i_scheme = 1:length(scheme_arr)
    % This cell arrays have three levels:
    % [exponent][interpolant][refinment level]
    % At the last level, the elements are arrays wherein the first row is
    % the time instant, the second row is the bound at that time instant
    % and the last row is the error at that time instant.
    cell_cell_arr = {};
    
    figgy = figure;
    for i_exponent = 1:length(exponent_arr)
        file_name = [scheme_arr{i_scheme},'_','sin_IC','Err_L2_','dt~dx',num2str(exponent_arr(i_exponent)),'.txt'];
        file_name_eoc = [scheme_arr{i_scheme},'_','sin_IC','Err_L2_','dt~dx',num2str(exponent_arr(i_exponent)),'_eoc.txt'];
        
        fileID = fopen(file_name,'w');
        fileID_eoc = fopen(file_name_eoc,'w');
        
        
        cell_arr_err_interpl = {};
        cell_arr_bound_interpl = {};
        cell_arr_eoc_err_interpl = {};
        cell_arr_eoc_bound_interpl = {};
        cell_arr_time = {};
        
        for i_interpolant = 1:length(intrpl_arr)
            
            
            % 'SSP3'   order 3 time, order 2 space central
            % 'SSP33S' order 3 time, order 3 space upwinding
            
            % 'pwl' the standard piecewise bilinear lagrange interpolant
            % 'qdr' quadratic (in space) interpolant
            % 'spl' cubic spline (in space)
            
            scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
            interpolant = intrpl_arr{i_interpolant};
            exponent =exponent_arr(i_exponent);
            IC_string = 'sin';
            
            T = 2; % final time
            maxit =5; % max refinement iterations
            
            cntr = .25;
            Lx=1;
            pow_exp=100;
            ratio_cTf=2;
            part_fine=.5;
            % ex = @(x,t) exp(-100*(x-cntr - t).^2); %smooth exact solution
%             ex = @(x,t) max(exp(-pow_exp*((mod(x-cntr-t,Lx)).^2)),exp(-pow_exp*((mod(x-cntr-t,-Lx)).^2))); %smooth exact solution
            ex = @(x,t) sin(2*pi*(x-t));

            %ex = @(x,t) (x-t > 0.25).*(x-t <= 0.75); %discontinuous exact solution
            %ex = @(x,t) (x-t > 0.25).*(x-t <= 0.75) .* (1-2*abs(2*(x-t)-1)); %C0 exact solution
            
            
            fprintf(1,'Using the %s scheme\n',scheme)
            
            ratio_cTf = 2;
            time1  = zeros(1,maxit);
            time2  = zeros(1,maxit);
            dx_arr = zeros(1,maxit);
            dt_arr = zeros(1,maxit);
            
            error_arr = zeros(1,maxit);
            bound_arr = zeros(1,maxit);
            EOC_error_arr = zeros(1,maxit);
            EOC_bound_arr = zeros(1,maxit);
            
            for m = 1:maxit
                if scheme_arr{i_scheme} == "FTBS"
                    m_plus = 3;
                else
                    m_plus =3;
                end
%                 h = 2^(-(m+m_plus)); % mesh size
                h = 2^(-(m+m_plus)); % mesh size

                dt = 0.1*h^exponent; % timestep size
                dx_arr(m) = h;
                dt_arr(m) = dt; % timestep size
                
                %     x = 0:h:(Lx-h); % spatial partition for uniform mesh
%                 x = create_grid(Lx, ratio_cTf, h, part_fine); % h is the step in the fine portion of the grid
                x = linspace(0,1,2^(m-1)*100+1);
                x=x(1:end-1);
                h= x(2)-x(1);
                dx=h;
                dt = .1*h;
                dist_x_pl = circshift(x,-1)-x;
                dist_x_pl(end)= Lx-x(end);
                
                dist_alpha_arr = [dist_x_pl; circshift(dist_x_pl,1);circshift(dist_x_pl,1); circshift(dist_x_pl,1) + circshift(dist_x_pl,2)];
                dist_alpha_coeffs = zeros(size(dist_alpha_arr));
                
                dist_x_min = x - circshift(x,1);
                dist_x_min(1)= Lx - x(end);
                
                uold = ex(x,0); % set initial condition
                uold_nu = ex(x,0); % IC for non-uniform mesh
                
                % Find the coefficients for the third_order scheme
                a_i_pl1 = zeros(size(x));
                a_i = zeros(size(x));
                a_i_min1 = zeros(size(x));
                a_i_min2 = zeros(size(x));
                
                for j = 1 :length(x)
                    alpha_mat = [1./dist_alpha_arr(:,j)';...
                        1 0 -1 -1; ...
                        [1,0,1,1].*dist_alpha_arr(:,j)'/2;  ...
                        [1,0,-1,-1].*(dist_alpha_arr(:,j).^2)'/6];
                    dist_alpha_coeffs(:,j) = alpha_mat\[0;1;0;0];
                end
                
                %     if showplot
                %         fig = figure();
                %         lnh = plot(x,uold); %approximate in blue
                %         hold on
                %         exh = plot(x,ex(x,0),'ro'); %exact in red
                %     end
                t = 0; %initialise time
                i = 0;
                L2L2R = 0; %L2L2 accumulation of ||R||^2
                tic;
                reps = 0;
                
                % arrays used for EOC for bound and error
                bound_arr_eoc = [];%zeros(1,ceil(T/dt));
                error_arr_eoc = [];% zeros(1,ceil(T/dt));
                time_arr = [];% zeros(1,ceil(T/dt));
                EI_index = [];
                it = 1;
                while t <= (T-dt/2)
                    tic;
                    
                    i = i + 1;
                    if scheme == "SSP3CS"
                        spat_disc = "4CS";
%                         uold_x_nu = 1./((dist_x_pl+dist_x_min)).*(-circshift(uold_nu,1) + circshift(uold_nu,-1));
                        uold_x_nu = -1./(12*h).*(-circshift(uold_nu,2)+8*circshift(uold_nu,1) - 8*circshift(uold_nu,-1) + circshift(uold_nu,-2));

                        if it ==1 
                            u_x_0 = uold_x_nu;
                            u0 = uold_nu;
                        end
                        % We don't actually use stage_0, which is the same as uold_nu i.e. ustage_0 = uold_nu;
                        
                        ustage_1 = uold_nu - dt*uold_x_nu;
%                         uold_x_stage_1 = 1./((dist_x_pl+dist_x_min)).*(-circshift(ustage_1,1) + circshift(ustage_1,-1));
                        uold_x_stage_1 = -1./(12*h).*(-circshift(ustage_1,2)+8*circshift(ustage_1,1) - 8*circshift(ustage_1,-1) + circshift(ustage_1,-2));
                        
                        ustage_2 = (.75)*uold_nu + .25*ustage_1 - .25 *dt *(uold_x_stage_1);
%                         uold_x_stage_2 = 1./((dist_x_pl+dist_x_min)).*(-circshift(ustage_2,1) + circshift(ustage_2,-1));
                        uold_x_stage_2 = -1./(12*h).*(-circshift(ustage_2,2)+8*circshift(ustage_2,1) - 8*circshift(ustage_2,-1) + circshift(ustage_2,-2));

                        ustage_3 = (1/3)*uold_nu +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
                        u_nu = ustage_3;
                        
%                         u_x_nu = 1./((dist_x_pl+dist_x_min)).*(-circshift(u_nu,1) + circshift(u_nu,-1));
                        u_x_nu = -1./(12*h).*(-circshift(u_nu,2)+8*circshift(u_nu,1) - 8*circshift(u_nu,-1) + circshift(u_nu,-2));

                   
                    elseif scheme == 'SSP3WENO'
                        spat_disc = '3S';
                        c=1;
%                         uold_x_nu = 1./(6*h)*( 2*circshift(uold_nu,-1)  + 3*circshift(uold_nu,0) - 6* circshift(uold_nu,1) +circshift(uold_nu,2));
                        uold_x_nu  = resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu);
%                         uold_x_nu = ( dist_alpha_coeffs(1,:).*circshift(uold_nu,-1)./dist_alpha_arr(1,:)...
%                             + dist_alpha_coeffs(2,:).*circshift(uold_nu,0)./dist_alpha_arr(2,:)...
%                             + dist_alpha_coeffs(3,:).* circshift(uold_nu,1)./dist_alpha_arr(3,:)...
%                             + dist_alpha_coeffs(4,:).*circshift(uold_nu,2)./dist_alpha_arr(4,:));
                            
                        if it ==1
                            u_x_0 = uold_x_nu;
                            u0 = uold_nu;
                        end
%                         ustage_0 = uold_nu;
                        ustage_0 = uold_nu;

%                         ustage_1 = uold_nu - dt*uold_x_nu;
                        ustage_1 = uold_nu - dt*uold_x_nu;
%                         uold_x_stage_1 = 1./(6*h)*( 2*circshift(ustage_1,-1)  + 3*circshift(ustage_1,0) - 6* circshift(ustage_1,1) +circshift(ustage_1,2));
                        uold_x_stage_1 = resWENO5(ustage_1,c,dx);% fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_1);
                        %                         uold_x_stage_1= ( dist_alpha_coeffs(1,:).*circshift(ustage_1,-1)./dist_alpha_arr(1,:)...
%                             + dist_alpha_coeffs(2,:).*circshift(ustage_1,0)./dist_alpha_arr(2,:)...
%                             + dist_alpha_coeffs(3,:).* circshift(ustage_1,1)./dist_alpha_arr(3,:)...
%                             + dist_alpha_coeffs(4,:).*circshift(ustage_1,2)./dist_alpha_arr(4,:));
                        
                        ustage_2 = (.75)*uold_nu +.25*ustage_1 - .25 *dt *(uold_x_stage_1);
%                         uold_x_stage_2 = 1./(6*h)*( 2*circshift(ustage_2,-1)  + 3*circshift(ustage_2,0) - 6* circshift(ustage_2,1) +circshift(ustage_2,2));
                        uold_x_stage_2 = resWENO5(ustage_2,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_2);
%                         uold_x_stage_2= ( dist_alpha_coeffs(1,:).*circshift(ustage_2,-1)./dist_alpha_arr(1,:)...
%                             + dist_alpha_coeffs(2,:).*circshift(ustage_2,0)./dist_alpha_arr(2,:)...
%                             + dist_alpha_coeffs(3,:).* circshift(ustage_2,1)./dist_alpha_arr(3,:)...
%                             + dist_alpha_coeffs(4,:).*circshift(ustage_2,2)./dist_alpha_arr(4,:));
%                         
                        
                        ustage_3 = (1/3)*uold_nu +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
                        
                        u_nu = ustage_3;
%                         u_x_nu = ( dist_alpha_coeffs(1,:).*circshift(u_nu,-1)./dist_alpha_arr(1,:)...
%                             + dist_alpha_coeffs(2,:).*circshift(u_nu,0)./dist_alpha_arr(2,:)...
%                             + dist_alpha_coeffs(3,:).* circshift(u_nu,1)./dist_alpha_arr(3,:)...
%                             + dist_alpha_coeffs(4,:).*circshift(u_nu,2)./dist_alpha_arr(4,:));
%                         u_x_nu = 1./(6*h)*( 2*circshift(u_nu,-1)  + 3*circshift(u_nu,0) - 6* circshift(u_nu,1) +circshift(u_nu,2));
                        u_x_nu = resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3);
                    else
                        fprintf(1,'Invalid scheme choice');
                        break
                    end
                    time1inc =toc;
                    time1(m) = time1(m) +time1inc;
                   
                    
                    %Loop over quad points in time
                    for iq = 1 : nq
                        tq(iq) = 0.5*dt*xq(iq) + t + dt/2; %iq-th temporal gauss point on [ti,ti+1]
                        
                        if interpolant == 'pwl'
                            %                 RL2iq = compute_R(x,uold,u,tq(iq),t,dt); %compute R at gauss point
                            RL2iq = compute_R_vector(x,dist_x_pl,dist_x_min,uold_nu,u_nu,tq(iq),t,dt); %compute R at gauss point
                        elseif interpolant == 'qdr'
                            %                 disp('crap1')
                            %                 RL2iq = compute_Rs(x,uold,u,tq(iq),t,dt); %compute R at gauss point
%                             RL2iq = compute_Rs_vector(x,dist_x_pl,dist_x_min,uold_nu,u_nu,tq(iq),t,dt, uold_x_nu, u_x_nu); %compute R at gauss point
                            [RL2iq, c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold_nu,u_nu,tq(iq),t,dt, uold_x_nu, u_x_nu,spat_disc);
                        elseif interpolant == 'spl'
%                              [RL2iq, c_0_coeff_arr] = compute_Rs3_vector(x,dist_x_pl,dist_x_min,uold_nu,u_nu,tq(iq),t,dt, uold_x_nu, u_x_nu); %compute R at gauss point
                             [RL2iq, c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs3_vector_temp_3_spatiotemp_3(x,dist_x_pl,dist_x_min,uold_nu,u_nu,tq(iq),t,dt, uold_x_nu, u_x_nu,spat_disc,dist_alpha_coeffs, dist_alpha_arr); %compute R at gauss point
                            
                        else
                            fprintf(1,'Invalid interpolant')
                            break
                        end
                        
                        L2L2R = L2L2R + wq(iq)*dt*(RL2iq); %quadrature formula
                        if (L2L2R<0)
                            disp('shit')
                        end
                    end
                    if it ==1
                        c_0_coeff_arr_0 = c_0_coeff_arr_old;
                    end
                    
                     if (showplot && mod(reps,fps)==0)
                        l_coef= length(c_0_coeff_arr_new);
                        IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl;dist_x_pl.^2;dist_x_pl.^3],1);
%                         IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl;dist_x_pl.^2],1);

                        %             set(lnh,'XData',x,'YData',u)
                        %             set(exh,'XData',x,'YData',ex(x,t+dt))
                        %             axis([x(1) x(end) -1e-3 1e-3])
                        %             drawnow
%                         plot(x,u_nu,'b') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
                        plot(x,u_nu,'b',x,IU,'r*') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')

                        axis([x(1) x(end) -1 1])
                        pause(0.01)
                    end
                    time2inc = toc-time1inc;
                    time2(m) = time2(m) +time2inc;
                    t = t + dt; %move in time
%                     uold = u;
                    uold_nu = u_nu;
                    reps = reps+1;
                    time_arr(it) = it*dt;
                    bound_arr_eoc(it) = sqrt(exp(it*dt) *(L2L2R+ space_int_vector(x,dist_x_pl,dist_x_min,ex(x,0),0,ex,c_0_coeff_arr_0))); % the second part of the expression is the initial error(e_0)
%                     bound_arr_eoc(it) = sqrt(L2L2R*exp(it*dt) ); % the second part of the expression is the initial error(e_0)

                    error_arr_eoc(it) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex,c_0_coeff_arr_new));
%                     error_arr_eoc_not_gp(it) = sqrt(space_int_vector_not_gp(x,dist_x_pl,dist_x_min,u_nu,it*dt,ex));

                    EI_index(it) = bound_arr_eoc(it)./error_arr_eoc(it);
                    it = it+1;
                    % Need the error and the bound at every time step for
                    % expected order of convergence (EOC) plots
                    % error :
                    
                end
                %     finalL2err(m) = sqrt(space_int(x,u,T,ex)); %compute final time L2 error
                %     R(m) = sqrt(exp(T)*(space_int(x,ex(x,0),0,ex) + L2L2R)); %bound
%                 finalL2err(m) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,T,ex)); %compute final time L2 error
                
                
                finalL2err(m) = sqrt(space_int_vector(x,dist_x_pl,dist_x_min,u_nu,T,ex,c_0_coeff_arr_new)); %compute final time L2 error
                R(m) = sqrt(exp(T)*(space_int_vector(x,dist_x_pl,dist_x_min,ex(x,0),0,ex,c_0_coeff_arr_0) + L2L2R)); %bound

%                 R(m) = sqrt(exp(T)*( L2L2R)); %bound
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
                
                error_arr(m) = finalL2err(m);
                bound_arr(m) = R(m);
                cell_cell_arr{i_exponent}{i_interpolant}{m}=[time_arr;bound_arr_eoc;error_arr_eoc;EI_index];
            end
            
            cell_arr_err_interpl{i_interpolant} = error_arr;
            cell_arr_bound_interpl{i_interpolant} = bound_arr;
            cell_arr_eoc_err_interpl{i_interpolant} = EOC_error_arr;
            cell_arr_eoc_bound_interpl{i_interpolant} = EOC_bound_arr;
            
            
            
        end
        pow_arr = exponent;
        subplot(1,length(exponent_arr),i_exponent)
%         errs = plot_errors_cumulative_Tristan(dx_arr, dt_arr, cell_arr_err_interpl, cell_arr_bound_interpl, cell_arr_eoc_err_interpl, cell_arr_eoc_bound_interpl, scheme, IC_string,pow_arr,intrpl_arr);
        % The next few lines create create a text file for every scheme,
        % for every exponent dt=dx^exponent.  The text file has a table
        % with the following form
        %
        % 1)dx_array, 2)pwl_bound, 3)qdr_bound, 4)cubic_bound, 5)error
        %
        % The number of lines in the text file is the same as the number of
        % spatial step sizes used (maxits)
        A = [dx_arr'];
        for i_A  = 1 : length(cell_arr_bound_interpl)
            A=[A,cell_arr_bound_interpl{i_A}'];
        end
        % For the error, we take the final column.  It would make no
        % difference if we took any of the others as they are the same,
        % because nothing essentially changes in the error computation:
        % only in the bound computation.
        A=[A, cell_arr_err_interpl{end}'];
        
        
        for ii = 1: size(A,1)
            fprintf(fileID,'%6.7f, %6.7f, %6.7f, %6.7f, %6.7f\n',A(ii,:));
        end
        fclose(fileID);
        
        % Do the EOC as well
        A_eoc = [dx_arr'];
        for i_A  = 1 : length(cell_arr_eoc_bound_interpl)
            A_eoc=[A_eoc,cell_arr_eoc_bound_interpl{i_A}'];
        end
        % For the error, we take the final column.  It would make no
        % difference if we took any of the others as they are the same,
        % because nothing essentially changes in the error computation:
        % only in the bound computation.
        A_eoc=[A_eoc, cell_arr_eoc_err_interpl{end}'];
        
        
        for ii = 1: size(A_eoc,1)
            fprintf(fileID_eoc,'%6.7f, %6.7f, %6.7f, %6.7f, %6.7f\n',A_eoc(ii,:));
        end
        
        fclose(fileID_eoc);
        %         Mread = dlmread(file_name);
    end
    set(gcf, 'Position',  [100, 100, 1000, 500]) % position position width height
    
%     saveas(figgy,['fig_',scheme_arr{i_scheme},' L2_err_bound_',IC_string,' ','_Tristan'],'png')
    save([scheme_arr{i_scheme},'_cell_arr_file_sin_IC.mat'],'cell_cell_arr')
end

function L2Rt = compute_Rs(x,uold,u,evalt,tj,dt,uold_x,u_x)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
h = x(2) - x(1);
%%%%%%%%%%%%%%%
% Each closed [x_{i-1},x_i] sub_interval has a piece-wise cubic Hermite
% spline interpolant in it.  This will be defined by four coeffiecients:
% IU(x) = c_0 +c_1(x-x_{i-1})^1 + c_2(x-x_{i-1})^2 +c_3(x-x_{i-1})^3 for
% x in [x_{i-1},x_i].   The way we will be figuring out this coefficients
% is by prescribing
% IU(x_i-1) = U(x_i-1), IU_(x_i) = U(x_i) , and for the derivatives of the
% interpolant we do IU_x(x_i-1) = (U_{x_i-2}-U_{x_i})/2h and
% IU_x(x_i) = (U_{x_i-1}-U_{x_i+1})/2h
% The coefficients of the interpolant can be calculated explicitly from the
% values of u and the difference quotient at the end points of the interval
% (see suli introdcution to numerical analysis page
%%%%%%%%%%%%%%%
% 
% uold_x = (circshift(uold,-1) - circshift(uold,1))/(2*h);
% u_x = (circshift(u,-1) - circshift(u,1))/(2*h);

% u = spline(x,u);
% [breaks,coefs,l,k,d] = unmkpp(u);
% du = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%
% uold = spline(x,uold);
% [breaks,coefs,l,k,d] = unmkpp(uold);
% duold = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%xx=0:h/100:1;
%yy = ppval(u, xx);
%plot(xx,yy)
% IUt_arr = (u - uold)/dt; %defined as a constant for any t and piece-wise quadratic in x
% IU = (evalt - tj)*u/dt + (dt+tj-evalt)*uold/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt

for j = 1 : length(x)-1
    % the values of the coefficients will be calculated on sub-interval by
    % sub-interval basis and are taken from suli page 301
    c_0_old = uold(j);
    c_1_old = uold_x(j);
    c_2_old = (uold(j+1) - (uold(j) + uold_x(j) * h))/(h^2);
    
    c_0_new = u(j);
    c_1_new = u_x(j);
    c_2_new = (u(j+1) - (u(j) + u_x(j) * h))/(h^2);
    
    xl= zeros(size(nq));
    for iq = 1:nq
        xl(iq) = 0.5*h*xq(iq) + x(j) + h/2;
        diff_x = xl(iq)-x(j);
        IUt = 1/dt * (c_0_new + c_1_new * diff_x + c_2_new * diff_x^2) - 1/dt * (c_0_old + c_1_old * diff_x + c_2_old * diff_x^2);
        %         IUt = (ppval(u,xl(iq)) - ppval(uold,xl(iq)))/dt; %defined as a constant for any t
        %         IUx = (evalt - tj)/dt * (ppval(duold,xl(iq))) + (dt+tj-evalt)/dt * (ppval(du,xl(iq)))/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
        IUx = (evalt - tj)/dt * (c_1_new + 2 * c_2_new * diff_x) + (dt+tj-evalt)/dt * (c_1_old + 2 * c_2_old * diff_x); %Lagrange interpolant of IU over the time interval evaluated at evalt
        
        %IUt(iq) = (xl(iq) - x(j))*IUt(j+1)/h + (x(j+1)-xl(iq))*IUt(j)/h;
        
        L2Rt = L2Rt + wq(iq)*h*(IUt + IUx)^2;
    end
end

end



function L2Rt = compute_Rs_vector(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
% h = x(2) - x(1);
%%%%%%%%%%%%%%%
% Each closed [x_{i-1},x_i] sub_interval has a piece-wise cubic Hermite
% spline interpolant in it.  This will be defined by four coeffiecients:
% IU(x) = c_0 +c_1(x-x_{i-1})^1 + c_2(x-x_{i-1})^2 +c_3(x-x_{i-1})^3 for
% x in [x_{i-1},x_i].   The way we will be figuring out this coefficients
% is by prescribing
% IU(x_i-1) = U(x_i-1), IU_(x_i) = U(x_i) , and for the derivatives of the
% interpolant we do IU_x(x_i-1) = (U_{x_i-2}-U_{x_i})/2h and
% IU_x(x_i) = (U_{x_i-1}-U_{x_i+1})/2h
% The coefficients of the interpolant can be calculated explicitly from the
% values of u and the difference quotient at the end points of the interval
% (see suli introdcution to numerical analysis page
%%%%%%%%%%%%%%%
% uold_x = (circshift(uold,-1) - circshift(uold,1))./(dist_x_pl+dist_x_min);
% u_x = (circshift(u,-1) - circshift(u,1))./(dist_x_pl+dist_x_min);

% u = spline(x,u);
% [breaks,coefs,l,k,d] = unmkpp(u);
% du = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%
% uold = spline(x,uold);
% [breaks,coefs,l,k,d] = unmkpp(uold);
% duold = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%xx=0:h/100:1;
%yy = ppval(u, xx);
%plot(xx,yy)
% IUt_arr = (u - uold)/dt; %defined as a constant for any t and piece-wise quadratic in x
% IU = (evalt - tj)*u/dt + (dt+tj-evalt)*uold/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt

% for j = 1 : length(x)-1
% the values of the coefficients will be calculated on sub-interval by
% sub-interval basis and are taken from suli page 301
c_0_old = uold;
c_1_old = uold_x;
c_2_old = (circshift(uold,-1) - (uold + uold_x .* dist_x_pl))./(dist_x_pl.^2);

c_0_new = u;
c_1_new = u_x;
c_2_new = (circshift(u,-1) - (u + u_x .* dist_x_pl))./(dist_x_pl.^2);

%     xl= zeros(size(nq));
for iq = 1:nq
    %         xl(iq) = 0.5*h*xq(iq) + x(j) + h/2;
    xiq = 0.5*dist_x_pl*xq(iq) + x + dist_x_pl/2;
    
    diff_x = xiq-x;
    IUt = 1/dt * (c_0_new + c_1_new .* diff_x + c_2_new .* diff_x.^2) - 1/dt * (c_0_old + c_1_old .* diff_x + c_2_old .* diff_x.^2);
    %         IUt = (ppval(u,xl(iq)) - ppval(uold,xl(iq)))/dt; %defined as a constant for any t
    %         IUx = (evalt - tj)/dt * (ppval(duold,xl(iq))) + (dt+tj-evalt)/dt * (ppval(du,xl(iq)))/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
    IUx = (evalt - tj)/dt * (c_1_new + 2 * c_2_new .* diff_x) + (dt+tj-evalt)/dt * (c_1_old + 2 * c_2_old .* diff_x); %Lagrange interpolant of IU over the time interval evaluated at evalt
    
    %IUt(iq) = (xl(iq) - x(j))*IUt(j+1)/h + (x(j+1)-xl(iq))*IUt(j)/h;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUt + IUx).^2);
    %         disp('crapvec')
end
% end

end


function L2Rt = compute_Rs3(x,uold,u,evalt,tj,dt,uold_x,u_x)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
h = x(2) - x(1);
%%%%%%%%%%%%%%%
% Each closed [x_{i-1},x_i] sub_interval has a piece-wise cubic Hermite
% spline interpolant in it.  This will be defined by four coeffiecients:
% IU(x) = c_0 +c_1(x-x_{i-1})^1 + c_2(x-x_{i-1})^2 +c_3(x-x_{i-1})^3 for
% x in [x_{i-1},x_i].   The way we will be figuring out this coefficients
% is by prescribing
% IU(x_i-1) = U(x_i-1), IU_(x_i) = U(x_i) , and for the derivatives of the
% interpolant we do IU_x(x_i-1) = (U_{x_i-2}-U_{x_i})/2h and
% IU_x(x_i) = (U_{x_i-1}-U_{x_i+1})/2h
% The coefficients of the interpolant can be calculated explicitly from the
% values of u and the difference quotient at the end points of the interval
% (see suli introdcution to numerical analysis page
%%%%%%%%%%%%%%%
% uold_x = (circshift(uold,-1) - circshift(uold,1))/(2*h);
% u_x = (circshift(u,-1) - circshift(u,1))/(2*h);

% u = spline(x,u);
% [breaks,coefs,l,k,d] = unmkpp(u);
% du = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%
% uold = spline(x,uold);
% [breaks,coefs,l,k,d] = unmkpp(uold);
% duold = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%xx=0:h/100:1;
%yy = ppval(u, xx);
%plot(xx,yy)
% IUt_arr = (u - uold)/dt; %defined as a constant for any t and piece-wise quadratic in x
% IU = (evalt - tj)*u/dt + (dt+tj-evalt)*uold/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt

for j = 1 : length(x)-1
    % the values of the coefficients will be calculated on sub-interval by
    % sub-interval basis and are taken from suli page 301
    c_0_old = uold(j);
    c_1_old = uold_x(j);
    c_2_old = 3/(h^2) * (uold(j+1)-uold(j)) - (1/h)*(uold_x(j+1) + 2 * uold_x(j));
    c_3_old = 1/(h^2) * (uold_x(j+1) + uold_x(j)) - 2/(h^3) * (uold(j+1) - uold(j));
    
    c_0_new = u(j);
    c_1_new = u_x(j);
    c_2_new = 3/(h^2) * (u(j+1)-u(j)) - (1/h)*(u_x(j+1) + 2 * u_x(j));
    c_3_new = 1/(h^2) * (u_x(j+1) + u_x(j)) - 2/(h^3) * (u(j+1) - u(j));
    
    xl= zeros(size(nq));
    for iq = 1:nq
        xl(iq) = 0.5*h*xq(iq) + x(j) + h/2;
        diff_x = xl(iq)-x(j);
        IU = (evalt - tj)/dt * (c_0_new + c_1_new * diff_x + c_2_new * diff_x^2 + c_3_new *diff_x^3) + (dt+tj-evalt)/dt * (c_0_old + c_1_old * diff_x + c_2_old * diff_x^2 + c_3_old*diff_x^3); %Lagrange interpolant of IU over the time interval evaluated at evalt
        IUt = 1/dt * (c_0_new + c_1_new * diff_x + c_2_new * diff_x^2 + c_3_new *diff_x^3) - 1/dt * (c_0_old + c_1_old * diff_x + c_2_old * diff_x^2 + c_3_old*diff_x^3);
        %         IUt = (ppval(u,xl(iq)) - ppval(uold,xl(iq)))/dt; %defined as a constant for any t
        %         IUx = (evalt - tj)/dt * (ppval(duold,xl(iq))) + (dt+tj-evalt)/dt * (ppval(du,xl(iq)))/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
        IUx = (evalt - tj)/dt * (c_1_new + 2 * c_2_new * diff_x + 3 * c_3_new*diff_x^2) + (dt+tj-evalt)/dt * (c_1_old + 2 * c_2_old * diff_x + 3* c_3_old*diff_x^2); %Lagrange interpolant of IU over the time interval evaluated at evalt
        
        %IUt(iq) = (xl(iq) - x(j))*IUt(j+1)/h + (x(j+1)-xl(iq))*IUt(j)/h;
        
        L2Rt = L2Rt + wq(iq)*h*(IUt + IUx)^2;
    end
end

end


function [L2Rt,c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs3_vector(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
% h = x(2) - x(1);
%%%%%%%%%%%%%%%
% Each closed [x_{i-1},x_i] sub_interval has a piece-wise cubic Hermite
% spline interpolant in it.  This will be defined by four coeffiecients:
% IU(x) = c_0 +c_1(x-x_{i-1})^1 + c_2(x-x_{i-1})^2 +c_3(x-x_{i-1})^3 for
% x in [x_{i-1},x_i].   The way we will be figuring out this coefficients
% is by prescribing
% IU(x_i-1) = U(x_i-1), IU_(x_i) = U(x_i) , and for the derivatives of the
% interpolant we do IU_x(x_i-1) = (U_{x_i-2}-U_{x_i})/2h and
% IU_x(x_i) = (U_{x_i-1}-U_{x_i+1})/2h
% The coefficients of the interpolant can be calculated explicitly from the
% values of u and the difference quotient at the end points of the interval
% (see suli introdcution to numerical analysis page
%%%%%%%%%%%%%%%
% uold_x = (circshift(uold,-1) - circshift(uold,1))./(dist_x_pl+dist_x_min);
% u_x = (circshift(u,-1) - circshift(u,1))./(dist_x_pl+dist_x_min);

% u = spline(x,u);
% [breaks,coefs,l,k,d] = unmkpp(u);
% du = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%
% uold = spline(x,uold);
% [breaks,coefs,l,k,d] = unmkpp(uold);
% duold = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
%xx=0:h/100:1;
%yy = ppval(u, xx);
%plot(xx,yy)
% IUt_arr = (u - uold)/dt; %defined as a constant for any t and piece-wise quadratic in x
% IU = (evalt - tj)*u/dt + (dt+tj-evalt)*uold/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt

% for j = 1 : length(x)-1
% the values of the coefficients will be calculated on sub-interval by
% sub-interval basis and are taken from suli page 301

c_0_old = uold;
c_1_old = uold_x;
c_2_old = 3./(dist_x_pl.^2) .* (circshift(uold,-1) - uold) - (1./dist_x_pl).*(circshift(uold_x,-1) + 2 * uold_x);
c_3_old = 1./(dist_x_pl.^2) .* (circshift(uold_x,-1) + uold_x) - (2./dist_x_pl.^3) .* (circshift(uold,-1) - uold);

c_0_new = u;
c_1_new = u_x;
c_2_new = 3./(dist_x_pl.^2) .* (circshift(u,-1)-u) - (1./dist_x_pl).*(circshift(u_x,-1) + 2 * u_x);
c_3_new = 1./(dist_x_pl.^2) .* (circshift(u_x,-1) + u_x) - (2./dist_x_pl.^3) .* (circshift(u,-1) - u);

c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old; c_3_old];
c_0_coeff_arr_new= [c_0_new; c_1_new; c_2_new; c_3_new];

%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5* dist_x_pl;
    diff_x = xiq-x;
    IU = (evalt - tj)/dt * (c_0_new + c_1_new .* diff_x + c_2_new .* diff_x.^2 + c_3_new .*diff_x.^3) + (dt+tj-evalt)/dt * (c_0_old + c_1_old .* diff_x + c_2_old .* diff_x.^2 + c_3_old.*diff_x.^3); %Lagrange interpolant of IU over the time interval evaluated at evalt
    
    IUt = 1/dt * (c_0_new + c_1_new .* diff_x + c_2_new .* diff_x.^2 + c_3_new .*diff_x.^3) - 1/dt * (c_0_old + c_1_old .* diff_x + c_2_old .* diff_x.^2 + c_3_old.*diff_x.^3);
    %         IUt = (ppval(u,xl(iq)) - ppval(uold,xl(iq)))/dt; %defined as a constant for any t
    %         IUx = (evalt - tj)/dt * (ppval(duold,xl(iq))) + (dt+tj-evalt)/dt * (ppval(du,xl(iq)))/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
    IUx = (evalt - tj)/dt * (c_1_new + 2 * c_2_new .* diff_x + 3 * c_3_new .* diff_x.^2) + (dt+tj-evalt)/dt * (c_1_old + 2 * c_2_old .* diff_x + 3 * c_3_old .* diff_x.^2); %Lagrange interpolant of IU over the time interval evaluated at evalt
    
    %IUt(iq) = (xl(iq) - x(j))*IUt(j+1)/h + (x(j+1)-xl(iq))*IUt(j)/h;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUt + IUx).^2);
end
% end

end

function [L2Rt,c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs_vector_temp_1_spatiotemp_3(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
%
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients
c_0_t = uold; % this is u at t_n
c_1_t = (1/dt)*(u - uold); % - f_h(U^n)


diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t *diff_t ;
Ut_t = c_1_t ;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
h=dist_x_pl(1);
if spatial_disc == 'CS'
    Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
    Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
elseif spatial_disc == 'BS'
    Ut_x = 1./(dist_x_min).*(-circshift(Ut,1) + circshift(Ut,0));
    Ut_xt = 1./(dist_x_min).*(-circshift(Ut_t,1) + circshift(Ut_t,0));
elseif spatial_disc == '2S'
    Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
    Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
elseif spatial_disc == '3S'
    Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2)); 
    Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
elseif spatial_disc == '4CS'
    Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));    
    Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
else
end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
c_0_ts = Ut; % this is u at t_n
c_1_ts = Ut_x; % - f_h(U^n)
c_2_ts = 3./(dist_x_pl.^2) .* (circshift(Ut,-1) - Ut) - (1./dist_x_pl) .* (circshift(Ut_x,-1) + 2 * Ut_x);
c_3_ts = 1./(dist_x_pl.^2) .* (circshift(Ut_x,-1) + Ut_x) - 2./(dist_x_pl.^3) .* (circshift(Ut,-1) - Ut);

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = Ut_xt; % - f_h(U^n)
c_2_ts_t = 3./(dist_x_pl.^2) .* (circshift(Ut_t,-1) - Ut_t) - (1./dist_x_pl) .* (circshift(Ut_xt,-1) + 2 * Ut_xt);
c_3_ts_t = 1./(dist_x_pl.^2) .* (circshift(Ut_xt,-1) + Ut_xt) - 2./(dist_x_pl.^3) .* (circshift(Ut_t,-1) - Ut_t);



c_0_old = uold;
c_1_old = uold_x;
c_2_old = 3./(dist_x_pl.^2) .* (circshift(uold,-1) - uold) - (1./dist_x_pl).*(circshift(uold_x,-1) + 2 * uold_x);
c_3_old = 1./(dist_x_pl.^2) .* (circshift(uold_x,-1) + uold_x) - 2./(dist_x_pl.^3) .* (circshift(uold,-1) - uold);


c_0_new = u;
c_1_new = u_x;
c_2_new = 3./(dist_x_pl.^2) .* (circshift(u,-1)-u) - (1./dist_x_pl).*(circshift(u_x,-1) + 2 * u_x);
c_3_new = 1./(dist_x_pl.^2) .* (circshift(u_x,-1) + u_x) - (2./dist_x_pl.^3) .* (circshift(u,-1) - u);

c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old; c_3_old];
c_0_coeff_arr_new=  [c_0_new; c_1_new; c_2_new; c_3_new];


%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU = c_0_ts + c_1_ts .* diff_x + c_2_ts .* diff_x .^2 + c_3_ts .* diff_x.^3;
    IUx = c_1_ts + 2* c_2_ts .*diff_x + 3* c_3_ts .*diff_x.^2;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x +  c_2_ts_t.*diff_x.^2 + c_3_ts_t.*diff_x.^3;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUt + IUx).^2);
end
% end

end

function [L2Rt,c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs3_vector_temp_3_spatiotemp_3(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc,dist_alpha_coeffs, dist_alpha_arr)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
c=1;
dx =x(2)- x(1);
%
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients
c_0_t = uold; % this is u at t_n
c_1_t = -uold_x; % - f_h(U^n)
c_2_t = 3/(dt^2) * (u - uold) + (1/dt) * (u_x + 2 * uold_x);
c_3_t = -1/(dt^2) .* (u_x + uold_x) - 2/(dt^3) .* (u - uold);

diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t *diff_t +c_2_t *diff_t^2 + c_3_t*diff_t^3;
Ut_t = c_1_t + 2*c_2_t*diff_t +3*c_3_t*diff_t.^2;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
h=dist_x_pl(1);
if spatial_disc == 'CS'
    Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
    Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
elseif spatial_disc == 'BS'
    Ut_x = 1./(dist_x_min).*(-circshift(Ut,1) + circshift(Ut,0));
    Ut_xt = 1./(dist_x_min).*(-circshift(Ut_t,1) + circshift(Ut_t,0));
elseif spatial_disc == '2S'
    Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
    Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
elseif spatial_disc == '3S'
%     Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2)); 
%     Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
    
    Ut_x = resWENO5(Ut,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut);
    Ut_xt = resWENO5(Ut_t,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,Ut_t);
elseif spatial_disc == '4CS'
    Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));    
    Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
else
end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
c_0_ts = Ut; % this is u at t_n
c_1_ts = Ut_x; % - f_h(U^n)
c_2_ts = 3./(dist_x_pl.^2) .* (circshift(Ut,-1) - Ut) - (1./dist_x_pl) .* (circshift(Ut_x,-1) + 2 * Ut_x);
c_3_ts = 1./(dist_x_pl.^2) .* (circshift(Ut_x,-1) + Ut_x) - 2./(dist_x_pl.^3) .* (circshift(Ut,-1) - Ut);

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = Ut_xt; % - f_h(U^n)
c_2_ts_t = 3./(dist_x_pl.^2) .* (circshift(Ut_t,-1) - Ut_t) - (1./dist_x_pl) .* (circshift(Ut_xt,-1) + 2 * Ut_xt);
c_3_ts_t = 1./(dist_x_pl.^2) .* (circshift(Ut_xt,-1) + Ut_xt) - 2./(dist_x_pl.^3) .* (circshift(Ut_t,-1) - Ut_t);



c_0_old = uold;
c_1_old = uold_x;
c_2_old = 3./(dist_x_pl.^2) .* (circshift(uold,-1) - uold) - (1./dist_x_pl).*(circshift(uold_x,-1) + 2 * uold_x);
c_3_old = 1./(dist_x_pl.^2) .* (circshift(uold_x,-1) + uold_x) - 2./(dist_x_pl.^3) .* (circshift(uold,-1) - uold);


c_0_new = u;
c_1_new = u_x;
c_2_new = 3./(dist_x_pl.^2) .* (circshift(u,-1)-u) - (1./dist_x_pl).*(circshift(u_x,-1) + 2 * u_x);
c_3_new = 1./(dist_x_pl.^2) .* (circshift(u_x,-1) + u_x) - (2./dist_x_pl.^3) .* (circshift(u,-1) - u);

c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old; c_3_old];
c_0_coeff_arr_new=  [c_0_new; c_1_new; c_2_new; c_3_new];



%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU = c_0_ts + c_1_ts .* diff_x + c_2_ts .* diff_x .^2 + c_3_ts .* diff_x.^3;
    IUx = c_1_ts + 2* c_2_ts .*diff_x + 3* c_3_ts .*diff_x.^2;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x +  c_2_ts_t.*diff_x.^2 + c_3_ts_t.*diff_x.^3;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUt + IUx).^2);
end
% end

end

function L2Rt = compute_R(x,uold,u,evalt,tj,dt) % RL2iq = compute_R(x,uold,u,tq(iq),t,dt); %compute R at gauss point
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
h = x(2) - x(1);

IUt = (u - uold)/dt; %defined as a constant for any t
IU = (evalt - tj)*u/dt + (dt+tj-evalt)*uold/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
for j = 1 : length(x)-1
    for iq = 1:nq
        xiq = 0.5*h*xq(iq) + x(j) + h/2;
        vl = IU(j);
        vr = IU(j+1);
        IUx = (vr - vl)/h;
        IUtx(iq) = (xiq - x(j))*IUt(j+1)/h + (x(j+1)-xiq)*IUt(j)/h;
        
        L2Rt = L2Rt + wq(iq)*h*(IUtx(iq) + IUx)^2;
    end
end

end

function L2Rt = compute_R_vector(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt) % RL2iq = compute_R(x,uold,u,tq(iq),t,dt); %compute R at gauss point
global nq xq wq
% Vectorial version of comput_R for linear lagrange interpolant
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
% h = x(2) - x(1);

IUt = (u - uold)/dt; %defined as a constant for any t
IU = (evalt - tj)*u/dt + (dt+tj-evalt)*uold/dt; %Lagrange interpolant of IU over the time interval evaluated at evalt
% for j = 1 : length(x)-1
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + dist_x_pl/2;
    vl = IU;
    vr = circshift(IU,-1);
    IUx = (vr - vl)./dist_x_pl;
    
    IUtx = (xiq - x) .* circshift(IUt,-1)./dist_x_pl + (x + dist_x_pl - xiq).*IUt./dist_x_pl;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUtx + IUx).^2);
end
% end

end

function [L2Rt,c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs_vector_temp_1_spatiotemp_1(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
%
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients
c_0_t = uold; % this is u at t_n
c_1_t = 1/dt(u-uold); % - f_h(U^n)


diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t *diff_t;
Ut_t = c_1_t ;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
h=dist_x_pl(1);
if spatial_disc == 'CS'
    Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
    Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
elseif spatial_disc == 'BS'
    Ut_x = 1./(dist_x_min).*(-circshift(Ut,1) + circshift(Ut,0));
    Ut_xt = 1./(dist_x_min).*(-circshift(Ut_t,1) + circshift(Ut_t,0));
elseif spatial_disc == '2S'
    Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
    Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
elseif spatial_disc == '3S'
    Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2)); 
    Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
elseif spatial_disc == '4CS'
    Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));    
    Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
else
end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
c_0_ts = Ut; % this is u at t_n
c_1_ts = 1./dist_x_pl.*(circshift(Ut,-1) - Ut); % - f_h(U^n)

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = 1./(dist_x_pl) .* (circshift(Ut_t,-1) - Ut_t) ; % - f_h(U^n)




c_0_old = uold;
c_1_old = (1./dist_x_pl).*(circshift(uold,-1)-uold);

c_0_new = u;
c_1_new = (1./dist_x_pl).*(circshift(u,-1) - u);
c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old; c_3_old];
c_0_coeff_arr_new=  [c_0_new; c_1_new];


%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU = c_0_ts + c_1_ts .* diff_x;
    IUx = c_1_ts ;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUt + IUx).^2);
end
% end

end

function [L2Rt,c_0_coeff_arr_new, c_0_coeff_arr_old] = compute_Rs_vector_temp_1_spatiotemp_2(x,dist_x_pl,dist_x_min,uold,u,evalt,tj,dt,uold_x,u_x,spatial_disc)
global nq xq wq
%uold defined at tn
%u defined at t
% R = IU_t + IU_x
% Take IU as a cubic spline
% for any t it can be represented as a pw linear function in space on each
% spatial element
L2Rt = 0;
%
% uold is u_{t_n}
% u is
% temporal reconstruction coefficients
c_0_t = uold; % this is u at t_n
c_1_t = (1/dt)*(u - uold); % - f_h(U^n)


diff_t=evalt-tj;
% Now we calculate the value of the spatial discretisation using the values
% of the spatial reconstruction at time t in the interval [t_n, t_{n+1}]
Ut = c_0_t + c_1_t *diff_t ;
Ut_t = c_1_t ;

% We calculate the spatial derivative discretisation using the spatial
% reconstruction at time t
h=dist_x_pl(1);
if spatial_disc == 'CS'
    Ut_x = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut,1) + circshift(Ut,-1));
    Ut_xt = 1./((dist_x_pl+dist_x_min)).*(-circshift(Ut_t,1) + circshift(Ut_t,-1));
elseif spatial_disc == 'BS'
    Ut_x = 1./(dist_x_min).*(-circshift(Ut,1) + circshift(Ut,0));
    Ut_xt = 1./(dist_x_min).*(-circshift(Ut_t,1) + circshift(Ut_t,0));
elseif spatial_disc == '2S'
    Ut_x= 1./(2*h)*( 3*circshift(Ut,0)  - 4*circshift(Ut,1) + circshift(Ut,2) );
    Ut_xt= 1./(2*h)*( 3*circshift(Ut_t,0)  - 4*circshift(Ut_t,1) + circshift(Ut_t,2) );
elseif spatial_disc == '3S'
    Ut_x = 1./(6*h)*( 2*circshift(Ut,-1)  + 3*circshift(Ut,0) - 6* circshift(Ut,1) +circshift(Ut,2)); 
    Ut_xt = 1./(6*h)*( 2*circshift(Ut_t,-1)  + 3*circshift(Ut_t,0) - 6* circshift(Ut_t,1) +circshift(Ut_t,2));
elseif spatial_disc == '4CS'
    Ut_x = -1./(12*h).*(-circshift(Ut,2)+8*circshift(Ut,1) - 8*circshift(Ut,-1) + circshift(Ut,-2));    
    Ut_xt = -1./(12*h).*(-circshift(Ut_t,2)+8*circshift(Ut_t,1) - 8*circshift(Ut_t,-1) + circshift(Ut_t,-2));
else
end

% We calculate the coefficients for the spatio-temporal discretisation and
% also the ones for the temporal derivative
c_0_ts = Ut; % this is u at t_n
c_1_ts = Ut_x; % - f_h(U^n)
c_2_ts = 1./(dist_x_pl.^2) .* (circshift(Ut,-1) - (1./dist_x_pl) .* (Ut + dist_x_pl.*Ut_x));

c_0_ts_t = Ut_t; % this is u at t_n
c_1_ts_t = Ut_xt; % - f_h(U^n)
c_2_ts_t =  1./(dist_x_pl.^2) .* (circshift(Ut_t,-1) - (1./dist_x_pl) .* (Ut_t + dist_x_pl.*Ut_xt));


c_0_old = uold;
c_1_old = uold_x;
c_2_old = 1./(dist_x_pl.^2).*(circshift(uold,-1) - (uold + uold_x .* dist_x_pl));

c_0_new = u;
c_1_new = u_x;
c_2_new = (1./(dist_x_pl.^2)).*(circshift(u,-1) - (u + u_x .* dist_x_pl));

c_0_coeff_arr_old= [c_0_old; c_1_old; c_2_old];
c_0_coeff_arr_new=  [c_0_new; c_1_new; c_2_new];


%     xl= zeros(size(nq));
for iq = 1:nq
    xiq = 0.5 * dist_x_pl * xq(iq) + x + dist_x_pl/2;
    diff_x = xiq-x;
    
    IU = c_0_ts + c_1_ts .* diff_x +c_2_ts.*diff_x.^2;
    IUx = c_1_ts + 2 * c_2_ts.*diff_x;
    IUt = c_0_ts_t + c_1_ts_t.*diff_x +  c_2_ts_t.*diff_x.^2;
    
    L2Rt = L2Rt + sum(wq(iq)*dist_x_pl.*(IUt + IUx).^2);
end
% end

end


%At a specific time level compute L2 spatial error use 2 point gauss quadrature elementwise
function int = space_int(x,u,evalt,ex) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq

h = x(2)-x(1);
int = 0;
for j = 1 : length(x)-1
    for iq = 1:nq
        xiq = 0.5*h*xq(iq) + x(j) + h/2;
        v(iq) = (xiq - x(j))*u(j+1)/h + (x(j+1)-xiq)*u(j)/h;
        
        int = int + wq(iq)*h*(v(iq) - ex(xiq,evalt))^2.;
    end
end

end


function int = space_int_vector(x,dist_x_pl,dist_x_min,u,evalt,ex,c_0_coeff_arr) % sqrt(space_int(x,u,T,ex)); %compute final time L2 error
global nq xq wq
% c_coefficients = [c_0; c_1; c_2; c_3];
% h = x(2)-x(1);
int = 0;
% for j = 1 : length(x)-1
l_coeffs = length(c_0_coeff_arr);
for iq = 1:nq
    xiq = 0.5*dist_x_pl*xq(iq) + x + .5*dist_x_pl;
    diff_x = xiq-x;
    IU = sum(c_0_coeff_arr.*[ones(1,l_coeffs); diff_x; diff_x.^2; diff_x.^3],1);
%     IU = sum(c_0_coeff_arr.*[ones(1,l_coeffs); diff_x; diff_x.^2],1);

%     v = (xiq - x).*circshift(u,-1)./dist_x_pl + (x + dist_x_pl - xiq).*u./dist_x_pl;
    
    int = int + sum(wq(iq).* dist_x_pl .* (IU - ex(xiq,evalt)).^2);
end
% end

end

function uold_x_nu  = fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu)
                        uold_x_nu = ( dist_alpha_coeffs(1,:).*circshift(uold_nu,-1)./dist_alpha_arr(1,:)...
                            + dist_alpha_coeffs(2,:).*circshift(uold_nu,0)./dist_alpha_arr(2,:)...
                            + dist_alpha_coeffs(3,:).* circshift(uold_nu,1)./dist_alpha_arr(3,:)...
                            + dist_alpha_coeffs(4,:).*circshift(uold_nu,2)./dist_alpha_arr(4,:));
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