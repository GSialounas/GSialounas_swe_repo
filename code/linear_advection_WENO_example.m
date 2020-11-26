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
    
    
    figgy = figure;
    for i_exponent = 1:length(exponent_arr)
        
       
        for i_interpolant = 1:length(intrpl_arr)
            
           
            
            scheme = scheme_arr{i_scheme}; % 'FTBS' order 1 time, order 1 space upwinding
            
            T = 2; % final time
            maxit =4; % max refinement iterations
            
            cntr = .25;
            Lx=1;
            pow_exp=100;
            ratio_cTf=2;
            part_fine=.5;
            ex = @(x,t) sin(2*pi*(x-t));
   
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
                    
                  
                    spat_disc = "WENO";
                    c=1;
                    uold_x_nu  = resWENO5(uold_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,uold_nu);
             
                    % Three stage
                    ustage_0 = uold_nu;
                    
                    ustage_1 = uold_nu - dt*uold_x_nu;
                    uold_x_stage_1 = resWENO5(ustage_1,c,dx);% fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_1);
                   
                    ustage_2 = (.75)*uold_nu +.25*ustage_1 - .25 *dt *(uold_x_stage_1);
                    uold_x_stage_2 = resWENO5(ustage_2,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_2);                   
                    
                    ustage_3 = (1/3)*uold_nu +(2/3)*ustage_2 - (2/3)*dt*(uold_x_stage_2);
                    
                    u_nu = ustage_3;
                    u_x_nu = resWENO5(u_nu,c,dx);%fn_B3S(dist_alpha_coeffs, dist_alpha_arr,ustage_3);
                    
                    
                    
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
                    
                     if (showplot && mod(reps,fps)==0)
                        l_coef= length(c_0_coeff_arr_new);
                        IU = sum(c_0_coeff_arr_new.*[ones(1,l_coef);dist_x_pl;dist_x_pl.^2;dist_x_pl.^3],1);
                        plot(x,u_nu,'b',x,IU,'r*') % plot(x,u_nu,'b',x,ex(x,t+dt),'r')
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
        end      
    end
%     set(gcf, 'Position',  [100, 100, 1000, 500]) % position position width height
    
    save([scheme_arr{i_scheme},'_cell_arr_file_sin_IC.mat'],'cell_cell_arr')
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