clc;
clear;
close all;

c={};
% scheme_arr = {scheme_name_str};
exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'qdr'};



c=struct2cell(load(['SSP3WENO3_cell_arr_file_sinIC.mat']));
%     c{i}=struct2cell('LeapFrog_cell_arr_file_sin_IC_uniform.mat');
% colour_arr =['y';'c';'g';"mo";"ro";"bo";"ko"];

colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);
which_scheme = 1;
xlim_plot = 40;
ylim_plot1 = 1e-12;
ylim_plot2 = 5e-6;
xlim_p = xlim_plot;
ylim_p = [ylim_plot1, ylim_plot2];
n_schemes= length(c);
n_plots = 4;
f_size= 14;
N_subplots = 5;
%% plots
pow=1;
figgy= figure();

l_refs = length(c{1});
for i_ref = 1: l_refs
    
    
    time_arr  =  c{1}{i_ref}(1,:);
    bound_arr =  c{1}{i_ref}(2,:);
    error_arr =  c{1}{i_ref}(3,:);
    EI_arr    =  c{1}{i_ref}(4,:);
   
    % Evolution of the error
    subplot(1,N_subplots,1)
    semilogy(time_arr,error_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
    ylabel(['$Error$'],'Interpreter','latex','FontSize',f_size)
    xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
    grid on;
    xlim([0, xlim_p])
    ylim(ylim_p)
    hold on;
    
    % Evolution of the indicator
    subplot(1,N_subplots,2)
    semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
    ylabel(['$Estimator$'],'Interpreter','latex','FontSize',f_size)
    xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
    grid on;
    xlim([0, xlim_p])
    ylim(ylim_p)
    hold on;
    
    
    if i_ref<l_refs
        time_arr_2  = c{1}{i_ref+1}(1,:);
        bound_arr_2 = c{1}{i_ref+1}(2,:);
        error_arr_2 = c{1}{i_ref+1}(3,:);
        
        time_arr_1  = c{1}{i_ref}(1,:);
        bound_arr_1 = c{1}{i_ref}(2,:);
        error_arr_1 = c{1}{i_ref}(3,:);
        
        if length(error_arr_2)>length(error_arr_1)
            eoc_time_arr = time_arr_2(1:2^(pow):end)./time_arr_1(1:end);
            EOC_error = log((error_arr_2(1:2^(pow):end))./(error_arr_1(1:end)))/log(0.5);
            EOC_bound = log((bound_arr_2(1:2^pow:end))./(bound_arr_1(1:end)))/log(0.5);
            

        else
            EOC_error = log((error_arr_2)./(error_arr_1))/log(0.5);
            EOC_bound = log((bound_arr_2)./(bound_arr_1))/log(0.5);
        end
        
       
            subplot(1,N_subplots,3)
            plot(time_arr_1,EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
%             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            ylabel(['$EOC\left( Error \right)$'],'Interpreter','latex')
            xlabel(['$t_i$'],'Interpreter','latex')
            grid on;
         %   title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            xlim([0,xlim_p])
            ylim([0,inf])
            hold on;
            
            subplot(1,N_subplots,4)
            
            plot(time_arr_1,EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            ylabel(['$EOC\left( Estimator \right)$'],'Interpreter','latex')
            xlabel('$t_i$','Interpreter','latex')
          %  title(['$\mathcal{P}^',num2str(3),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([0,xlim_p])
            ylim([0,inf])
            hold on;
        
        
    end
    
      
    subplot(1,N_subplots,5)
    plot(time_arr,EI_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
    
    ylabel(['$EI$'],'Interpreter','latex')
    
    xlabel('$t_i$','Interpreter','latex')
    xlim([0,xlim_p])
    ylim([0,inf])
    grid on;
    hold on;
    
    
end
hold off;
set(gcf, 'Position',  [100, 100, 1050, 200])
%          print(gcf,['/Users/gs1511/Desktop/GSialounas_repo_temporal_rec/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_harmonic_u',num2str(round(10*u_init)),'_v',num2str(round(10*v_init)),'_paperrec_poly_tristan.png'],'-dpng','-r600')

saveas(figgy,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_SSP3WENO5_sinIC_plots_1x5_linadvect'],'png');%char(Scheme_string)


hold off;