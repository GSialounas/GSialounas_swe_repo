clc;
clear;
close all;


c={};
scheme_arr = {'SSP3WENO'};
exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'spl'};



for i = 1 : length(scheme_arr)
    c{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_sin_IC.mat']));
end

colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);
which_scheme = 1;

n_schemes= length(c);
n_plots = 4;
xlim_p = 2;

%% plots

i_int=1;
for i_scheme = 1:n_schemes
    figgy= figure();
    
    l_int = length(interpolant_arr); % we will only present dt~dx^2
    %     for  i_int = 1:l_int
    l_refs = length(c{i_scheme}{1}{1}{i_int});
    for i_ref = 1: l_refs
   
                
        time_arr = c{i_scheme}{1}{1}{i_int}{i_ref}(1,:);
        bound_arr = c{i_scheme}{1}{1}{i_int}{i_ref}(2,:);
        error_arr = c{i_scheme}{1}{1}{i_int}{i_ref}(3,:);
        EI_arr = c{i_scheme}{1}{1}{i_int}{i_ref}(4,:);
              subplot(1,5,1)
        %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
        
        semilogy(time_arr,error_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t_i$','Interpreter','Latex')
        ylabel(['$ Error $'],'Interpreter','latex')
        
        grid on;
        xlim([0,xlim_p])
        ylim([1e-11,1e-4])
        hold on;
        
        
        subplot(1,5,2)
        semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
        
        ylabel(['$ Estimator $'],'Interpreter','latex')
        
        %title(['$\mathcal{P}^',num2str(3),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        xlabel('$t_i$','Interpreter','latex')
        grid on;
        xlim([0,xlim_p])
        ylim([1e-11,1e-4])
        
        hold on;
       
        
        if i_ref<l_refs
            time_arr_2 = c{i_scheme}{1}{1}{i_int}{i_ref+1}(1,:);
            bound_arr_2 = c{i_scheme}{1}{1}{i_int}{i_ref+1}(2,:);
            error_arr_2 = c{i_scheme}{1}{1}{i_int}{i_ref+1}(3,:);
            
            time_arr_1 = c{i_scheme}{1}{1}{i_int}{i_ref}(1,:);
            bound_arr_1 = c{i_scheme}{1}{1}{i_int}{i_ref}(2,:);
            error_arr_1 = c{i_scheme}{1}{1}{i_int}{i_ref}(3,:);
            
            EOC_error = log((error_arr_2(1:(2^(i_ref))^(1):end))./(error_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            EOC_bound = log((bound_arr_2(1:(2^(i_ref))^1:end))./(bound_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            
            subplot(1,5,3)
            plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            ylabel(['$EOC\left( Error \right)$'],'Interpreter','latex')
            xlabel(['$t_i$'],'Interpreter','latex')
            grid on;
            %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            xlim([0,xlim_p])
            ylim([0,6])
            hold on;
            
            subplot(1,5,4)
            
            plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            ylabel(['$EOC\left( Estimator\right)$'],'Interpreter','latex')
            xlabel('$t_i$','Interpreter','latex')
            %title(['$\mathcal{P}^',num2str(3),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([0,xlim_p])
            ylim([0,4])
            hold on;
            %                 ylim([0,3])
            
            
            
            
        end
       
        
         
        
        subplot(1,5,5)
        plot(time_arr,EI_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
        
        ylabel(['$EI$'],'Interpreter','latex')
        
        %title(['$\mathcal{P}^',num2str(3),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        xlabel('$t_i$','Interpreter','latex')
        xlim([0,xlim_p])
        ylim([0,inf])
        grid on;
        hold on;
    end
    %     end
    
    %     sg%title([scheme_arr{i_scheme},' bound'],'Color','red')
    %     sgt.FontSize = 20;
    hold off;
    set(gcf, 'Position',  [100, 100, 1050, 200])
    
    saveas(figgy,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC'],'png');%char(Scheme_string)
    
end

hold off;


