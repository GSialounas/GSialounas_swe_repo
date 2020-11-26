function var_out = fn_plot_WENO_burger(init_conds)

c={};
scheme_arr = {'SSP3WENO3'};
exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'qdr'};



for i = 1 : length(scheme_arr)
    c{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_sinIC_burger.mat']));
end


colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);
which_scheme = 1;

n_schemes= length(c);
n_plots = 4;
f_size= 14;
%% plots

i_int=1;
for i_scheme = 1:n_schemes
    figgy= figure();
    
    l_int = length(interpolant_arr); % we will only present dt~dx^2
    %     for  i_int = 1:l_int
    l_refs = length(c{i_scheme}{1});
    for i_ref = 1: l_refs
        
        time_arr  =  c{i_scheme}{1}{i_ref}(1,:);
        bound_arr =  c{i_scheme}{1}{i_ref}(2,:);
        error_arr =  c{i_scheme}{1}{i_ref}(3,:);
        EI_arr    =  c{i_scheme}{1}{i_ref}(4,:);

        
        
        % evolution of the error
        subplot(1,5,1)
%         title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')       
        semilogy(time_arr,error_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
        xlabel('$t_i$','Interpreter','Latex','Fontsize',f_size)
%         ylabel(['$ \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')   
        ylabel(['$Error$'],'Interpreter','latex','Fontsize',f_size)
        grid on;
        xlim([0,0.5])
        ylim([1e-10,1e-2])
        hold on;
        
        % Evolution of the indicator
        subplot(1,5,2)
        semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)       
%         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')        
%         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        xlabel('$t_i$','Interpreter','latex','Fontsize',f_size)
        ylabel('$Estimator$','Interpreter','latex','Fontsize',f_size)

        grid on;
        xlim([0,0.5])
        ylim([1e-10,1e-2])
        hold on;
        
        
        if i_ref<l_refs
            time_arr_2  = c{i_scheme}{1}{i_ref+1}(1,:);
            bound_arr_2 = c{i_scheme}{1}{i_ref+1}(2,:);
            error_arr_2 = c{i_scheme}{1}{i_ref+1}(3,:);
            
            time_arr_1  = c{i_scheme}{1}{i_ref}(1,:);
            bound_arr_1 = c{i_scheme}{1}{i_ref}(2,:);
            error_arr_1 = c{i_scheme}{1}{i_ref}(3,:);
%             length(time_arr_2(1:2:end))
%             length(time_arr_1)
            
            if length(error_arr_2)>length(error_arr_1)
                EOC_error = log((error_arr_2(1:(2^(i_ref))^(1):end))./(error_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                EOC_bound = log((bound_arr_2(1:(2^(i_ref))^1:end))./(bound_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
            else
                EOC_error = log((error_arr_2)./(error_arr_1))/log(0.5);
                EOC_bound = log((bound_arr_2)./(bound_arr_1))/log(0.5);
            end

            subplot(1,5,3)
            if length(time_arr_1)>length(EOC_error)
            plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            else
                plot(time_arr_1,EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
            end

%             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left(Error\right)$'],'Interpreter','latex','Fontsize',f_size)

            xlabel(['$t_i$'],'Interpreter','latex','Fontsize',f_size)
            grid on;
%             title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
            xlim([0,0.5])
            ylim([0,6])
            hold on;
            
            subplot(1,5,4)
              if length(time_arr_1)>length(EOC_bound)
            plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
              else
                  plot(time_arr_1,EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
              end

%             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
            ylabel(['$EOC\left( Estimator\right)$'],'Interpreter','latex','Fontsize',f_size)

            xlabel('$t_i$','Interpreter','latex','Fontsize',f_size)
%             title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            grid on;
            xlim([0,0.5])
            ylim([0,6])
            hold on;
            %                 ylim([0,3])
            
            
            
            
        end
        
        
        
        
        subplot(1,5,5)
        plot(time_arr,EI_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
        
        ylabel(['$EI\left(t_i\right)$'],'Interpreter','latex','Fontsize',f_size)
        
%         title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
        xlabel('$t_i$','Interpreter','latex','Fontsize',f_size)
%         xlim([0,0.5])
        ylim([0,10])
%         ylim([0,10])

        grid on;
        hold on;
    end
    %     end
    
    %     sgtitle([scheme_arr{i_scheme},' bound'],'Color','red')
    %     sgt.FontSize = 20;
    hold off;
    set(gcf, 'Position',  [100, 100, 1050, 200])
saveas(figgy,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_SSP3WENO3_',init_conds,'_plots_1x5_burgers'],'png');%char(Scheme_string)


hold off;
var_out =1;
end
