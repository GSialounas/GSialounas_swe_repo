function var_out = fn_plot_shw_periodic(u_init,v_init, ylim_plot1, ylim_plot2, xlim_plot, EI_plot_y_lim,R_plot_lims,pow,scheme_name_str)


c={};
scheme_arr = {scheme_name_str};
exponent_arr = [1];%[1,2];

% interpolant_arr = {'pwl','qdr','spl'};
interpolant_arr = {'qdr'};



for i = 1 : length(scheme_arr)
    c{i}=struct2cell(load([scheme_arr{i},'_cell_arr_file_shw_periodic.mat']));
    %     c{i}=struct2cell('LeapFrog_cell_arr_file_sin_IC_uniform.mat');
    
end
% colour_arr =['y';'c';'g';"mo";"ro";"bo";"ko"];

colour_arr =['y';'c';'g';'m';'r';'b';'k'];
n_entries = length(c);
which_scheme = 1;

xlim_p = xlim_plot;
ylim_p = [ylim_plot1, ylim_plot2];
n_schemes= length(c);
n_plots = 4;
f_size= 14;
N_subplots=4;
%% plots

i_int=1;
for i_scheme = 1:n_schemes
    %
   
        figgy= figure();
        
        l_refs = length(c{i_scheme}{1});
        for i_ref = 1: l_refs
            
            
            time_arr  =  c{i_scheme}{1}{i_ref}(1,:);
            bound_arr =  c{i_scheme}{1}{i_ref}(2,:);
            error_arr =  c{i_scheme}{1}{i_ref}(3,:);
            EI_arr    =  c{i_scheme}{1}{i_ref}(4,:);
            dofs_arr  =  c{i_scheme}{1}{i_ref}(5,:);
            bound_arr_h =  c{i_scheme}{1}{i_ref}(6,:);
            bound_arr_hv =   c{i_scheme}{1}{i_ref}(7,:);
            
%             error_arr_el =  c{i_scheme}{1}{i_ref}(5,:);
%             error_arr_u =  c{i_scheme}{1}{i_ref}(6,:);
%             error_arr_v =  c{i_scheme}{1}{i_ref}(7,:);
            
            
            % evolution of the error
%             subplot(1,N_subplots,1)
%             %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%             %             error_arr(1)=1e-9;
%             semilogy(time_arr,error_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
%             xlabel('$t_i$','Interpreter','Latex','FontSize',f_size)
%             %         ylabel(['$ \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
%             ylabel(['$Error$'],'Interpreter','latex','FontSize',f_size)
%             
%             grid on;
%             %             xlim([0,1])
%             xlim([0, xlim_p])
% %             ylim(ylim_p)
%             hold on;
            
            % Evolution of the indicator
            subplot(1,N_subplots,1)
            semilogy(time_arr,bound_arr,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
            %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
            ylabel(['$Estimator$'],'Interpreter','latex','FontSize',f_size)
            
            %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
            grid on;
            xlim([0, xlim_p])
%             ylim(ylim_p)
            hold on;
%             
% 
%             subplot(1,N_subplots,3)
%             %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%             %             error_arr(1)=1e-9;
%             semilogy(time_arr,error_arr_el,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
%             xlabel('$t_i$','Interpreter','Latex','FontSize',f_size)
%             %         ylabel(['$ \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
%             ylabel(['$Error\, e_L$'],'Interpreter','latex','FontSize',f_size)
%             
%             grid on;
%             %             xlim([0,1])
%             xlim([0, xlim_p])
%             ylim(ylim_p)
%             hold on;
%             
            if i_ref<l_refs
                time_arr_2  = c{i_scheme}{1}{i_ref+1}(1,:);
                bound_arr_2 = c{i_scheme}{1}{i_ref+1}(2,:);
                error_arr_2 = c{i_scheme}{1}{i_ref+1}(3,:);

                time_arr_1  = c{i_scheme}{1}{i_ref}(1,:);
                bound_arr_1 = c{i_scheme}{1}{i_ref}(2,:);
                error_arr_1 = c{i_scheme}{1}{i_ref}(3,:);
                
                if length(error_arr_2)>length(error_arr_1)
                    eoc_time_arr = time_arr_2(1:2^(pow):end)./time_arr_1(1:end);
                    EOC_error = log((error_arr_2(1:2^(pow):end))./(error_arr_1(1:end)))/log(0.5);
                    EOC_bound = log((bound_arr_2(1:2^pow:end))./(bound_arr_1(1:end)))/log(0.5);
                    
%                     EOC_R1 = log((EOC_R1_arr_2(1:(2^(i_ref))^(1):end))./(EOC_R1_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
%                     EOC_R2 = log((EOC_R2_arr_2(1:(2^(i_ref))^(1):end))./(EOC_R2_arr_1(1:(2^(i_ref-1))^(1):end)))/log(0.5);
                    
                else
                    EOC_error = log((error_arr_2)./(error_arr_1))/log(0.5);
                    EOC_bound = log((bound_arr_2)./(bound_arr_1))/log(0.5);
                end
                
%                 subplot(1,N_subplots,3)
%                 if length(time_arr_1)>length(EOC_error)
%                     plot(time_arr_1,EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
%                 else
%                     plot(time_arr_1,EOC_error,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
%                 end
%                 
%                 %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%                 %             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
%                 ylabel(['$EOC\left(Error\right)$'],'Interpreter','latex','FontSize',f_size)
%                 
%                 xlabel(['$t_i$'],'Interpreter','latex','FontSize',f_size)
%                 grid on;
%                 %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%                 xlim([0, xlim_p])
%                 ylim([0,4])
%                 hold on;
                
                subplot(1,N_subplots,2)
                if length(time_arr_1)>length(EOC_bound)
                    plot(time_arr_1(1:2^(pow):end),EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
                else
                    plot(time_arr_1,EOC_bound,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
                end
                ylabel(['$EOC\left(Estimator\right)$'],'Interpreter','latex','FontSize',f_size)
                
                %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
                xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
                %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
                grid on;
                xlim([0, xlim_p])
                ylim([0,5])
                hold on;
                %                 ylim([0,3])
                
                 
%                 subplot(1,N_subplots,7)
%                 if length(time_arr_1)>length(EOC_R1)
%                     plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_R1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
%                 else
%                     plot(time_arr_1,EOC_R1,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
%                 end
%                 
%                 %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%                 %             ylabel(['$EOC\left( \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
%                 ylabel(['$EOC\left(R1\right)$'],'Interpreter','latex','FontSize',f_size)
%                 
%                 xlabel(['$t_i$'],'Interpreter','latex','FontSize',f_size)
%                 grid on;
%                 %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%                 xlim([0, xlim_p])
%                 ylim([0,4])
%                 hold on;
%                 
%                 subplot(1,N_subplots,8)
%                 if length(time_arr_1)>length(EOC_R2)
%                     plot(time_arr_1(1:(2^(i_ref-1))^(1):end),EOC_R2,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
%                 else
%                     plot(time_arr_1,EOC_R2,colour_arr(end-(l_refs-i_ref)+1),'LineWidth',2)
%                 end
%                 ylabel(['$EOC\left(R2\right)$'],'Interpreter','latex','FontSize',f_size)
%                 
%                 %             ylabel(['$EOC\left( \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} \right)$'],'Interpreter','latex')
%                 xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
%                 %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
%                 grid on;
%                 xlim([0, xlim_p])
%                 ylim([0,4])
%                 hold on;
                %                 ylim([0,3])
                
                
                
                
            end
            
            
%             
%             
%             subplot(1,N_subplots,5)
% %                         plot(time_arr,EI_arr+(l_refs-i_ref)*.05,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
% 
%             plot(time_arr(1:end),EI_arr(1:end),colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
%             
%             ylabel(['$EI\left(t_i\right)$'],'Interpreter','latex','FontSize',f_size)
%             
%             %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
%             xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
%             xlim([0, xlim_p])
%                         ylim([0,EI_plot_y_lim])
%             %                     ylim([0,5])
%             
%             grid on;
%             hold on;
%             

%   subplot(1,N_subplots,7)
% %                         plot(time_arr,EI_arr+(l_refs-i_ref)*.05,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
% 
%             semilogy(time_arr(1:end),error_arr_u(1:end),colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
%             
%             ylabel(['$error u\left(t_i\right)$'],'Interpreter','latex','FontSize',f_size)
%             
%             %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
%             xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
%             xlim([0, xlim_p])
%             ylim(ylim_p)
% 
% %                         ylim([0,max(EI_arr)])
%             %                     ylim([0,5])
%             
%             grid on;
%             hold on;
% 
%   subplot(1,N_subplots,8)
% %                         plot(time_arr,EI_arr+(l_refs-i_ref)*.05,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
% 
%             semilogy(time_arr(1:end),error_arr_v(1:end),colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
%             
%             ylabel(['$error v\left(t_i\right)$'],'Interpreter','latex','FontSize',f_size)
%             
%             %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
%             xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
%             xlim([0, xlim_p])
%             ylim(ylim_p)
% 
% %                         ylim([0,max(EI_arr)])
%             %                     ylim([0,5])
%             
%             grid on;
%             hold on;
%             subplot(1,N_subplots,6)
%             %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%             %             error_arr(1)=1e-9;
%             semilogy(time_arr,error_arr_num,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
%             xlabel('$t_i$','Interpreter','Latex','FontSize',f_size)
%             %         ylabel(['$ \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
%             ylabel(['$Numerical\, Error$'],'Interpreter','latex','FontSize',f_size)
% 
%             grid on;
%             %             xlim([0,1])
%             xlim([0, xlim_p])
%             ylim(ylim_p)
%             hold on;
            
%             subplot(1,N_subplots,7)
%             %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%             %             error_arr(1)=1e-9;
%             semilogy(time_arr,R1_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
%             xlabel('$t_i$','Interpreter','Latex','FontSize',f_size)
%             %         ylabel(['$ \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
%             ylabel(['$R1\,T$'],'Interpreter','latex','FontSize',f_size)
% 
%             grid on;
%             %             xlim([0,1])
%             xlim([0, xlim_p])
%             ylim(R_plot_lims)
%             hold on;
%             
%              subplot(1,N_subplots,8)
%             %title(['$dt=0.1dx^', num2str(exponent_arr(1)),'$'],'Interpreter','latex')
%             %             error_arr(1)=1e-9;
%             semilogy(time_arr,R2_arr,[colour_arr(end-(l_refs-i_ref))],'LineWidth',2)
%             xlabel('$t_i$','Interpreter','Latex','FontSize',f_size)
%             %         ylabel(['$ \left|\left|e\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
%             ylabel(['$R2\,T$'],'Interpreter','latex','FontSize',f_size)
% 
%             grid on;
%             %             xlim([0,1])
%             xlim([0, xlim_p])
%             ylim(R_plot_lims)
%             hold on;
            subplot(1,N_subplots,3)
            semilogy(time_arr,bound_arr_h,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
            %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
            ylabel(['$Estimator R_h$'],'Interpreter','latex','FontSize',f_size)
            
            %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
            grid on;
            xlim([0, xlim_p])
            ylim([1e-4 1e-0])
            hold on;
            
             subplot(1,N_subplots,4)
            semilogy(time_arr,bound_arr_hv,colour_arr(end-(l_refs-i_ref)),'LineWidth',2)
            %         ylabel(['$ \left|\left|\mathcal{E}\left(t_i\right)\right|\right|_{L_2\left(0,t_i\right)} $'],'Interpreter','latex')
            ylabel(['$Estimator R_{q}$'],'Interpreter','latex','FontSize',f_size)
            
            %title(['$\mathcal{P}^',num2str(2),'$',': $dt=0.1dx^', num2str(1),'$'],'Interpreter','latex')
            xlabel('$t_i$','Interpreter','latex','FontSize',f_size)
            grid on;
            xlim([0, xlim_p])
            ylim([1e-4 1e-0])
            hold on;
            
        end
        hold off;
        set(gcf, 'Position',  [100, 100, 1050, 200])
%          print(gcf,['/Users/gs1511/Desktop/GSialounas_repo_temporal_rec/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_harmonic_u',num2str(round(10*u_init)),'_v',num2str(round(10*v_init)),'_paperrec_poly_tristan.png'],'-dpng','-r600')

            saveas(figgy,['/Users/gs1511/Desktop/GSialounas_swe_repo/figures/fig_',scheme_arr{i_scheme},'plots_1x5_shw_periodic'],'png');%char(Scheme_string)
%             saveas(figgy,['/Users/gs1511/Desktop/GSialounas_BBucket_ApostFD_repo/apostfd/paper/fig_',scheme_arr{i_scheme},'plots_1x5_sin_IC_uniform'],'png');%char(Scheme_string)
        
   
    
    %     sg%title([scheme_arr{i_scheme},' bound'],'Color','red')
    %     sgt.FontSize = 20;
    
end

hold off;

var_out =1;end