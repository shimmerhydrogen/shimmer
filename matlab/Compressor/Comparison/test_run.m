function output = test_run(Ne)

Ne
FDM_name = 'FDM-code';
DGM_name = 'DGM-Argo';

demise_file = strcat("../demise_simulations/Ne",int2str(Ne),"/delta_vel_t_all.txt");
T1 = load(demise_file);    
T2 = load("../argo_output/delta_argo.txt");    
% T2 = load("delta_argo_processed_even.txt");    
% T3 = load("delta_argo_processed_odd.txt");    

disp('Loading fine!');

makeFigure(T1, T2, 1, 2, FDM_name, DGM_name,'t', 'delta')
export_fig -dpdf 'Delta_argo_vs_python.pdf';

makeFigure(T1, T2, 1, 3, FDM_name, DGM_name, 't', 'v_w');
export_fig -dpdf 'Vw_argo_vs_python.pdf';


% makeFigure(T1, T2, 1, 3, FDM_name, DGM_name, 't', 'v_w');
% export_fig -dpdf 'Vw_argo_vs_python_even.pdf';
% 
% makeFigure(T1, T2, 1, 2, FDM_name, DGM_name,'t', 'delta')
% export_fig -dpdf 'Delta_argo_vs_python_even.pdf';
% 
% makeFigure(T1, T3, 1, 3, FDM_name, DGM_name, 't', 'v_w')
% export_fig -dpdf 'Vw_argo_vs_python_odd.pdf';
% 
% makeFigure(T1, T3, 1, 2, FDM_name, DGM_name, 't', 'delta')
% export_fig -dpdf 'Delta_argo_vs_python_odd.pdf';

output = 1;
end
function fig_output = makeFigure(T1, T2, idx1, idx2, legend_1, legend_2,x_label, y_label)
close all;
figure;
hold on;
plot( T1(:,idx1), T1(:,idx2),  'r' )%, 'LineWidth',1);
plot( T2(:,idx1), T2(:,idx2),  'b' , 'LineWidth',1);
% set(gca, 'YScale', 'log');
set(gca,'fontsize',16);
set(gcf,'color' ,'w');
box on; grid on;
xlabel(x_label); 
ylabel(y_label);
legend(legend_1, legend_2);
hold off;
fig_output = 1;
end