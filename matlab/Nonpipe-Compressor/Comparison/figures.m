% Tp = load("../python_output/Temperature_t0.1.txt");    
% Ta = load("../argo_output/myDomain_Heat_Immersed/Temperature_bottom_processed.t100.dat");    
% 
% Tp2 = load("../python_output/Temperature_t0.18.txt");    
% Ta2 = load("../argo_output/myDomain_Heat_Immersed/Temperature_bottom_processed.t180.dat");    
% 
% figure
% hold on;
% plot( Tp(:,1),  Tp(:,3),  'k'  , 'LineWidth', 2);
% plot( Ta(:,1),  Ta(:,3),  '+r' , 'MarkerSize', 8);
% plot( Tp2(:,1), Tp2(:,3), '-.k', 'LineWidth', 2); 
% plot( Ta2(:,1), Ta2(:,3), '*r' , 'MarkerSize', 8);
% 
% legend('t = 0.10 with 1D-FDM', 't = 0.10 with Argo-DGM', 't = 0.18 with 1D-FDM','t = 0.18 with Argo-DGM');
% xlabel('x','FontSize',16); ylabel("T",'FontSize', 16);
% axis([0 0.005 200 1800]);
% set(gcf,'color' ,'w');
% set(gca,'fontsize',16);
% box on; grid on;
% export_fig -dpdf 'XvsTemperature.pdf'

Dp = load("../python_output/delta_t_all.txt");    
Da = load("../argo_output/myDomain_Heat_Immersed/delta_argo.dat");   
figure
hold on;
plot(Dp(:,1), Dp(:,2), 'k', 'LineWidth', 2);
plot(Da(:,1), Da(:,2), '+r','MarkerSize', 8);
%plot( Da(:,1), Da(:,2), '+r', 'MarkerSize', 8);
axis([0 0.2 0 3.5e-5]);
legend('1D-FDM','Argo-DGM');
set(gcf,'color' ,'w')
xlabel('t','FontSize',14); ylabel("\delta",'FontSize', 14);
set(gca,'fontsize',16);
box on; grid on;
export_fig -dpdf 'tvsDelta.pdf'