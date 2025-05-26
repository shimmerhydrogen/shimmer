
classdef makeFigClass
   properties
      figNo = 1;
      name = 'sthg';
      fig;
      x_label = 'x_label';
      y_label = 'y_label';

   end
   
   %METHODS
   methods        
   
       function obj = makeFigClass(No, name, x_label,y_label )
           if nargin == 4
               
                obj.figNo  = No;
                obj.name   = name;
                obj.x_label = x_label;
                obj.y_label = y_label;
                
                %close(findobj('type','figure','number',No));
                obj.fig = figure(obj.figNo);
                hold on; set(gcf,'color' ,'w'); 
                %set(gca,'xscale', 'log','yscale', 'log')
                
           else 
            disp('Error constructor needs 3 arguments');
           end
       end
    
       function o = add(obj, T, idx_x, idx_y, legend_name, lineSets)
            figure(obj.figNo);           
            if ~exist('lineSets','var')
                plot(T(:, idx_x), T(:, idx_y), 'DisplayName', legend_name, 'Color', rand(1,3), 'LineWidth', 1); 
            else
                plot(T(:, idx_x), T(:, idx_y), lineSets, 'DisplayName', legend_name,'LineWidth', 1); %, 'Color', color); 
            end 

            o = 1;           
       end
      function o = add2(obj, T1,T2, legend_name, lineSets)
            figure(obj.figNo);
            if ~exist('lineSets','var')
                plot(T1, T2, 'DisplayName', legend_name, 'Color', rand(1,3),'LineWidth', 1, 'MarkerFaceColor', 'r', 'MarkerSize', 5);  
            else
                 plot(T1, T2, lineSets, 'DisplayName', legend_name,'LineWidth', 1, 'MarkerFaceColor', 'r', 'MarkerSize', 5); 
            end 
            o = 1;           
      end
      function o = plot3D(obj, T, idx_x, idx_y, idx_z, legend_name, lineSets)
            figure(obj.figNo);           
            if ~exist('lineSets','var')
                plot3(T(:, idx_x), T(:, idx_y), T(:, idx_z), 'DisplayName', legend_name, 'Color', rand(1,3), 'LineWidth', 1); 
            else
                plot3(T(:, idx_x), T(:, idx_y), T(:, idx_z), lineSets, 'DisplayName', legend_name,'LineWidth', 1); %, 'Color', color); 
            end 

            o = 1;           
       end

      function o = scat3D(obj, T, iX,iY,iZ, legend_name)
            figure(obj.figNo);
            
            h = scatter3(T(:,iX),T(:,iY),T(:,iZ),[], T(:,iZ),'filled', 'DisplayName', legend_name); 
            %C = jet(numel(T(:,iZ)));
            %h.CData = fliplr(C);
            % set(gca, 'Colormap', C);
            % colorbar();
            o = 1;           
      end

           function o = plot3D2(obj, TX,TY,TZ, legend_name, lineSets)
            figure(obj.figNo);           
            if ~exist('lineSets','var')
                plot3(TX,TY,TZ, 'DisplayName', legend_name, 'Color', rand(1,3), 'LineWidth', 1); 
            else
                plot3(TX,TY,TZ, lineSets, 'DisplayName', legend_name,'LineWidth', 1); %, 'Color', color); 
            end 

            o = 1;           
       end

     function o = scat3D2(obj, TX,TY,TZ, legend_name)
            figure(obj.figNo);
   
            h = scatter3(TX,TY,TZ,[],'filled', 'DisplayName', legend_name);  
            C = jet(numel(T(:,iZ)));
            h.CData = fliplr(C);
%             set(gca, 'Colormap', C)
%             colorbar()
            o = 1;           
      end

      
      function o = change2log(obj, setLogX, setLogY)
            figure(obj.figNo);
            if(setLogX)
                disp('Seeting log sclae in X')
                set(gca, 'XScale', 'log')
            elseif(setLogY)
                disp('Seeting log sclae in Y')
                set(gca, 'YScale', 'log')
            else
                disp('Warning: Log Scale not set, check inputs')
            end 
            o = 1;
       end
       function o = endFig(obj)
            figure(obj.figNo);
            legend('-DynamicLegend');
            legend('AutoUpdate','off');            
            set(gca,'xgrid', 'on', 'ygrid', 'on','box','on');
            set(gca,'fontsize',16, 'FontName', 'Arial');
            set(gca,'XMinorTick','on','YMinorTick','on');
            xlabel(obj.x_label, 'fontsize',16, 'FontName', 'Arial');
            ylabel(obj.y_label, 'fontsize',16, 'FontName', 'Arial');
            o = 1;
       end
       
       function export(obj, root)
          figure(obj.figNo);
            path = strcat(root, obj.name);
            disp(path)
           export_fig(path, '-pdf');
           close(obj.fig);
           
       end
   end

end


