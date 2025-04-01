function pass = compare()

fp = makeFigClass(1, "V_h40", "x" , "V"); 

readAndDraw("flux"); 

pass = true;
end

function pass = readAndDraw(control)
    close all;
    %data
    t = [0:24]; 
    cmp_node_in = 14;
    cmp_node_out = 15; 
    inj_node = 12;
    remi_node = 1;
    num_nodes = 15;
    num_pipes = 17; 
    comp_num = 16 + num_nodes;

    % Load results
    cpp_dir = "/home/karol/Documents/UNIVERSITA/POLITO/shimmer/shimmer++/build/unit_tests";

    disma = load(strcat(cpp_dir , "/var_time_DISMA_", control, ".dat"))';
    denerg = load(strcat("var_time_DENERG_", control, ".dat"));

    disp(" * Relative Error")
    (abs(denerg - disma))./denerg
    
    %Built figure
    fp = makeFigClass(1, strcat(control, "_COMPRESSOR_PRESSURE") , "Time [h]" , "Pressure [bar]") 
    fg = makeFigClass(2, strcat(control, "_COMPRESSOR_FLOW_RATE"), "Time [h]" , "Flow rate [kg/s]") 
    
    % Pressure P_in
    fp.add2(t, denerg(cmp_node_in, :)/1e5, "Ref. Inlet", 'or');
    fp.add2(t, denerg(cmp_node_out, :)/1e5, "Ref. Outlet", '*r');
    fp.add2(t, disma(cmp_node_in,:)/1e5, "shimmer++. Inlet", 'b');
    fp.add2(t, disma(cmp_node_out,:)/1e5, "shimmer++. Outlet", 'k');
    
    fp.endFig();
    fp.export("");
    
    % Flow rate P_in
    fg.add2(t, denerg(comp_num, :), "Ref.", '*r');
    fg.add2(t, disma(comp_num,:), "shimmer++", 'b');
    
    fg.endFig();
    fg.export("");

    fip = makeFigClass(3, strcat(control, "_INLET_ST_PRESSURE") , "Time [h]" , "Pressure [bar]")     

    fip.add2(t, denerg( inj_node, :)/1e5, "Ref. Entry-L", 'or');
    fip.add2(t, denerg(remi_node, :)/1e5, "Ref. Entry-p", 'or');
    fip.add2(t, disma(  inj_node, :)/1e5, "shimmer++. Entry-l", 'b');
    fip.add2(t, disma( remi_node, :)/1e5, "shimmer++. Entry-p", 'k');
    fip.endFig();
    fip.export("");

    offset = num_nodes + num_pipes;
    fig = makeFigClass(3, strcat(control, "_INLET_ST_FLOW_RATE") , "Time [h]" ,  "Flow rate [kg/s]")    
    fig.add2(t, denerg(inj_node + offset, :), "Ref. Entry-L", 'or');
    fig.add2(t, denerg(remi_node+ offset, :), "Ref. Entry-p", 'or');
    fig.add2(t, disma(inj_node  + offset, :), "shimmer++. Entry-l", 'b');
    fig.add2(t, disma(remi_node + offset, :), "shimmer++. Entry-p", 'k');
    fig.endFig();
    fig.export("");


    fe = makeFigClass(3, strcat(control, "_ERROR") , "Time [h]" ,  "Relative Error")    

    error = abs(denerg - disma)./denerg; 
    error(find(isnan(error)>0)) = 0;
    fe.add2(t, sqrt(sum(error(1:(num_pipes+num_nodes),:).^2,1)), "Error", 'or');
    fe.endFig();
    fe.change2log(false, true);
    fe.export("");
    pass = true;
end 