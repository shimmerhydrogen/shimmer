function [] = sql_populate_stations_gass_molar_fractions(db_path, graph)

num_nodes = size(graph.Nodes.Nodes_ID, 1);

if (num_nodes < 2)
    return;
end

num_variables = 2;
stations_tab_variables_name = ["s_number", "frac_CH4"];

%% Gas molar Fraction
% Now all the gas fractions are 100% CH4

stations_gas_molar_fractions = cell(num_nodes, num_variables);
for s = 1:num_nodes
    stations_gas_molar_fractions{s, 1} = graph.Nodes.Nodes_ID(s);
    stations_gas_molar_fractions{s, 2} = 1.0;
end

stations_gas_molar_fractions_tab = cell2table(stations_gas_molar_fractions, ...
    "VariableNames", stations_tab_variables_name );

conn = sqlite(db_path, 'connect');
sqlwrite(conn, "gas_molar_fractions", stations_gas_molar_fractions_tab);
close(conn);

end