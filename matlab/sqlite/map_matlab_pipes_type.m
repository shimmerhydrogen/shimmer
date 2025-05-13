function [matlab_pipes_type] = map_matlab_pipes_type(pipe_tab)

num_pipes = size(pipe_tab, 1);

if (num_pipes < 1)
    matlab_pipes_type = [];
    return;
end

matlab_pipes_type = zeros(1, num_pipes);

% pipes_type(p) = 1; Not implemented - No, resistor not supported yet

for p = 1:num_pipes
    if pipe_tab.p_type(p) == 0 % pipe
        matlab_pipes_type(p) = 0;
    elseif pipe_tab.p_type(p) == 1 % compressor
        matlab_pipes_type(p) = 1;
    elseif pipe_tab.p_type(p) == 3 % valve
        matlab_pipes_type(p) = 3;
    elseif pipe_tab.p_type(p) == 2 % reduction and regulation station
        matlab_pipes_type(p) = 2;
    else
        matlab_pipes_type(p) = 0; % DEFAULT is pipe
    end
end

end