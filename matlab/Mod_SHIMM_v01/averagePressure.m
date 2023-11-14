function pm =  averagePressure(p_in, p_out);
% Eq pag 28 pm = (p_in^2 + p_in*p_out + p_out^2) /(p_in+p_out)
	pm =(p_in.^2 + p_in*p_out + p_out.^2) /(p_in + p_out);
end

