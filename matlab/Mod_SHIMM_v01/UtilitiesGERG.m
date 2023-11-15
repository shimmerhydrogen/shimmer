classdef UtilitiesGERG

	properties
	Tr
	Dr
	Tcx
	Dcx
	Vcx
	end

	methods
		function obj = UtilitiesGERG(x, dim)
			[obj.Tr,obj.Dr] = ReducingParametersGERG(x);
			[obj.Tcx,obj.Dcx,obj.Vcx] = PseudoCriticalPointGERG(x,dim);
		end
	end

end
