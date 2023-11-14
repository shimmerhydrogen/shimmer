function gerg = UtilitiesGERG( x, dim )
	[gerg.Tr,gerg.Dr] = ReducingParametersGERG(x);
	[gerg.Tcx,gerg.Dcx,gerg.Vcx] = PseudoCriticalPointGERG(x,dim);
end
