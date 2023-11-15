function ADP = matrixADP(Aminus, Aplus, cc2b, delH, dimn)
	% A.2 FIND ADP
	s = 2 * 9.81 * delH./cc2b;
	Aminus_s  = Aminus.*repmat(exp(s),1,dimn)';
	Aminus_s1 = Aminus.*repmat(exp(s/2),1,dimn)';
	ADP = (-Aminus_s1 + Aplus)'; % KCommit:  Ag't pag 27-28

end
