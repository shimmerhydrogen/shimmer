function [f Re viscosity]=friction(Tm,epsi,G,D,MOLEFRAC)
G=abs(G);                           % Absolute value mass flow rate
viscosity=Viscositycalculator(MOLEFRAC,Tm);         % Dynamic viscosity [Pa*s]
Re=G.*D./ (pi.*(D/2).^2)./viscosity;       % Reynolds number
 for j=1:length(Re)                  % Friction factor for each pipe: numerical solution of Colebrook-White equation
	%     if Re(j)<2000
	%         f(j,1)=64/Re(j);
	%     elseif Re(j)>=2000 && Re(j)<= 4000
	%           f_lam(j,1)=64/Re(j)  ;% Fow low Re numbers
	%           f_turb(j,1)=fzero(@(f)1/(f^0.5)+2.0*log10(epsi(j)/3.7/D(j)+2.51/Re(j)/(f^0.5)),0.05);
	%
	%           f(j,1)=(f_lam(j,1)*(4000-Re(j))+f_turb(j,1)*(Re(j)-2000))/2000;
	%     else
	% %          f(j,1)=fzero(@(f)1/(f^0.5)+0.86*log(epsi(j)/3.7/D(j)+2.51/Re(j)/(f^0.5)),0.05);   % For high Re numbers
	%          f(j,1)=fzero(@(f)1/(f^0.5)+2.0*log10(epsi(j)/3.7/D(j)+2.51/Re(j)/(f^0.5)),0.05);   % For high Re numbers
	%
	%     end
	%
	% friction factor numeric
	a=1./(1+(Re(j)./2720).^9);
	b=1./(1+(Re(j)./(160.*D(j)./epsi(j))).^2);

	f(j,1)=((Re(j)./64).^a*(1.8.*log10(Re(j)./6.8)).^(2.*(1-a).*b).*(2.0.*log10(3.7*D(j)./epsi(j))).^(2.*(1-a).*(1-b))).^-1;

end
