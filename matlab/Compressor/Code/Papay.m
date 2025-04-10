function ZZ=Papay(p,T)
T_cr=190.6;
p_cr=4.595*10^6;
ZZ=1-3.52.*(p./p_cr).*exp(-2.260.*(T./T_cr))+0.274.*(p./p_cr).^2.*exp(-1.878.*(T./T_cr));
end
