function Omega = Resistance(dxe, dt, pm, p, G, AA, ADP, MOLEFRAC, Tb, epsi, DD, cc2b)

Omega.Ri = 2.0*dxe.*pm./(AA*dt)./(abs(ADP)*p); % Inertia Resistance

[lambda Re viscosity] =  friction(Tb,epsi,G,DD, MOLEFRAC);
Omega.Rf = 16.*lambda.*cc2b.*dxe./(DD.^5.*pi.^2)./(abs(ADP)*p); % Fluid-dynamic Resistance

rr = (2*Omega.Rf .*(abs(G)) + Omega.Ri); % composite resistance linearized problem (R)
Omega.R = sparse(diag(rr));       % transformed into sparse diagonal matrix
end
