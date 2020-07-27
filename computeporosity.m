function [phi,dphi] = computeporosity(P,cr,phi0,P0)
phi = phi0.*exp(cr.*(P-P0));
dphi= cr*phi0.*exp(cr.*(P-P0));

end