function [rho,drho] = computedensity(P,cf,rho0,P0)
rho = rho0.*exp(cf.*(P-P0));
drho= cf*rho0.*exp(cf.*(P-P0));

end