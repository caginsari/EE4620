function [Jxpws] = Current_PWS(th,phi,k0,L,W,er)
% Assignment 3 EE4620: Current_Uniform(phi,ksw,L,W)

kx = k0 .*sin(th).*cos(phi) ;
ky = k0.* sin(th).*sin(phi) ;
k1 = k0 .* sqrt(er) ;

keq = (k0 + k1) ./2 ;
Jxpws = 2.*keq .* (cos(kx .* L./2)-cos(keq.*L./2) )./((keq.^2 - kx.^2).* sin(keq.*L./2)) .* sinc(ky .* W ./ 2) ; 


end