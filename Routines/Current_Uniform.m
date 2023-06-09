function [JxUniform] = Current_Uniform(th,phi,k0,L,W)
% Assignment 3 EE4620: Current_Uniform(phi,ksw,L,W)

kx = k0 .*sin(th).*cos(phi) ;
ky = k0.* sin(th) .*sin(phi) ;
% k1 = k0 .* sqrt(er) ;

JxUniform = L.*sinc(kx.*L ./ 2) .* sinc(ky .* W ./ 2) ;


end