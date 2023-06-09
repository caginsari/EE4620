function [I_phi] = PhiInt_Uniform(ksw,L,W)
%% [I_phi] = PhiInt_Uniform(phi,er,k0,ksw,l,w)
% Computes the I_phi for uniform current distribution.
% Note: I have doubts about the PhiInt_PWS function, that's why it is not
% implemented.
phi = linspace(eps,2.*pi,202) ;
dphi = phi(2) - phi(1) ;

kx = ksw .*cos(phi) ;
ky = ksw.* sin(phi) ;
% k1 = k0 .* sqrt(er) ;

Jx = L.*sinc(kx.*L ./ 2) .* sinc(ky .* W ./ 2) ;

% keq = (k0 + k1) ./2 ;
% Jx = FTCurrent( keq ,kx, ky, l, w ) ;

I_phi = sum(abs(Jx).^2 .* cos(phi).^2 ) .*dphi ;


end
