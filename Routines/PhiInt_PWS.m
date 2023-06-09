function [I_phi] = PhiInt_PWS(er,k0,ksw,L,W)
%% [I_phi] = PhiInt_Uniform(phi,er,k0,ksw,l,w)
% Computes the I_phi for uniform current distribution.
% Note: I have doubts about the PhiInt_PWS function, that's why it is not
% implemented.
phi = linspace(eps,2.*pi,202) ;
dphi = phi(2) - phi(1) ;

kx = ksw .*cos(phi) ;
ky = ksw.* sin(phi) ;
k1 = k0 .* sqrt(er) ;

keq = (k0 + k1) ./2 ;
Jx = FTCurrent(keq,kx,ky,L,W) ;
% % Jx = 2.*keq .* (cos(kx .* L./2)-cos(keq.*L./2) )./((keq.^2 - kx.^2).* sin(keq.*L./2)) .* sinc(ky .* W ./ 2) ; 

I_phi = sum(abs(Jx).^2 .* cos(phi).^2 ) .*dphi ;


end