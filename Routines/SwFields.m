function [Erho,Hphi,Ez,h] = SwFields(k0,ksw,er,Vr,Ir,rho,phi,l,w)
%% EE4620 Assignment 3 : [Erho,Ephi,Ez] = SwFields(k0,ksw,er,Vr,Ir,TE_TM_flag,rho,phi)
% Calculates the field components of the surface waves
% Vr and Ir could be both te and tm. For field components in TE and TM
% basically provide Vr and Ir calculated by using TE or TM components
% er: (relative permittivity) according to the er the field is calculated either in air or other
% dielectric material

zeta = 120.*pi / sqrt(er) ;
k = k0 * sqrt(er) ;

keq = (k + k0)/2 ;

kx = ksw .* cos(phi) ;
ky = ksw .* sin(phi) ;
JxFT = FTCurrent(keq,kx,ky,l,w) ; 
C = 1i.*sqrt(ksw./(2.*pi)).*exp(1i.*pi./4) ;

Erho = Vr .* JxFT .* C .*cos(phi) .*exp(-1i.* ksw .*rho)./ sqrt(rho) ;
Ez = -zeta.*ksw./k .* Ir .* JxFT .* C .* cos(phi) .*exp(-1i.* ksw .*rho)./ sqrt(rho) ;
Hphi = Ir .* JxFT .* C .*cos(phi) .*exp(-1i.* ksw .*rho)./ sqrt(rho) ;

% Electrical Length h
h = JxFT./(2.*keq) .* cos(phi) .* 2.* keq;

end
