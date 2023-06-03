function [vte,vtm,ite,itm,ks] = trxline_semi_inf_Superstrate(k0,zeta0,er,h,krho,z,freq,Layer)
% EE4620 Assignment 1:  semi infinite trxline voltage and current
% parameters

lam0 = 3e8/freq ;
lams = lam0/sqrt(er) ;
ks = 2.*pi./lams ;

kz0 = -1i .*sqrt(-(k0.^2-krho.^2) ) ; 
kzs = -1i .*sqrt(-(ks.^2-krho.^2) ) ; 
zetas = zeta0 / sqrt(er) ;

[Z0tm,Z0te] = TxImpedance(zeta0, kz0, k0) ; % Z0 in the instruction slides
[Z1tm,Z1te] = TxImpedance(zetas, kzs, ks) ; % Zs in the instruction slides

[gamma1te,gamma1tm] = TxLineReflectionCoeff(Z1te,Z1tm,Z0te,Z0tm) ;

L1 = strcmpi('Layer1',Layer) ;
L2 = strcmpi('Layer2',Layer) ;

if L1 == 1
    [vte,vtm,ite,itm] = Layer1SemiInfSuperstrate(gamma1te,gamma1tm,kz0,h,z,Z0te,Z0tm) ; % Layer 1 is air 

elseif L2 == 1
    [vte,vtm,ite,itm] = Layer2SemiInfSuperstrate(gamma1te,gamma1tm,kz0,h,z,Z1te,Z1tm,kzs) ; % Layer 2 is semi infinite dielectric 
else
    error(['\nInput must be one of three superstrate layers: \n Layer1 : air \n Layer2 : Semi infinite dielectric \n ' ...
        'Layer3 : air \n your input: %s.'],Layer) ;
end

end