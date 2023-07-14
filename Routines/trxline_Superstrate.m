function [vte,vtm,ite,itm,ks,kz0] = trxline_Superstrate(k0,zeta0,er,h,hs,krho,Layer,z,freq)
% EE4620 Assignment 1: te and tm parameter calculation function for a
% transmission line involving 3 layers; case for homework : 
% ground plane-> air->thin substrate->air
lam0 = 3e8./freq ;
lams = lam0./sqrt(er) ;
ks = 2.*pi./lams ;
kz0 = -1i .*sqrt(-(k0.^2-krho.^2) ) ; 
kzs = -1i .*sqrt(-(ks.^2-krho.^2) ) ; 
zetas = zeta0 / sqrt(er) ;

[Z3tm,Z3te] = TxImpedance(zeta0, kz0, k0) ; % Z0 in the instruction slides
[Z2tm,Z2te] = TxImpedance(zetas, kzs, ks) ; % Zs in the instruction slides
Z1in_te =  TxLineInputImpedance(Z2te,Z3te,kzs,hs) ;
Z1in_tm =  TxLineInputImpedance(Z2tm,Z3tm,kzs,hs) ;

[gamma1te,gamma1tm] = TxLineReflectionCoeff(Z1in_te,Z1in_tm,Z3te,Z3tm) ;

[gamma2te,gamma2tm] = TxLineReflectionCoeff(Z3te,Z3tm,Z2te,Z2tm) ;

L1 = strcmpi('Layer1',Layer) ;
L2 = strcmpi('Layer2',Layer) ;
L3 = strcmpi('Layer3',Layer) ;

if L1 == 1
    [vte,vtm,ite,itm] = Layer1Superstrate(gamma1te,gamma1tm,h,z,kz0,Z3te,Z3tm) ;

elseif L2 == 1
    [vte,vtm,ite,itm] = Layer2Superstrate(kz0,kzs,h,hs,Z2te,Z2tm,gamma1te,gamma1tm,gamma2te,gamma2tm) ;

elseif L3 == 1
    [vte,vtm,ite,itm] = Layer3Superstrate(gamma1te,gamma1tm,gamma2te,gamma2tm,kz0,kzs,h,hs,z,Z3te,Z3tm) ;
else
    error(['\nInput must be one of three superstrate layers: \n Layer1 : air \n Layer2 : thin dielectric \n ' ...
        'Layer3 : air \n your input: %s.'],Layer) ;
end


end