function [Vr,Ir] = Residue_GroundSlab_Air(k0,er,h,krho,z,TE_TM_flag,Dprime)
%% EE4620 Assignment 3: [Vr,Ir] = Residue_GroundSlab_Air(k0,er,h,z,TE_TM_flag,Dprime)
% Dprime must be the one obtained from the Newton-Raphson Method in function IterativeMethod 
% or it can be calculated using ksw.
zeta0 = 120*pi ;
% ks = 2.*pi./lams ;
ks = k0 * sqrt(er) ;% medium 1 is dielectric 
kz0 = -1i .*sqrt(-(k0.^2-krho.^2) ) ; % kz in air is given as kz0 
kzs = -1i .*sqrt(-(ks.^2-krho.^2) ) ; 
zetas = zeta0 / sqrt(er) ;

[ZsTM,ZsTE] = TxImpedance(zetas, kzs, ks) ;
[ZuTM,ZuTE] = TxImpedance(zeta0, kz0, k0) ;

ZdTM = 1i .* ZsTM .* tan( kzs .* h ) ;
ZdTE = 1i .* ZsTE .* tan( kzs .* h ) ;

if 1 == strcmpi(TE_TM_flag,'TE')
    Vr = ZuTE.* ZdTE ./ Dprime .* exp(1i.*kz0.*h).*exp(-1i.*kz0.*z) ;
    Ir = ZuTE.* ZdTE ./ ( Dprime.* ZuTE ) .* exp(1i.*kz0.*h).*exp(-1i.*kz0.*z) ;
elseif 1 == strcmpi(TE_TM_flag,'TM')
    Vr = ZuTM.* ZdTM ./ Dprime .* exp(1i.*kz0.*h).*exp(-1i.*kz0.*z) ;
    Ir = ZuTM.* ZdTM ./ ( Dprime.* ZuTM ) .* exp(1i.*kz0.*h).*exp(-1i.*kz0.*z) ;
else
    error('The input must be: \n 1.TE \n 2.TM \n Your input is %s',TE_TM_flag) ;
end

end