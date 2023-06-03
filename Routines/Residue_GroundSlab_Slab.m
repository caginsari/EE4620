function [Vr,Ir] = Residue_GroundSlab_Slab(k0,krho,er,h,z,TE_TM_flag,Dprime)
%% EE4620 Assignment 3: [Vr,Ir] = Residue_GroundSlab_Slab(k0,er,h,z,TE_TM_flag,Dprime)
% Dprime must be the one obtained from the Newton-Raphson Method in function IterativeMethod  
zeta0 = 120*pi ;
% ks = 2.*pi./lams ;
ks = k0 * sqrt(er) ;% medium 1 is dielectric 
kz0 = -1i .*sqrt(-(k0.^2-krho.^2) ) ; % kz in air is given as kz0 
kzs = -1i .*sqrt(-(ks.^2-krho.^2) ) ; 
zetas = zeta0 / sqrt(er) ;

[ZsTM,ZsTE] = TxImpedance(zetas, kzs, ks) ;
[ZuTM,ZuTE] = TxImpedance(zeta0, kz0, k0) ;
ZuTE(isinf(ZuTE) ) = 0 ;

ZdTM = 1i .* ZsTM .* tan( kzs .* h ) ;
ZdTE = 1i .* ZsTE .* tan( kzs .* h ) ;

if 1 == strcmpi(TE_TM_flag,'TE')
    Vr = ZuTE.* ZdTE ./ Dprime .* sin(kzs.*z)./sin(kzs.*h) ;
    Ir = ZuTE.* ZdTE ./ ( Dprime.* ZsTE ) .* 1i.*cos(kzs.*z)./sin(kzs.*h) ;
elseif 1 == strcmpi(TE_TM_flag,'TM')
    Vr = ZuTM.* ZdTM ./ Dprime .* sin(kzs.*z)./sin(kzs.*h) ;
    Ir = ZuTM.* ZdTM ./ ( Dprime.* ZsTM ) .* 1i.*cos(kzs.*z)./sin(kzs.*h) ;
else
    error('The input must be: \n 1.TE \n 2.TM \n Your input is %s',TE_TM_flag) ;
end

end