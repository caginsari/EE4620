function [Psw] = PswTMelem(k0,er,h,ksw,TE_TM_flag)
%% EE4620 Assignment 3: [Psw] = PswTMelem(k0,er,h,ksw)
% Calculates the phi component of the current necessary for numerical
% computation of the power in the surface wave.
% Iphi_delt: is the phi component of the current from elementary current
% source

zeta0 = 120*pi ;
ks = k0 * sqrt(er) ;% medium 1 is dielectric 
kz0 = -1i .*sqrt(-(k0.^2-ksw.^2) ) ; % kz in air is given as kz0 
kzs = -1i .*sqrt(-(ks.^2-ksw.^2) ) ; 
zetas = zeta0 / sqrt(er) ;

[ZsTM,ZsTE] = TxImpedance(zetas, kzs, ks) ;
[ZuTM,ZuTE] = TxImpedance(zeta0, kz0, k0) ;
ZuTE(isinf(ZuTE) ) = 0 ;

ZdTM = 1i .* ZsTM .* tan( kzs .* h ) ;
ZdTE = 1i .* ZsTE .* tan( kzs .* h ) ;

deltk = k0 ./500 ;
Dkrhop = ksw + deltk./2 ;
Dkrhom = ksw - deltk./2 ;

Dp = Den_GroundSlab(k0,er,h,Dkrhop,f,TE_TM_flag) ;
Dm = Den_GroundSlab(k0,er,h,Dkrhom,f,TE_TM_flag) ;
Dprime = (Dp-Dm) ./ deltk ;

if 1 == strcmpi('TM',TE_TM_flag)
    Zs = ZsTM ;
    Zu = ZuTM ; 
    Zd = ZdTM ; 

elseif 1 == strcmpi('TE',TE_TM_flag)
    Zs = ZsTE ;
    Zu = ZuTE ; 
    Zd = ZdTE ; 
else
    error('The input must be one of the following:\ n 1. TE \n 2. TM \n Your input is %s',TE_TM_flag)
end 

Cons = zeta0 .* Zu.* Zd ./(k0.* Dprime) ;
Iz_slab = Cons .* h./(2.*er.*Zs .*sin(kzs.*h) ) .*h./2 .*(1+sinc(2.*kzs.*h./pi) ) ;





end
