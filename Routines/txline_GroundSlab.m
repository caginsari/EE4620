function [vtm, vte, itm, ite, kz1,kz2] = txline_GroundSlab(k0,zeta0,er,z,h,krho,freq)
% EE4620 Assignment 1 [vtm, vte, itm, ite] = txline_GroundSlab(k0,zeta0,er,z,h,krho)

lam0 = 3e8/freq ;
lams = lam0/sqrt(er) ;
k1 = 2.*pi./lams ;
k2 = k0 ;
kz2 = -1i .*sqrt(-(k2.^2-krho.^2) ) ; 
kz1 = -1i .*sqrt(-(k1.^2-krho.^2) ) ; 
zeta2 = zeta0 ;
zeta1 = zeta0 / sqrt(er) ;

[ZsTM,ZsTE] = TxImpedance(zeta1, kz1, k1) ;
[ZuTM,ZuTE] = TxImpedance(zeta2, kz2, k2) ;

ZdTM = 1i .* ZsTM .* tan( kz1 .* h ) ;
ZdTE = 1i .* ZsTE .* tan( kz1 .* h ) ;

if z<h
    % in slab
    [vtm,itm] = VandI_Slab(ZuTM,ZdTM,z,h,kz1) ;
    [vte,ite] = VandI_Slab(ZuTE,ZdTE,z,h,kz1) ;
else
    % in air
    [vtm,itm] = VandI_Air(ZuTM,ZdTM,  ZuTM , z,h,kz2) ;
    [vte,ite] = VandI_Air(ZuTE,ZdTE, ZuTE , z,h,kz2 ) ;
end


end