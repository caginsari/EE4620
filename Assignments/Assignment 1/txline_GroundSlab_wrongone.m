function [vtm, vte, itm, ite, k1] = txline_GroundSlab_wrongone(k0,zeta0,er,z,h,krho,freq)
% EE4620 Assignment 1 [vtm, vte, itm, ite] = txline_GroundSlab(k0,zeta0,er,z,h,krho)
lam0 = 3e8/freq ;
lams = lam0/(er) ;
k1 = 2.*pi./lams ;
k2 = k0 ;
kz2 = -1i .*sqrt(-(k2.^2-krho.^2)) ; 
kz1 = -1i .*sqrt(-(k1.^2-krho.^2)) ; 
zeta2 = zeta0 ;
zeta1 = sqrt((pi*4e-7)/(8.85e-12 *er)) ;

ZiTM =@(zetai,kzi,ki) zetai .* kzi ./ki ;
ZiTE =@(zetai,kzi,ki) zetai .* ki ./kzi ;
Zd =@(kzi,h,Z) 1i .* Z .* tan(kzi .* h) ;

Vs =@(Zu,Zd,kzs,z,h) ( Zu.*Zd ./ (Zu+Zd) ) .* sin(kzs.*z) ./ sin(kzs.*h) ;

Is =@(Zu,Zd,Zs,kzs,z,h) 1./ Zs .* ( Zu.*Zd ./ (Zu+Zd) ) .* 1i.*cos(kzs.*z) ./ sin(kzs.*h) ;

Va =@(Zu,Zd,kz0,z,h) ( Zu.*Zd ./ (Zu+Zd) ) .*exp(1i.*kz0.*h).*exp(-1i.*kz0.*z) ;

Ia =@(Zu,Zd,Zs,kz0,z,h) 1./ Zs .* (  Zu.*Zd ./ (Zu+Zd) ) .*exp(1i.*kz0.*h).*exp(-1i.*kz0.*z) ;


ZuTE = ZiTE(zeta2,kz2,k2) ;
ZuTM = ZiTM(zeta2,kz2,k2) ;

ZsTE = ZiTE(zeta1,kz1,k1) ;
ZsTM = ZiTM(zeta1,kz1,k1) ;
ZdTE = Zd(kz1,h,ZsTE) ;
ZdTM = Zd(kz1,h,ZsTM) ;


if z<h
    % in slab
    vte = Vs(ZuTE,ZdTE,kz1,z,h) ;
    ite = Is(ZuTE,ZdTE,ZsTE,kz1,z,h) ;
    vtm = Vs(ZuTM,ZdTM,kz1,z,h) ;
    itm = Is(ZuTM,ZdTM,ZsTM,kz1,z,h) ;
else
    % in air
    vte = Va(ZuTE,ZdTE,kz2,z,h) ;
    ite = Ia(ZuTE,ZdTE,ZuTE,kz2,z,h) ;
    vtm = Va(ZuTM,ZdTM,kz2,z,h) ;
    itm = Ia(ZuTM,ZdTM,ZuTM,kz2,z,h) ;
end


end