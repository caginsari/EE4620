function [Va,Ia] = VandI_Air(Zu,Zd,Zs,z,h,kz0)
% Assignment 1 EE4620 : Ground slab 
Va = ( Zu.*Zd ./ (Zu+Zd) ) .* exp(1i.*kz0.*h) .* exp(-1i.*kz0.*z) ;

Ia = 1./ Zs .* (  Zu.*Zd ./ (Zu+Zd) ) .* exp(1i.*kz0.*h).* exp(-1i.*kz0.*z) ;

end