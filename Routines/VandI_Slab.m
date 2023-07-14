function [Vs, Is ] = VandI_Slab(Zu,Zd,z,h,kzs)
% Assignment 1 EE4620 : Ground slab 
Vs = ( Zu.*Zd ./ (Zu+Zd) ) .* sin(kzs.*z) ./ sin(kzs.*h) ;

Is =( Zu.*Zd ./ (Zu+Zd) ) .* 1i.*cos(kzs.*z) ./ sin(kzs.*h) ;

end