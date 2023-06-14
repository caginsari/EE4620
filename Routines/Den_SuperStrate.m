function [D] = Den_SuperStrate(k0,er,h,krho)
%% EE4620 Assignment 4: [D] = Den_SuperStrate(k0,er,h,krho,freq,TE_TM_flag)
% krho needs to be the guess point of the propagation constant for TE1 or
% TM1 mode.
% Dispersion equation of the SuperStrate
zeta0 = 120*pi ;
kz0 = -1i .*sqrt(-(k0.^2-krho.^2) ) ; 
n = 1 ;


D = zeta0./ er .* 1i .* zeta0 .* (kz0.*h - n.*pi) ;

end