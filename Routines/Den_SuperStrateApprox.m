function [D] = Den_SuperStrateApprox(k0,er,h,krho,TE_TM_flag)
%% EE4620 Assignment 4: [D] = Den_SuperStrate(k0,er,h,krho,freq,TE_TM_flag)
% Approximated Dispersion Equation. 
% krho needs to be the guess point of the propagation constant for TE1 or
% TM1 mode.
% Dispersion equation of the SuperStrate
zeta0 = 120*pi ;
kz0 = -1i .*sqrt(-(k0.^2-krho.^2) ) ; 
n = 1 ;

switch TE_TM_flag
    case 'TM'
        % Z0_TM
        [Z0,~] = TxImpedance(zeta0,kz0,k0) ;
    case 'TE'
        % Z0_TE
        [~,Z0] = TxImpedance(zeta0,kz0,k0) ;
end


D = zeta0./ er + 1i .* Z0 .* (kz0.*h - n.*pi) ;

end