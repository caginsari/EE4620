function [D] = Den_SemiInfNumeric(k0,er,h,krho,TE_TM_flag)
%% EE4620 Assignment 4: [D] = Den_SemiInfNumeric(k0,er,krho,TE_TM_flag)
% Numerical Dispersion equation of the Semi Infinite superstrate
lam0 = 2*pi./k0 ;
lams = lam0/sqrt(er) ;
ks = 2.*pi./lams ;
zeta0 = 120*pi ;

kz0 = 1i .*sqrt(-(k0.^2-krho.^2) ) ; 
kzs = 1i .*sqrt(-(ks.^2-krho.^2) ) ; 
zetas = zeta0 / sqrt(er) ;

[Z0TM,Z0TE] = TxImpedance(zeta0, kz0, k0) ; % Z0 in the instruction slides
[Z1TM,Z1TE] = TxImpedance(zetas, kzs, ks) ; % Zs in the instruction slides

switch TE_TM_flag 
    case 'TE'
    D = Z1TE + 1i .* Z0TE.* tan(kz0.*h);

    case 'TM' 
    D = Z1TM + 1i .* Z0TM.* tan(kz0.*h);

    case 'TM0' 
    D = Z1TM + 1i .* Z0TM.* tan(kz0.*h);
end

end

