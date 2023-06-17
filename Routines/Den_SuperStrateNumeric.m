function [D] = Den_SuperStrateNumeric(k0,er,h,hs,krho,TE_TM_flag)
%% EE4620 Assignment 4: [D] = Den_SuperStrateNumeric(k0,er,hs,krho,TE_TM_flag)
% Numerical Dispersion equation of the Super Strate

zeta0 = 120*pi ;
ks = k0 * sqrt(er) ;
kz0 = 1i .*sqrt(-(k0.^2-krho.^2) ) ; 
kzs = 1i .*sqrt(-(ks.^2-krho.^2) ) ; 
zetas = zeta0 ./ sqrt(er) ;

[ZsTM,ZsTE] = TxImpedance(zetas, kzs, ks) ;
[Z0TM,Z0TE] = TxImpedance(zeta0, kz0, k0) ;

Z1in_te =  TxLineInputImpedance(ZsTE,Z0TE,kzs,hs) ;
Z1in_tm =  TxLineInputImpedance(ZsTM,Z0TM,kzs,hs) ;

logical = strcmpi('TE',TE_TM_flag) ;

if logical == 1 % if TE 1
    D = Z1in_te + 1i .* Z0TE.* tan(kz0.*h);

elseif logical == 0
    D = Z1in_tm + 1i .* Z0TM.* tan(kz0.*h);
end

end

