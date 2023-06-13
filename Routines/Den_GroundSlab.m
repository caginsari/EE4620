function [D] = Den_GroundSlab(k0,er,h,krho,freq,TE_TM_flag)
%% EE4620 Assignment 3: [D] = Den_GroundSlab(k0,er,h,krho,freq,TE_TM_flag)
% Dispersion equation of the Grounded Slab
zeta0 = 120*pi ;
lam0 = 3e8/freq ;
lams = lam0/sqrt(er) ;
% k1 = 2.*pi./lams ;
k1 = k0 * sqrt(er) ;
k2 = k0 ;
kz2 = -1i .*sqrt(-(k2.^2-krho.^2) ) ; 
kz1 = -1i .*sqrt(-(k1.^2-krho.^2) ) ; 
zeta2 = zeta0 ;
zeta1 = zeta0 / sqrt(er) ;

[ZsTM,ZsTE] = TxImpedance(zeta1, kz1, k1) ;
[ZuTM,ZuTE] = TxImpedance(zeta2, kz2, k2) ;

ZdTM = 1i .* ZsTM .* tan( kz1 .* h ) ;
ZdTE = 1i .* ZsTE .* tan( kz1 .* h ) ;

logical = strcmpi('TE',TE_TM_flag) ;

if logical == 1 % if TE 1
    D = ZuTE + ZdTE ;

elseif logical == 0
    D = ZuTM + ZdTM ;
end

end
