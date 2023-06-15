function [klw] =krhogSuperStrate(k0, er, h, TE_TM_flag)
%% EE4620 Assignment 4: [klw] =krhogSuperStrate(k0, er, h, TE_TM_flag)
% Approximated first guess for the leaky wave propagation constant

lambda0 = 2*pi./k0 ;
hbar = h./lambda0 ;

switch TE_TM_flag
    case 'TE'
        % kz0 TE
        kz0 = k0./(pi.*er.*(2.*hbar).^2) .* (2.*pi.*hbar.*er + 1i)./(1 + 1./((2.*pi.*hbar.*er).^2 ) ) ;
    case 'TM'
        % kz0 TM 
        kz0 = k0./(4.*hbar) .* (1+sqrt(1+8.*1i .* hbar./(pi.*er) ) ) ;
end

klw = sqrt(k0.^2-kz0.^2) ;

end