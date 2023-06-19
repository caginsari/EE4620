function [FTMeq] = FTCurrentDoubleSlot(kx,ky,keq,W,L,d,slotType)
%% EE4620 Assignment 4 Question 3.
% FT of Double slot antenna with spacing d.
switch slotType
    case 'single'
        AF = 1 ;
    case 'double'
        AF = 2 .* cos(kx .* d./2) ;
end

% FTMeq = -4.*keq.*cos(ky.*L./2)./(ky.^2 - keq.^2) .* sinc(kx.*W./(2.*pi) ) .* AF;
FTMeq = 2.*keq.*(cos(ky.*L./2) - cos(keq.*L./2) )./((keq.^2-ky.^2).*sin(keq.*L./2)).* sinc(kx.*W./(2.*pi) ) .* AF;

end