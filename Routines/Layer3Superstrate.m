function [vte,vtm,ite,itm]=Layer3Superstrate(gamma1te,gamma1tm,gamma2te,gamma2tm,kz0,kzs,h,hs,z,Z0te,Z0tm)
%EE4620 Assignment 1: Calculates the transmission line parameters for the
%3rd layer of the superstrate

vteprogte = (1+gamma1te).*(1+gamma2te) .*exp(1i.*kz0.*h).*exp(-1i.*kzs.*hs).*exp(1i.*kz0.*(h+hs) ) ./...
    ( (gamma1te + exp(2.*1i.*kz0.*h) ).* (1 + gamma2te .* exp(-2.*1i.*kzs.*hs) ) ) ;

vtmprogtm = (1+gamma1tm).*(1+gamma2tm) .*exp(1i.*kz0.*h).*exp(-1i.*kzs.*hs).*exp(1i.*kz0.*(h+hs) ) ./...
    ( (gamma1tm + exp(2.*1i.*kz0.*h) ).* (1 + gamma2tm .* exp(-2.*1i.*kzs.*hs) ) ) ;

ite = ( vteprogte./ Z0te ) .*exp(-1i.*kz0.*z) ;

itm = ( vtmprogtm./ Z0tm ) .*exp(-1i.*kz0.*z) ;

vte = vteprogte .*exp(-1i.*kz0.*z) ;

vtm = vtmprogtm .*exp(-1i.*kz0.*z) ;

end