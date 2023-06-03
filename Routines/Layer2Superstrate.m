function [vte,vtm,ite,itm]=Layer2Superstrate(kz0,kzs,h,hs,Zste,Zstm,gamma1te,gamma1tm,gamma2te,gamma2tm)
% EE4620 Assignment 1: Layer 2 of the superstrate

vte =( exp(1i.*kz0.*h) .* exp(1i.*kzs.*h) .* (1+gamma1te) ) ./ ( (gamma1te + exp(2.*1i.*kz0.*h)) .*(1+gamma2te.* exp(-2.*1i.*kzs.*hs) ) )...
    .*( exp(-1i.*kzs.*z) + exp(1i.*kzs.*z).*gamma2te .* exp(-2.*1i.*kzs.*(h+hs) ) )  ;

vtm =( exp(1i.*kz0.*h) .* exp(1i.*kzs.*h) .* (1+gamma1tm) ) ./ ( (gamma1tm + exp(2.*1i.*kz0.*h)) .*(1+gamma2tm.* exp(-2.*1i.*kzs.*hs) ) )...
    .*( exp(-1i.*kzs.*z) + exp(1i.*kzs.*z).*gamma2tm .* exp(-2.*1i.*kzs.*(h+hs) ) )  ;

ite =( exp(1i.*kz0.*h) .* exp(1i.*kzs.*h) .* (1+gamma1te) ) ./ ( Zste .* (gamma1te - exp(2.*1i.*kz0.*h)) .*(1+gamma2te.* exp(-2.*1i.*kzs.*hs) ) )...
    .*( exp(-1i.*kzs.*z) + exp(1i.*kzs.*z).*gamma2te .* exp(-2.*1i.*kzs.*(h+hs) ) )  ;

itm =( exp(1i.*kz0.*h) .* exp(1i.*kzs.*h) .* (1+gamma1tm) ) ./ ( Zstm .* (gamma1tm - exp(2.*1i.*kz0.*h)) .*(1+gamma2tm.* exp(-2.*1i.*kzs.*hs) ) )...
    .*( exp(-1i.*kzs.*z) + exp(1i.*kzs.*z).*gamma2tm .* exp(-2.*1i.*kzs.*(h+hs) ) )  ;

end 
