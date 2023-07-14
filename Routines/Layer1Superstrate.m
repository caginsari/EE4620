function [vte,vtm,ite,itm] = Layer1Superstrate(gamma1te,gamma1tm,h,z,kz0,Z0te,Z0tm)
% first layer of the superstrate

vte = (exp(2.*1i.*kz0.*h) ./ (gamma1te + exp(2.*1i.*kz0.*h) ) ).*exp(-1i.*kz0.*z).* ...
    ( 1 + gamma1te.* exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z) ) ; 

vtm = (exp(2.*1i.*kz0.*h) ./ (gamma1tm + exp(2.*1i.*kz0.*h) ) ).*exp(-1i.*kz0.*z).* ...
    ( 1 + gamma1tm .* exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z) ) ; 

ite = (exp(2.*1i.*kz0.*h) ./ ( Z0te .* (gamma1te - exp(2.*1i.*kz0.*h) ) ) ).*exp(-1i.*kz0.*z).* ...
    ( 1 - gamma1te.* exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z) ) ; 

itm = (exp(2.*1i.*kz0.*h) ./ ( Z0tm .* (gamma1tm - exp(2.*1i.*kz0.*h) ) ) ).*exp(-1i.*kz0.*z).* ...
    ( 1 - gamma1tm.* exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z) ) ; 

end