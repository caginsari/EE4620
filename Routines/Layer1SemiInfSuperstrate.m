function [vte,vtm,ite,itm] = Layer1SemiInfSuperstrate(gamma1te,gamma1tm,kz0,h,z,Z0te,Z0tm)

vte = exp(-1i.*kz0.*z) ./(1+gamma1te.*exp(-2.*1i.*kz0.*h)) .* ...
    (1+gamma1te.*exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z));

vtm = exp(-1i.*kz0.*z) ./(1+gamma1tm.*exp(-2.*1i.*kz0.*h)) .* ...
    (1+gamma1tm.*exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z));

ite = exp(-1i.*kz0.*z) ./(Z0te.*(1+gamma1te.*exp(-2.*1i.*kz0.*h) ) ) .* ...
    (1-gamma1te.*exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z));

itm = exp(-1i.*kz0.*z) ./(Z0tm .* (1+gamma1tm.*exp(-2.*1i.*kz0.*h) ) ) .* ...
    (1-gamma1tm.*exp(-2.*1i.*kz0.*h).*exp(2.*1i.*kz0.*z));

end