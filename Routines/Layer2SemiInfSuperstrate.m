function [vte,vtm,ite,itm] = Layer2SemiInfSuperstrate(gamma1te,gamma1tm,kz0,h,z,Zste,Zstm,kzs)
% EE4620 Assignment 1: semiinf superstrate layer 2(slab)
% [vte,vtm,ite,itm] = Layer2SemiInfSuperstrate(gamma1te,gamma1tm,kz0,h,z,Zste,Zstm,kzs)

vte = ( exp(-1i.*kz0.*h) .*exp(1i.*kzs.*h)./(1+gamma1te.*exp(-2.*1i.*kz0.*h) ) ) .*(1+gamma1te).*exp(-1i.*kzs.*z) ;

vtm = ( exp(-1i.*kz0.*h) .*exp(1i.*kzs.*h)./(1+gamma1tm.*exp(-2.*1i.*kz0.*h) ) ) .*(1+gamma1tm).*exp(-1i.*kzs.*z) ;

ite = exp(-1i.*kz0.*h) .*exp(1i.*kzs.*h)./( Zste .*(1+gamma1te.*exp(-2.*1i.*kz0.*h) ) ) .*(1+gamma1te).*exp(-1i.*kzs.*z) ;

itm = exp(-1i.*kz0.*h) .*exp(1i.*kzs.*h)./(Zstm .*(1+gamma1tm.*exp(-2.*1i.*kz0.*h) ) ) .*(1+gamma1tm).*exp(-1i.*kzs.*z) ;

end