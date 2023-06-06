function [Eff] = delt_source(f,TH,PHI,k0,r) 
omega = 2.* pi .*f ;

% Field from elementary dipole
Eth = -1i.* omega .* pi.*4e-7 .* cos(TH) .*cos(PHI) .* exp(-1i .* k0 .* r);
Eph = 1i.* omega .* pi.*4e-7 .*sin(PHI) .* exp(-1i .* k0 .* r) ;
Etot = sqrt(abs(Eth).^2 + abs(Eph).^2) ;

Eff(:,:,1) = Eth ;
Eff(:,:,2) = Eph ;
Eff(:,:,3) = Etot ;

end