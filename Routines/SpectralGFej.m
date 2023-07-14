function [ej_SGF] = SpectralGFej(k0,ks,er,kx,ky,vtm,vte,itm,ite,Zeta0,krho,z,h)
% EE4620 Assignemnt 1 : Numerically calculates dyadic SGF. required input arguements are relative permittivity er, wavenumber k, propagation vector in x direction kx and y direction ky. 
%aplicable only to upper hemisphere

e0 = 8.85e-12 ;
mu = pi * 4e-7; 

Zetas = sqrt(mu./ (er .*e0) ) ;

% kz = -1i*sqrt(- (k.^2-kx.^2-ky.^2) );

Dxx = -(vtm .* kx.^2 + vte.*ky.^2)./(krho.^2);
Dxy = 1./krho.^2 .*(vte -vtm) .* (kx.*ky) ;

Dyx = 1./krho.^2 .*(vte -vtm) .* (kx.*ky) ;
Dyy = -(vte.*kx.^2+vtm.*ky.^2)./krho.^2;

if z<h
    % in substrate
    zeta = Zetas ;
    k = ks ;
else
    % in air
    zeta = Zeta0 ;
    k = k0 ;
end

Dzx = zeta .* kx./k .*itm ;
Dzy = zeta .* ky./k .*itm ;

%Gx
ej_SGF(:,:,1,1) =  Dxx;
ej_SGF(:,:,1,2) =  Dxy;

%Gy
ej_SGF(:,:,2,1) =  Dyx;
ej_SGF(:,:,2,2) =  Dyy;

%Gz
ej_SGF(:,:,3,1) =  Dzx;
ej_SGF(:,:,3,2) =  Dzy;

end