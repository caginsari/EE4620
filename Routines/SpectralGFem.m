function [em_SGF] = SpectralGFem(k0,ks,er,kx,ky,vtm,vte,itm,Layer,Zeta0,krho)
% EE4620 Assignemnt 1 : Numerically calculates dyadic SGF. required input arguements are relative permittivity er, wavenumber k, propagation vector in x direction kx and y direction ky. 
%aplicable only to upper hemisphere

e0 = 8.85e-12 ;
mu = pi * 4e-7; 

Zetas = sqrt(mu./ (er .*e0) ) ;

% kz = -1i*sqrt(- (k.^2-kx.^2-ky.^2) ) ;

Dxx = (kx.*ky./krho.^2).*(vtm - vte) ;
Dxy = -1./krho.^2 .*(vtm.*kx.^2 +vte.*ky.^2) ;

Dyx = (vtm.*ky.^2 + vte.*kx.^2) ./ krho.^2 ;
Dyy = (kx.*ky./krho.^2).*(vte - vtm) ;

L1 = strcmpi('Layer1',Layer) ;
L2 = strcmpi('Layer2',Layer) ;

if L1 == 1
    % in substrate
    zeta = Zetas ;
    k = ks ;
elseif L2 == 1
    % in air
    zeta = Zeta0 ;
    k = k0 ;
else
    error(['\nInput must be one of two superstrate layers: \n Layer1 : dielectric \n Layer2 : air \n Your input: %s.'],Layer) ;
end

Dzx = -zeta .* ky./k .*itm ;
Dzy = zeta .* kx./k .*itm ;

%Gx
em_SGF(:,:,1,1) =  Dxx;
em_SGF(:,:,1,2) =  Dxy;

%Gy
em_SGF(:,:,2,1) =  Dyx;
em_SGF(:,:,2,2) =  Dyy;

%Gz
em_SGF(:,:,3,1) =  Dzx;
em_SGF(:,:,3,2) =  Dzy;

end