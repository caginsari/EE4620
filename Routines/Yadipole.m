function [Ya_dipole] = Yadipole(f,V0,m,dx,dy,TH,PHI)

k0 = 2.*pi.* 3e8 ./ f ;
mx= -m:m ;
my = -m:m; 
[MX,MY] = meshgrid(mx,my) ;
[kxm,kym] = FloquetModes(k0,TH,PHI,MX,MY,dx,dy) ;

Ya_dipole = 1./dx .* (-sinc(kxm(ii)s .*delt_d./(2.*pi) ) ) .* D_inf ;
[kx0,ky0,~] = PropVectors( k0, TH, PHI) ;

Spectral_ej = EJ_SGF(er,k0,kx0,kym) ;
Gxx = Spectral_ej(:,:,1,1) ;

D_inf = 1./dy .* Gxx .* besselj(0,kym.*wd./2) ;
Ya_dipole = 1./dx .* (-sinc(kxm(ii) .*delt_d./(2.*pi) ) ) .* D_inf ;
ix = Ya_dipole .*V0 .* exp(-1i .* kxm(ii).* x ) ;

for

end

end