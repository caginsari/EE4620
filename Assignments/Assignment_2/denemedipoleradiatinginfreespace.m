clear

numth=101;
numph= 303;
l = 0.25e-3 ;%15e-3 ;
w = 0.25e-3 ;%0.5e-3 ;
f= linspace(10,40,100)*1e9;%1.3748e+10 ;
k0 = 2.*pi.* f./3e8 ;
no_ofpt = 1001 ;
phi1 = linspace(eps,2*pi,303) ;
theta = linspace(eps,pi/2,numth);
phi = linspace(eps,2*pi,numph);
[TH,PHI] = meshgrid(theta,phi1);

% Erhotm=zeros(size(RHO))
for ii = 1:length(f)

    [kx,ky,kz] = PropVectors(k0(ii),TH,PHI) ;

    % calculate voltages and currents
    ej_SGF = EJ_SGF(1, k0(ii), kx, ky);
    Gxx = ej_SGF(:,:,1,1);
    Gyx = ej_SGF(:,:,2,1);
    Gzx = ej_SGF(:,:,3,1);

    dth = theta(2) - theta(1);
    dph = phi(2) - phi(1);

    Jx = FTCurrent(k0(ii), kx , ky , l , w);

    % calculate far field of the free space dipole
    [Eth, Eph] = farfield(TH, PHI, kz, Gxx, Gyx, Gzx, Jx,0,1,0,k0(ii));
    E_tot_element = sqrt(abs(Eth).^2 + abs(Eph).^2);
    [~, Pradiated_element(ii)] = Direc(E_tot_element, TH, dth, dph, 1); 


end

figure;
hold on 
PradP = plot(f./1e9,abs(prad)./abs(Pradiated_element) .*0.5 ,'DisplayName','Radiated Power');