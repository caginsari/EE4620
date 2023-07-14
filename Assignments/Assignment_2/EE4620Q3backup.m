% Constant
clc
clear all
close all

% FF parameters
R_FF = 1;
phi  = (eps:2:360) * pi / 180;

theta = linspace(-90, 90, 303) * pi / 180;
dth = theta(2) - theta(1);
dph = phi(2) - phi(1);
[TH, PH] = meshgrid(theta, phi);

% calculate voltages and currents
% txline_groundslab(k0, er, h, KRHO, Z);
zeta0 = 120*pi ;
% [vtm, vte, itm, ite,kz1] = txline_GroundSlab(k0,zeta0,er,Z,h,KRHO,f) ;
%  figure

h = 5.4e-3 ;
er = 12 ;
freq = 28e9 ;


lambda = 3e8./freq ;
k0 = 2.*pi./lambda ;
W = lambda ./ 20 ;
L = lambda./2 ;
k1 = k0./sqrt(er) ;
keq = k1;

[KX,KY,KZ] = PropVectors(k1,TH,PH) ;
KRHO = sqrt(KX.^2 + KY.^2) ;
Z    = R_FF * cos(TH) ;

% calculate voltages and currents
[vtm, vte, itm, ite, ks] = trxline_semi_inf_Superstrate(k0,zeta0,er,h,KRHO,Z,freq,'Layer2') ;
% calculate Green's function
[em_sgf] = SpectralGFem(k0,ks,er,KX,KY,vtm,vte,itm,'Layer1',zeta0,KRHO) ;
Gxx = em_sgf(:,:,1,1) ;
Gyx = em_sgf(:,:,2,1) ;
Gzx = em_sgf(:,:,3,1) ;
% calculate FT of current distribution
Mx = FTCurrent( keq, KX, KY, L, W ) ;
% calculate far field
[Eth, Eph] = farfield( TH, PH, KZ, Gxx, Gyx, Gzx, Mx, Z,R_FF,h,k1) ;
Etot = sqrt( abs(Eth).^2 +abs(Eph).^2 ) ;

figure
hold on
plot(rad2deg(theta), 20*log10(abs(Etot(1,:,1)))-20*log10( abs(max(max(max(Etot(1,:,1) ) ) ) ) ),'DisplayName',sprintf('f=%.1f GHz, $\\phi=0^\\circ$',freq./1e9 ) ) ;
plot(rad2deg(theta), 20*log10(abs(Etot(round(length(phi)./4),:,1)))-20*log10(abs(max(max(max(Etot(1,:,1) ) ) ) ) ),'DisplayName',sprintf('f=%.1f GHz, $\\phi=90^\\circ$',freq./1e9 ) ) ;
hold off
legend('Interpreter','latex');
ylim([-40 0]);
title(['Resonant Leaky Wave Antenna Far-Field'],'Interpreter','latex')
xlabel('$\theta$[deg]','Interpreter','latex');
ylabel('$\left|E\right|[dB]$ ','Interpreter','latex') ;
legend('Location','best') ; grid on ; grid minor; 
hold off
% calculate voltages and currents
% [vtm, vte, itm, ite, ks] = trxline_semi_inf_Superstrate(k_0,zeta0,er(ii),h,krho,z,freq) ;