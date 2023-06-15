clear 
close all
%% EE4620 Assignment 4
% Leaky wave antenna characterisation

%% Question 1: Leaky Wave propagation constants
% Constants
h = 15e-3 ;
hs = 2.1e-3 ;
er = 12 ;
freq = 11.*1e9 ;
freqn = linspace(9e9,11e9,1001);
k0 = 2*pi .*freq./3e8 ;
lambda = 3e8./freq ;
w = lambda./20 ;
lst= k0 ;
no_ofpt = 1001 ;
krho = linspace(eps,lst,no_ofpt) ;
krho_norm = linspace(1,sqrt(er),no_ofpt) ;

klwTM = zeros(size(freq)) ;
klwTE = zeros(size(freq)) ;
for ii =length(freqn):-1:1

    [klwTMn(ii),klwTM(ii)] = IterativeMethod(h,hs,er,freqn(ii),'SuperStrate','TM',krho);
    klwTMn(isnan(klwTMn)) = 0 ;
    [klwTEn(ii),klwTE(ii)] = IterativeMethod(h,hs,er,freqn(ii),'SuperStrate','TE',krho);
    klwTEn(isnan(klwTEn)) = 0 ;

end

k0 = 2*pi .*freqn./3e8 ;
krhoTM = krhogSuperStrate(k0, er, h, 'TM') ;
krhoTE = krhogSuperStrate(k0, er, h, 'TE') ;

figure;
hold on
plot(freqn./1e9,real(klwTMn),'k','DisplayName','Real,Numeric,TM') ;
plot(freqn./1e9,imag(klwTMn),'k--','DisplayName','Imag,Numeric,TM') ;
plot(freqn./1e9,real(krhoTM)./k0,'b','DisplayName','Real,Approx,TM') ;
plot(freqn./1e9,imag(krhoTM)./k0,'b--','DisplayName','Imag,Approx,TM') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_0$','Interpreter','latex')
xlabel('Frequency[GHz]','Interpreter','latex')
grid on;grid minor;
ylim([-0.53 0.44])
hold off

figure;
hold on
plot(freqn./1e9,real(klwTEn),'k','DisplayName','Real,Numeric,TE')  ;
plot(freqn./1e9,imag(klwTEn),'k--','DisplayName','Imag,Numeric,TE')  ;
plot(freqn./1e9,real(krhoTE)./k0,'b','DisplayName','Real,Approx,TE') ;
plot(freqn./1e9,imag(krhoTE)./k0,'b--','DisplayName','Imag,Approx,TE') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_0$','Interpreter','latex')
xlabel('Frequency[GHz]','Interpreter','latex')
grid on;grid minor;
ylim([-0.53 0.44])
hold off


%% Change in  k_rho with er.
clear

h = 15e-3 ;
freq = 10.*1e9 ;
er =linspace(1-2*eps,25,1001);
k0 = 2*pi .*freq./3e8 ;
lambda = 3e8./freq ;
hs = lambda./(4 .*sqrt(er) ) ;
w = lambda./20 ;
lst= k0 ;
no_ofpt = 1001 ;
krho = linspace(eps,lst,no_ofpt) ;
% krho_norm = linspace(1,sqrt(er),no_ofpt) ;

klwTM = zeros(size(freq)) ;
klwTE = zeros(size(freq)) ;
for ii =1:length(er)

    [klwTMn(ii),klwTM(ii)] = IterativeMethod(h,hs(ii),er(ii),freq,'SuperStrate','TM',krho);
    klwTMn(isnan(klwTMn)) = 0 ;
    [klwTEn(ii),klwTE(ii)] = IterativeMethod(h,hs(ii),er(ii),freq,'SuperStrate','TE',krho);
    klwTEn(isnan(klwTEn)) = 0 ;

end

k0 = 2*pi .*freq./3e8 ;
krhoTM = krhogSuperStrate(k0, er, h, 'TM') ;
krhoTE = krhogSuperStrate(k0, er, h, 'TE') ;

figure;
hold on
plot(er,real(klwTMn),'k','DisplayName','Real,Numeric,TM') ;
plot(er,imag(klwTMn),'k--','DisplayName','Imag,Numeric,TM') ;
plot(er,real(krhoTM)./k0,'b','DisplayName','Real,Approx,TM') ;
plot(er,imag(krhoTM)./k0,'b--','DisplayName','Imag,Approx,TM') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_0$','Interpreter','latex')
xlabel('$\varepsilon_r$','Interpreter','latex')
grid on;grid minor;
ylim([-0.5 0.5])
hold off

figure;
hold on
plot(er,real(klwTEn),'k','DisplayName','Real,Numeric,TE')  ;
plot(er,imag(klwTEn),'k--','DisplayName','Imag,Numeric,TE')  ;
plot(er,real(krhoTE)./k0,'b','DisplayName','Real,Approx,TE') ;
plot(er,imag(krhoTE)./k0,'b--','DisplayName','Imag,Approx,TE') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_0$','Interpreter','latex')
xlabel('$\varepsilon_r$','Interpreter','latex')
grid on;grid minor;
ylim([-0.5 0.5])
hold off

%% Question 1 Bandwidth

% Constant
clc
clear all
close all

% FF parameters
freq = 28e9 ;
R_FF = 1;
phi  = (eps:2:360) * pi / 180;
h = 5.4e-3 ;
lambda = 3e8 / freq;
k0 = 2 * pi / lambda;

theta = linspace(eps, 89.9, 303) * pi / 180;
dth = theta(2) - theta(1);
dph = phi(2) - phi(1);
[TH, PH] = meshgrid(theta, phi);
zeta0 = 120*pi ;
er = linspace(1,25,100) ;
W = lambda ./ 20 ; 
L = lambda./2 ;

[KX,KY,KZ] =PropVectors(k0,TH,PH) ;
KRHO = sqrt(KX.^2 + KY.^2);
Z    = R_FF * cos(TH);

for ff = 1:length(er) 

    
    hs = lambda./(4 .*sqrt(er(ff) ) ) ;
    k1 = k0.*sqrt(er(ff)) ;
    keq = (k0+k0)./2 ;

    [vtm, vte, itm, ~, ks,kz0] = trxline_Superstrate(k0, zeta0, er(ff), h, hs, KRHO, 'Layer3' ,Z, freq) ;
    % calculate Green's function
    [em_sgf] = SpectralGFem(k0,ks,er(ff),KX,KY,vtm,vte,itm,'Layer2',zeta0,KRHO) ;
    Gxx = em_sgf(:,:,1,1) ;
    Gyx = em_sgf(:,:,2,1) ;
    Gzx = em_sgf(:,:,3,1) ;
    % calculate FT of current distribution
    Mx = FTCurrent( keq, KX, KY, L, W ) ;
    % calculate far field
    [Eth, Eph] = farfield( TH, PH, KZ, Gxx, Gyx, Gzx, Mx, Z,R_FF,h,k0) ;
    Etot(:,:,ff) = sqrt( abs(Eth).^2 +abs(Eph).^2 ) ;
    Etot(isnan(Etot(:,:,ff) ) ) = 0 ;
    [Dir, prad(ff)] = Direc(Etot(:,:,ff), TH, dth, dph, R_FF);
    D(ff) = Dir(1,1) ;

end
figure 
plot(er, 10*log10(abs(D))) ;
xlabel('$\varepsilon_r$','Interpreter','latex');
ylabel('Directivity','Interpreter','latex');
title('Directivity with incresing relative permittivity')
grid on; grid minor ;
