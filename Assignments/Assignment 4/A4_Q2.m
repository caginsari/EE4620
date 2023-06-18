%% Question 2 Leaky-wave propagation constants for Semi-inf superstrate
%% Change in  k_rho with er.
close all
clear

h = 10e-3 ;
freq = 15.*1e9 ;
er =linspace(1.049,25,1001);
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

    [klwTMn(ii),klwTM(ii)] = IterativeMethod(h,hs(ii),er(ii),freq,'SemiInf','TM',krho);
    klwTMn(isnan(klwTMn)) = 0 ;
    [klwTEn(ii),klwTE(ii)] = IterativeMethod(h,hs(ii),er(ii),freq,'SemiInf','TE',krho);
    klwTEn(isnan(klwTEn)) = 0 ;

end

k0 = 2*pi .*freq./3e8 ;
krhoTM = krhogSemiInf(k0, er, h, 'TM') ;
krhoTE = krhogSemiInf(k0, er, h, 'TE') ;
ks = k0.*sqrt(er) ;

figure;
hold on
plot(er,real(klwTMn),'k','DisplayName','Real,Numeric,TM') ;
plot(er,imag(klwTMn),'k--','DisplayName','Imag,Numeric,TM') ;
plot(er,real(krhoTM)./ks,'b','DisplayName','Real,Approx,TM') ;
plot(er,imag(krhoTM)./ks,'b--','DisplayName','Imag,Approx,TM') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_d$','Interpreter','latex')
xlabel('$\varepsilon_r$','Interpreter','latex')
title('Semi Infinite Superstrate','Interpreter','latex');
grid on;grid minor;
ylim([-0.5 0.5])
hold off

figure;
hold on
plot(er,real(klwTEn),'k','DisplayName','Real,Numeric,TE')  ;
plot(er,imag(klwTEn),'k--','DisplayName','Imag,Numeric,TE')  ;
plot(er,real(krhoTE)./ks,'b','DisplayName','Real,Approx,TE') ;
plot(er,imag(krhoTE)./ks,'b--','DisplayName','Imag,Approx,TE') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_d$','Interpreter','latex')
xlabel('$\varepsilon_r$','Interpreter','latex')
title('Semi Infinite Superstrate','Interpreter','latex');
grid on;grid minor;
ylim([-0.5 0.5])
hold off

%% Question 2 Change in krho with frequency at er=  12
% Constants
clear
h = 10e-3 ;
er = 12 ;
freq = 11.*1e9 ;
freqn = linspace(13e9,17e9,1001);
k0 = 2*pi .*freq./3e8 ;
lambda = 3e8./freq ;
w = lambda./20 ;
lst= k0 ;
no_ofpt = 1001 ;
krho = linspace(eps,lst,no_ofpt) ;
krho_norm = linspace(1,sqrt(er),no_ofpt) ;

klwTM = zeros(size(freqn)) ;
klwTE = zeros(size(freqn)) ;
for ii =length(freqn):-1:1

    [klwTMn(ii),klwTM(ii)] = IterativeMethod(h,0,er,freqn(ii),'SemiInf','TM',krho);
    klwTMn(isnan(klwTMn)) = 0 ;
    [klwTEn(ii),klwTE(ii)] = IterativeMethod(h,0,er,freqn(ii),'SemiInf','TE',krho);
    klwTEn(isnan(klwTEn)) = 0 ;

end

k0 = 2*pi .*freqn./3e8 ;
ks = k0.*sqrt(er) ;
krhoTM = krhogSemiInf(k0, er, h, 'TM') ;
krhoTE = krhogSemiInf(k0, er, h, 'TE') ;

figure;
hold on
plot(freqn./1e9,real(klwTMn),'k','DisplayName','Real,Numeric,TM') ;
plot(freqn./1e9,imag(klwTMn),'k--','DisplayName','Imag,Numeric,TM') ;
plot(freqn./1e9,real(krhoTM)./ks,'b','DisplayName','Real,Approx,TM') ;
plot(freqn./1e9,imag(krhoTM)./ks,'b--','DisplayName','Imag,Approx,TM') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_d$','Interpreter','latex')
xlabel('Frequency[GHz]','Interpreter','latex')
title('Semi Infinite Superstrate','Interpreter','latex');
grid on;grid minor;
ylim([-0.2 0.15])
hold off

figure;
hold on
plot(freqn./1e9,real(klwTEn),'k','DisplayName','Real,Numeric,TE')  ;
plot(freqn./1e9,imag(klwTEn),'k--','DisplayName','Imag,Numeric,TE')  ;
plot(freqn./1e9,real(krhoTE)./ks,'b','DisplayName','Real,Approx,TE') ;
plot(freqn./1e9,imag(krhoTE)./ks,'b--','DisplayName','Imag,Approx,TE') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_d$','Interpreter','latex')
xlabel('Frequency[GHz]','Interpreter','latex')
grid on;grid minor;
ylim([-0.2 0.15])
title('Semi Infinite Superstrate','Interpreter','latex');
hold off

%% Bandwidth and Directivity
clc
clear all

% FF parameters
freq = linspace(10e9,20e9,103); 
R_FF = 1;
phi  = (eps:2:360) * pi / 180;
h = 10e-3 ;
lambda = 3e8 ./ freq;
k0 = 2 .* pi ./ lambda;

theta = linspace(eps, 89.9, 303) * pi / 180;
dth = theta(2) - theta(1);
dph = phi(2) - phi(1);
[TH, PH] = meshgrid(theta, phi);
zeta0 = 120*pi ;
er = 1:2:25 ;
W = lambda ./ 20 ; 
L = lambda./2 ;


for jj = 1:length(er)
    for ff = 1:length(freq)   
        ks = k0(ff).*sqrt(er(jj)) ;
        [KX,KY,KZ] =PropVectors(ks,TH,PH) ;
        KRHO = sqrt(KX.^2 + KY.^2);
        Z    = R_FF * cos(TH);
        hs = lambda(ff)./(4.*sqrt(er(jj))) ;
         
        [vte,vtm,~,itm,~] = trxline_semi_inf_Superstrate(k0(ff), er(jj), h, KRHO ,Z, 'Layer2') ;
        % equivalent dielectric constant
        keq = (ks+k0(ff))./2 ;


        % calculate Green's function
        [em_sgf] = SpectralGFem(k0(ff),ks,er(jj),KX,KY,vtm,vte,itm,'Layer1',zeta0,KRHO) ;
        Gxx = em_sgf(:,:,1,1) ;
        Gyx = em_sgf(:,:,2,1) ;
        Gzx = em_sgf(:,:,3,1) ;
        % calculate FT of current distribution
        Mx = FTCurrent( keq, KX, KY, L(ff), W(ff) ) ;
        % calculate far field
        [Eth, Eph] = farfield( TH, PH, KZ, Gxx, Gyx, Gzx, Mx, Z,R_FF,0,ks) ;
        Etot(:,:,ff) = sqrt( abs(Eth).^2 +abs(Eph).^2 ) ;
        Etot(isnan(Etot(:,:,ff) ) ) = 0 ;
        [Dir, prad(ff)] = Direc(Etot(:,:,ff), TH, dth, dph, R_FF);
        D(ff,jj) = Dir(1,1) ;

    end
end
figure 
plot(freq./1e9, 10*log10(abs(D))) ;
xlabel('$Freq[GHz]$','Interpreter','latex');
ylabel('Directivity','Interpreter','latex');
title('Directivity With Incresing Frequency')
title('Semi-Infinite Superstrate','Interpreter','latex');
grid on; grid minor ;
DdB = 10.*log10(D);
%%
% Calculate Bandwidth
fH=zeros(size(er)) ;
fL = zeros(size(er)) ; 
for ii = 1:length(er)
    [r,c]=findpeaks(abs(DdB(:,ii))) ;
    rf = 0;
    if length(r)==1
        tolerance = 0.00005;
        while length(rf)<=3           
            [rf,cf] =find(abs(r-3-DdB(:,ii))<tolerance) ;
            tolerance=tolerance+0.001;
            if tolerance>0.8
                rf = [1 1 1 1] ;
                warning('higher than tolerance')
            end
        end
        fH(ii) = freq(max(rf)) ;
        fL(ii) = freq(min(rf)) ;
    end
end

BW = 200 .* (fH-fL)./(fH+fL) ;
BW(isnan(BW))= 0;
figure;
plot(er,BW);
xlabel('$\varepsilon_r$','Interpreter','latex')
ylabel('BW[%]','Interpreter','latex')