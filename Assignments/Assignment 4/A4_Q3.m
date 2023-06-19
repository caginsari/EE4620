%% Question 3 Change in krho with frequency at er=12
% Constants
clear
h = 10e-3 ;
er = 10 ;
freq = 11.*1e9 ;
freqn = linspace(13e9,17e9,202);
k0 = 2*pi .*freq./3e8 ;
lambda = 3e8./freq ;
w = lambda./20 ;
lst= k0 ;
no_ofpt = 202 ;
krho = linspace(eps,lst,no_ofpt) ;
krho_norm = linspace(1,sqrt(er),no_ofpt) ;

klwTM = zeros(size(freqn)) ;
klwTE = zeros(size(freqn)) ;
for ii =length(freqn):-1:1

    [klwTMn(ii),klwTM(ii)] = IterativeMethod(h,0,er,freqn(ii),'SemiInf','TM0',krho);
    klwTMn(isnan(klwTMn)) = 0 ;

end

k0 = 2*pi .*freqn./3e8 ;
ks = k0.*sqrt(er) ;
krhoTM = krhogSemiInf(k0, er, h, 'TM0') ;

figure;
hold on
plot(freqn./1e9,real(klwTMn),'k','DisplayName','Real,Numeric,TM0') ;
plot(freqn./1e9,imag(klwTMn),'k--','DisplayName','Imag,Numeric,TM0') ;
plot(freqn./1e9,real(krhoTM)./ks,'b','DisplayName','Real,Approx,TM0') ;
plot(freqn./1e9,imag(krhoTM)./ks,'b--','DisplayName','Imag,Approx,TM0') ;
legend('Location','best','Interpreter','latex')
ylabel('$k_{\rho}/k_d$','Interpreter','latex')
xlabel('Frequency[GHz]','Interpreter','latex')
title('Semi Infinite Superstrate','Interpreter','latex');
grid on;grid minor;
ylim([-0.05 0.33])
hold off

d = pi ./ real(klwTM(101)) ;

% Radiation from double slot antenna and single slot antenna
%% Bandwidth and Directivity
clc
clear all

% FF parameters
freq = linspace(10e9,20e9,203);
R_FF = 1;
phi  = (eps:1:360) * pi / 180;
h = 10e-3 ;
lambda = 3e8 ./ freq;
k0 = 2 .* pi ./ lambda;

theta = linspace(eps, 89.9, 303) * pi / 180;
dth = theta(2) - theta(1);
dph = phi(2) - phi(1);
[TH, PH] = meshgrid(theta, phi);
zeta0 = 120*pi ;
er = 10 ;
W = lambda ./ 20 ; 
L = lambda./2 ;

for jj = 1:length(er)
    for ff = 1:length(freq)   
        ks = k0(ff).*sqrt(er(jj)) ;
        [KX,KY,KZ] =PropVectors(ks,TH,PH) ;
        KRHO = sqrt(KX.^2 + KY.^2);
        Z    = R_FF * cos(TH);
        k_comp(:, :, :, 1) = KX;
        k_comp(:, :, :, 2) = KY;
        k_comp(:, :, :, 3) = KZ;
         
        [vte,vtm,~,itm,~] = trxline_semi_inf_Superstrate(k0(ff), er(jj), h, KRHO ,Z, 'Layer2') ;
        % equivalent dielectric constant
        keq = (k0(ff)+k0(ff))./2 ;

        % calculate Green's function
        [em_sgf] = SpectralGFem(k0(ff),ks,er(jj),KX,KY,vtm,vte,itm,'Layer1',zeta0,KRHO) ;
        Gxy = em_sgf(:,:,1,2) ;
        Gyy = em_sgf(:,:,2,2) ;
        Gzy = em_sgf(:,:,3,2) ;

        % calculate FT of current distribution
        Mx = FTCurrentDoubleSlot(KX,KY,keq,W(ff),L(ff),0.5.*lambda(ff),'single') ;

        % calculate far field
        [Eth, Eph] = farfield( TH, PH, KZ, Gxy, Gyy, Gzy, Mx, Z,R_FF,0,ks) ;
        Etotsingle(:,:,ff) = sqrt( abs(Eth).^2 +abs(Eph).^2 ) ;
        Etotsingle(isnan(Etotsingle(:,:,ff) ) ) = 0 ;
        [Dir, prad(ff)] = Direc(Etotsingle(:,:,ff), TH, dth, dph, R_FF);
        Ddouble(ff,jj) = Dir(1,1) ;

        % calculate FT of current distribution
        Mx = FTCurrentDoubleSlot(KX,KY,keq,W(ff),L(ff),0.5.*lambda(ff),'double') ;
        % calculate far field
        [Eth, Eph] = farfield( TH, PH, KZ, Gxy, Gyy, Gzy, Mx, Z,R_FF,0,ks) ;
        Etotdouble(:,:,ff) = sqrt( abs(Eth).^2 +abs(Eph).^2 ) ;
        Etotdouble(isnan(Etotdouble(:,:,ff) ) ) = 0 ;
        [Dir, prad(ff)] = Direc(Etotdouble(:,:,ff), TH, dth, dph, R_FF);
        Dsingle(ff,jj) = Dir(1,1) ;

    end
end

figure 
plot(freq./1e9, 10*log10(abs(Dsingle))) ;
xlabel('$Freq[GHz]$','Interpreter','latex');
ylabel('Directivity','Interpreter','latex');
title('Directivity With Incresing Frequency')
title('Semi-Infinite Superstrate','Interpreter','latex');
grid on; grid minor ;

figure 
plot(freq./1e9, 10*log10(abs(Ddouble))) ;
xlabel('$Freq[GHz]$','Interpreter','latex');
ylabel('Directivity','Interpreter','latex');
title('Directivity With Incresing Frequency')
title('Semi-Infinite Superstrate','Interpreter','latex');
grid on; grid minor ;

figure
hold on
plot(rad2deg(theta),20.*log10(abs(Etotsingle(1,:,102) ) ) -20.*log10(max(max(abs(Etotsingle(:,:,102) )))),'DisplayName','phi = 0')
plot(rad2deg(theta),20.*log10(abs(Etotsingle(46,:,102) ) ) -20.*log10(max(max(abs(Etotsingle(:,:,102) )))),'DisplayName','phi = 45')
plot(rad2deg(theta),20.*log10(abs(Etotsingle(91,:,102) ) ) -20.*log10(max(max(abs(Etotsingle(:,:,102) )))),'DisplayName','phi = 90')
hold off
xlim([0 35])
ylim([-40 0])
title('Single')
legend

figure
hold on
plot(rad2deg(theta),20.*log10(abs(Etotdouble(1,:,102) ) ) -20.*log10(max(max(abs(Etotdouble(:,:,102) )))),'DisplayName','phi = 0')
plot(rad2deg(theta),20.*log10(abs(Etotdouble(46,:,102) ) ) -20.*log10(max(max(abs(Etotdouble(:,:,102) )))),'DisplayName','phi = 45')
plot(rad2deg(theta),20.*log10(abs(Etotdouble(91,:,102) ) ) -20.*log10(max(max(abs(Etotdouble(:,:,102) )))),'DisplayName','phi = 90')
hold off
xlim([0 35])
ylim([-40 0])
title('Double')
legend

figure
hold on
plot(rad2deg(theta),20.*log10(abs(Etotsingle(1,:,102) ) ) -20.*log10(max(max(abs(Etotsingle(:,:,102) )))),'DisplayName','phi = 0')
plot(rad2deg(theta),20.*log10(abs(Etotsingle(46,:,102) ) ) -20.*log10(max(max(abs(Etotsingle(:,:,102) )))),'DisplayName','phi = 45')
plot(rad2deg(theta),20.*log10(abs(Etotsingle(91,:,102) ) ) -20.*log10(max(max(abs(Etotsingle(:,:,102) )))),'DisplayName','phi = 90')
hold off
xlim([0 35])
ylim([-40 0])
title('Single')
legend

