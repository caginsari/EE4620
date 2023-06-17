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

%% Question 2 Change in krho with frequency at er=  
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