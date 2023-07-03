%% EE4620 Exam Project: Integrated Communications Network Using Surface Waves
% 27/06/2023
% Exam Date -> 03/07/2023
% Figure Defaults 
set(0,'DefaultLineLineWidth',2)
set(0,'defaultAxesFontSize',18)
set(0,'defaultAxesLinewidth',2)
set(0,'defaultfigureposition',[100 100 600 600])
%---------------------------------------------------------------------------

L = 3e-3 ;
W = 1e-3 ;
freq = 10e9 ;
clear all
h = 2.2e-3;
z = h+eps ;
er = 12 ;
freq = 10.*1e9 ;
freqn = linspace(1e9,15e9,1001);
lambdan=3e8./freqn ;
lambda = 3e8./freq ;
k0 = 2*pi./lambda ;
ky = 0 ;
zeta0 = 120*pi ;
no_ofpt = 1001 ;
h_lambdas = h.*sqrt(er).*freqn./3e8 ;

lst= k0.*sqrt(er) ;
krho = linspace(eps,lst,no_ofpt) ;
krho_norm = linspace(1,sqrt(er),no_ofpt) ;

% figure
% plot(krho./k0,abs(D))
% figure
% plot(krho./k0,abs(Dn))
% plot(freqn,abs(krholoop))
for ii = length(freqn):-1:1
    k0 = 2.*pi.* freqn(ii)./3e8 ;
    krho = krho_norm .* k0 ;
    [krhonte(ii),krho_newte(ii)] = IterativeMethod(h,1,er,freqn(ii),'GroundSlab','TE',krho) ;
    [krhontm(ii),krho_newtm(ii)] = IterativeMethod(h,1,er,freqn(ii),'GroundSlab','TM',krho) ;
end
figure
hold on
plot(freqn./1e9,abs(krhonte),'DisplayName','TE1') ;
plot(freqn./1e9,abs(krhontm),'DisplayName','TM0') ;
ylim([1 2.7]) 
title(sprintf('$\\left|k_{\\rho}^{norm}\\right|$ vs Frequency, $\\varepsilon_r=%.f$',er),'Interpreter','latex') ;
xlabel('Frequency[GHz]','Interpreter','latex');
ylabel('$\left|k_{\rho} ^{norm}\right|$','Interpreter','latex') ;
legend('Interpreter','latex','Location','best')
grid on;

figure
hold on
plot(h_lambdas,abs(krhonte),'DisplayName','TE1') ;
plot(h_lambdas,abs(krhontm),'DisplayName','TM0') ;
ylim([1 2.7]) 
xlim([0 0.38])
title(sprintf('$\\left|k_{\\rho}^{norm}\\right|$ vs $\\frac{h}{\\lambda_d}$, $\\varepsilon_r=%.f$',er),'Interpreter','latex') ;
xlabel('$\frac{h}{\lambda_d}$','Interpreter','latex');
ylabel('$\left|k_{\rho} ^{norm}\right|$','Interpreter','latex') ;
legend('Interpreter','latex','Location','best')
grid on;
