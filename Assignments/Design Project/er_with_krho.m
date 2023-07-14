%% EE4620 Design Project: Effect of er on the krho
% 30/06/2023
% Exam Date -> 03/07/2023
% Figure Defaults 
close 
set(0,'DefaultLineLineWidth',2)
set(0,'defaultAxesFontSize',18)
set(0,'defaultAxesLinewidth',2)
set(0,'defaultfigureposition',[100 100 600 600])
%---------------------------------------------------------------------------

L = 3e-3 ;
W = 1e-3 ;
freq = 25e9 ;
clear all
h = 2.5e-3;
z = h+eps ;
er = 12 ;
freq = 10.*1e9 ;
er = linspace(1,25,1001);
lambda = 3e8./freq ;
k0 = 2*pi./lambda ;
ky = 0 ;
zeta0 = 120*pi ;
no_ofpt = 1001 ;

lst= k0.*sqrt(25) ;
krho = linspace(eps,lst,no_ofpt) ;
krho_norm = linspace(1,sqrt(25),no_ofpt) ;

% figure
% plot(krho./k0,abs(D))
% figure
% plot(krho./k0,abs(Dn))
% plot(freqn,abs(krholoop))
for ii = length(er):-1:1
    k0 = 2.*pi.* freq./3e8 ;
    krho = krho_norm .* k0 ;
    [krhonte(ii),krho_newte(ii)] = IterativeMethod(h,1,er(ii),freq,'GroundSlab','TE',krho) ;
    [krhontm(ii),krho_newtm(ii)] = IterativeMethod(h,1,er(ii),freq,'GroundSlab','TM',krho) ;
end
figure
hold on
TE = plot(er,abs(krhonte),'DisplayName','TE1') ;
TM = plot(er,abs(krhontm),'DisplayName','TM0') ;
ylim([1 2.7]) 
title(sprintf('Propagation Constant vs Permittivity , $f=%.fGHz$',freq./1e9),'Interpreter','latex') ;
xlabel('$\varepsilon_r$','Interpreter','latex');
ylabel('$\left|k_\rho ^{norm}\right|$','Interpreter','latex') ;
legend('Interpreter','latex','Location','best')
grid on;

%% krho with increasing slab height
clear all

L = 3e-3 ;
W = 1e-3 ;
h = linspace(eps,3e-3,1001);
z = h+eps ;
freq = 10.*1e9 ;
er = 15 ;%linspace(1,25,1001);
lambda = 3e8./freq ;
k0 = 2*pi./lambda ;
ky = 0 ;
zeta0 = 120*pi ;
no_ofpt = 1001 ;

lst= k0.*sqrt(15) ;
krho = linspace(eps,lst,no_ofpt) ;
krho_norm = linspace(1,sqrt(15),no_ofpt) ;

% figure
% plot(krho./k0,abs(D))
% figure
% plot(krho./k0,abs(Dn))
% plot(freqn,abs(krholoop))
for ii = length(h):-1:1
    k0 = 2.*pi.* freq./3e8 ;
    krho = krho_norm .* k0 ;
    [krhonte(ii),krho_newte(ii)] = IterativeMethod(h(ii),1,er,freq,'GroundSlab','TE',krho) ;
    [krhontm(ii),krho_newtm(ii)] = IterativeMethod(h(ii),1,er,freq,'GroundSlab','TM',krho) ;
end

figure
hold on
TE = plot(h./1e-3,abs(krhonte),'DisplayName','TE1') ;
TM = plot(h./1e-3,abs(krhontm),'DisplayName','TM0') ;
ylim([1 2.7]) 
title(sprintf('Propagation Constant vs Slab Height, $f=%.fGHz,\\varepsilon=%.f$',freq./1e9,er),'Interpreter','latex') ;
xlabel('$h[mm]$','Interpreter','latex');
ylabel('$\left|k_\rho ^{norm}\right|$','Interpreter','latex') ;
legend('Interpreter','latex','Location','best')
grid on;
