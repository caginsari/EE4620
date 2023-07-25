% 01/07/2023
% Exam Date -> 03/07/2023
% Figure Defaults 
set(0,'DefaultLineLineWidth',2)
set(0,'defaultAxesFontSize',18)
set(0,'defaultAxesLinewidth',2)
set(0,'defaultfigureposition',[100 100 600 600])
clear 
l = 3e-3 ;
w = 1e-3 ;
h = 2.2e-3 ;
er = 12; 
f= 10e9 ;
k0 = 2.*pi.* f./3e8 ;
no_ofpt = 1001 ;
phi = linspace(-pi/2,pi/2,303) ;
lambdas = 3e8./(f.*sqrt(er)) ;
rho = linspace(10,15,500)*1e-2; % linspace(eps,1,200) ; 
[RHO,PHI] = meshgrid(rho,phi) ;
z = h ;
Ra = 50 ;

lst= k0.*sqrt(er) ;
krho = linspace(eps,lst,no_ofpt) ;
% Erhotm=zeros(size(RHO))
Erhotm = zeros(length(phi),length(rho)) ;
Prx = zeros(length(phi),length(rho)) ;
h2 = zeros(length(phi),length(rho)) ;
for rr = 1:length(rho)
    for ii = 1:length(phi)
        [~,kswTM] = IterativeMethod(h,1,er,f,'GroundSlab','TM',krho) ;
        [VrTM,IrTM] = Residue_GroundSlab(k0,er,h,kswTM,z,f,'TM') ;
        [Erhotm(ii,rr),~,~,h2(ii,rr)] = SwFields(k0,kswTM,er,VrTM,IrTM,rho(rr),phi(ii),l,w) ;

        Voc = abs(Erhotm(ii,rr).*h2(ii,rr) ) ;
        Prx(ii,rr) = Voc.^2 ./ (8.*Ra) ;
    end
end
% [~,~,ResTM,~] = residue_stratified(k0, kswTE, kswTM, z, 'GroundSlab',h,er)

figure 
plot(rad2deg(phi),10.*log10(Prx(:,1)) ) ;
title('Power Received With Changing Incident Angle','Interpreter','latex') ;
xlim([-90 90]);ylim([-60 0]);
ylabel('$P_{rx}^{norm}[dB]$','Interpreter','latex') ;
xlabel('$\phi_i$[Degree]','Interpreter','latex') ;
grid on;

figure 
plot(rad2deg(phi),10.*log10(Prx(:,1)) - 10.*log10(max(Prx(:,1)) ) ) ;
title('Power Received With Changing Incident Angle','Interpreter','latex') ;
xlim([-90 90]);ylim([-60 0]);
ylabel('$P_{rx}^{norm}[dB]$','Interpreter','latex') ;
xlabel('$\phi_i$[Degree]','Interpreter','latex') ;
grid on;

figure 
hold on
plot(rho/1e-2, 10.*log10( Prx(round( length(phi)/2 ), :) )  ,'k','DisplayName','Prx[dB]' );
% plot(rho./1e-2, 10.*log10(1./rho)-10.*log10(max(1./rho)),'r--',DisplayName='$1/\rho$[dB]')
% plot(rho./1e-2,10.*log10(1./sqrt(rho))-10.*log10(max(1./sqrt(rho)) ) ,'r--', DisplayName='1/sqrt(rho)')
title('Power Received With Radial Distance','Interpreter','latex')
xlabel('$\rho[cm]$','Interpreter','latex')
ylabel('$P_{rx}^{norm}[dB]$','Interpreter','latex')
legend('Interpreter','latex','Location','best')
grid on;
