clear all
theta_lw =@(k0,betalw) asin(sqrt(12).*betalw./k0);

hs = 0.77e-3 ;
h = 5.4e-3 ;
er = 12 ;
freqn = [26,28,30] *1e9 ;
k0 = 2*pi .*freqn./3e8 ;
lambda = 3e8./freqn ;
w = lambda./20 ;
lst= k0(1) ;
no_ofpt = 1001 ;
krho = linspace(eps,lst,no_ofpt) ;
krho_norm = linspace(1,sqrt(er),no_ofpt) ;

klwTM = zeros(size(freqn)) ;
klwTE = zeros(size(freqn)) ;
for ii = 1:length(freqn)
    lst= k0(ii) ;
    no_ofpt = 1001 ;
    krho = linspace(eps,lst,no_ofpt) ;

    h = lambda(ii)/2 ;
    hs = lambda(ii)/(4*sqrt(er));
    [klwTMn(ii),klwTM(ii)] = IterativeMethod(h,hs,er,freqn(ii),'SuperStrate','TM',krho);
    klwTMn(isnan(klwTMn)) = 0 ;
    [klwTEn(ii),klwTE(ii)] = IterativeMethod(h,hs,er,freqn(ii),'SuperStrate','TE',krho);
    klwTEn(isnan(klwTEn)) = 0 ;

end

beta_tm = real(klwTM) ;
beta_te = real(klwTE) ;

theta_lw0_tm = rad2deg(theta_lw(k0,beta_tm)) 
theta_lw0_te = rad2deg(theta_lw(k0,beta_te)) 



