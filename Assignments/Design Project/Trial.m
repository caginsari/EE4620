theta = linspace(eps,pi,101) ;
f = 10e9 ;
k0 = 2*pi *(f/3e8) ;
h = 2./k0 .* cos(pi./2 .* cos(theta) )./sin(theta) ;

figure
plot(rad2deg(theta),h)