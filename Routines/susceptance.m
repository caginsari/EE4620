function [B] = susceptance(dy,w,m,f)

m = [-m:-1 1:m ]  ;

omg = 2.*pi.*f ; 
cons = omg.*8.85e-12.*dy./pi ;
arg = abs(sinc(pi.*m.*w./(dy.*pi))).^2./abs(m) ;
B =cons.* sum(arg) ;
end