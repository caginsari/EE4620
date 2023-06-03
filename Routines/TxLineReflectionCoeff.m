function [gammate,gammatm] = TxLineReflectionCoeff(Z1te,Z1tm,Z2te,Z2tm)
%EE4620 Assigment 1: reflection coefficient calculation of a transmission line 

Gamma =@(Z1,Z2) (Z1-Z2) ./ (Z1+Z2) ;

gammate = Gamma(Z1te,Z2te) ;
gammatm = Gamma(Z1tm,Z2tm) ;


end
