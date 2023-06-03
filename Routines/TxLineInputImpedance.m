function Zin = TxLineInputImpedance(Z1,Z2,k1,h)
%EE4620 Assignment 1: Transmission line input impendance function.
% Zin = TxLineInputImpedance(Z1,Z2,k1,h)

num = Z2 + 1i.* Z1 .*tan(k1.*h) ;
denum = Z1 + 1i.* Z2 .*tan(k1.*h) ;

Zin = Z1 .* num ./denum ;

end
