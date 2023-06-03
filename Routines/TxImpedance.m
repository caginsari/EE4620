function [ZiTM,ZiTE] = TxImpedance(zetai,kzi,ki) 
% Calculates the transmission line impedance of the equivalent transmission
% line
% [ZiTM,ZiTE] = TxImpedance(zetai,kzi,ki) 

ZiTM = zetai .* kzi ./ki ;
ZiTE = zetai .* ki ./kzi ;

end