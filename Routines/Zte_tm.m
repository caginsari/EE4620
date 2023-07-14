function [Zte,Ztm] = Zte_tm(B,theta)

%% Zte part
Zte = -1i./B .* 1./(1-sin(theta).^2./2);

Ztm = -1i./B;

end