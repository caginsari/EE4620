function [krho] = findpropc(freqn,er,h,krhog,TE_TM,Strat_type) 
% [krho] = findprop(freqn,er,h,krhog,TE_TM,Strat_type) 
%Strat_type: Stratification type

lambda = 3e8./freqn ;
k0= 2*pi./lambda ;
deltk = k0/500 ;
Dkrhop = krhog.*k0 + deltk./2 ;
Dkrhom = krhog.*k0 - deltk./2 ;

T = strcmpi('GroundSlab',Strat_type) ;
if T == 1
    Dp = Den_GroundSlab(k0,er,h,Dkrhop,freqn,TE_TM) ;
    Dm = Den_GroundSlab(k0,er,h,Dkrhom,freqn,TE_TM) ;
end

deltk = k0./500 ;

Dprime = (Dp-Dm) ./ deltk ;

krho = krhog - Dkrhog./Dprime ;

end