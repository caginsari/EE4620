function [krhogn,krho] = IterativeMethod(h,hs,er,f,stratification_type,TE_TM_flag,krho)
%% EE4620 Assignment 3: [krhogn,krho] = IterativeMethod(h,er,f,type,TE_TM_flag,krho)
% f is frequency in Hz a single value not an array.
% Solution: solution could be Numerical or Analytical.

k0 = @(freq) 2.*pi.* freq./3e8 ;
switch stratification_type
    case 'GroundSlab'
        [D] = Den_GroundSlab(k0(f),er,h,krho,f,TE_TM_flag) ;
        [y, peak] = findpeaks( abs(1 ./ D) ) ;

        if length(peak)>1
            smallPeakIndexes = y < 0.8*max(y);
            y(smallPeakIndexes) = [] ; %Reject Y value of peaks below this threshold
            peak(smallPeakIndexes) = [] ; %Reject X value of peaks below this threshold

        end
        krho = krho(peak);

        %% No Peaks
        if isempty(krho)
            krho_old = k0(f) ;
            krho = k0(f) ;
            krhogn = k0(f) ./ k0(f) ;
        else
            krho_old = 0 ;
        end
    case 'SuperStrate'
        krho = krhogSuperStrate(k0(f), er, h, TE_TM_flag ) ;
        krho_old = 0 ;
    case 'SemiInf'
        krho = krhogSemiInf(k0(f), er, h, TE_TM_flag ) ;
        krho_old = 0;
end



%% Propagation Constant at highest frequency
% krho = krhogn_0 .* k0(f) ;

%% Iteration part
if 1 == strcmpi('GroundSlab',stratification_type)
    while abs(krho-krho_old) > 0.001
        k0_it = k0(f) ;

        deltk = k0_it ./500 ;
        Dkrhop = krho + deltk./2 ;
        Dkrhom = krho - deltk./2 ;

        Dp = Den_GroundSlab(k0_it,er,h,Dkrhop,f,TE_TM_flag) ;
        Dm = Den_GroundSlab(k0_it,er,h,Dkrhom,f,TE_TM_flag) ;
        Dkrhog = Den_GroundSlab(k0_it,er,h,krho,f,TE_TM_flag) ;
        Dprime = (Dp-Dm) ./ deltk ;

        krho_old = krho ;
        krho = krho - Dkrhog./Dprime ;
        krhogn = krho ./ k0(f) ;
        if imag(krho) ~= 0
            krho = k0_it;
            break;
        end
    end
elseif 1 == strcmpi('SuperStrate',stratification_type)
     while abs(krho-krho_old) > 0.0001
        k0_it = k0(f) ;

        deltk = k0_it ./500 ;
        Dkrhop = krho + deltk./2 ;
        Dkrhom = krho - deltk./2 ;

        Dp = Den_SuperStrateNumeric(k0_it,er,h,hs,Dkrhop,TE_TM_flag) ;
        Dm = Den_SuperStrateNumeric(k0_it,er,h,hs,Dkrhom,TE_TM_flag) ;
        Dkrhog = Den_SuperStrateNumeric(k0_it,er,h,hs,krho,TE_TM_flag) ;
        Dprime = (Dp-Dm) ./ deltk ;

        krho_old = krho ;
        krho = krho - Dkrhog./Dprime ;
        krhogn = krho ./ k0(f) ;
    end
elseif 1 == strcmpi('SemiInf',stratification_type)
     while abs(krho-krho_old) > 0.0001
        k0_it = k0(f) ;

        deltk = k0_it ./500 ;
        Dkrhop = krho + deltk./2 ;
        Dkrhom = krho - deltk./2 ;

        Dp = Den_SemiInfNumeric(k0_it,er,h,Dkrhop,TE_TM_flag) ;
        Dm = Den_SemiInfNumeric(k0_it,er,h,Dkrhom,TE_TM_flag) ;
        Dkrhog = Den_SemiInfNumeric(k0_it,er,h,krho,TE_TM_flag) ;
        Dprime = (Dp-Dm) ./ deltk ;

        krho_old = krho ;
        krho = krho - Dkrhog./Dprime ;
        krhogn = krho ./ (k0(f).*sqrt(er) ) ;
     end
else
    error('Error. Input must be one of the following: \n 1. GroundSlab \n 2. SuperStrate\n 3. SemiInf \n Your input is %s',stratification_type) ;
end


end