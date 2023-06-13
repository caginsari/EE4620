function [krhogn,krho] = IterativeMethod(h,er,f,type,TE_TM_flag,krho)
%% EE4620 Assignment 3: [krhogn,krho,Dprime] = IterativeMethod(h,er,f,type,TE_TM_flag,krho)
% f is frequency in Hz a single value not an array.

k0 = @(freq) 2.*pi.* freq./3e8 ;

if 1 == strcmpi('GroundSlab',type)
    [D] = Den_GroundSlab(k0(f),er,h,krho,f,TE_TM_flag) ;
    [y, peak] = findpeaks( abs(1 ./ D) );

    if length(peak)>1
        smallPeakIndexes = y < 0.5*max(y);
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

elseif 1 == strcmpi('SuperStrate',type)
    krho = krhogSuperStrate(k0, er, h, TE_TM_flag) ;
else
    error('Error. Input must be one of the following: \n 1. GroundSlab \n 2. SuperStrate \n Your input is %s',type) ;
end


%% Propagation Constant at highest frequency
% krho = krhogn_0 .* k0(f) ;

%% Iteration part
if 1 == strcmpi('GroundSlab',type)
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
elseif 1 == strcmpi('SuperStrate',type)
     while abs(krho-krho_old) > 0.001
        k0_it = k0(f) ;

        deltk = k0_it ./500 ;
        Dkrhop = krho + deltk./2 ;
        Dkrhom = krho - deltk./2 ;

        Dp = Den_SuperStrate(k0_it,er,h,Dkrhop,f,TE_TM_flag) ;
        Dm = Den_SuperStrate(k0_it,er,h,Dkrhom,f,TE_TM_flag) ;
        Dkrhog = Den_SuperStrate(k0_it,er,h,krho,f,TE_TM_flag) ;
        Dprime = (Dp-Dm) ./ deltk ;

        krho_old = krho ;
        krho = krho - Dkrhog./Dprime ;
        krhogn = krho ./ k0(f) ;
        if imag(krho) ~= 0
            krho = k0_it;
            break;
        end
    end
else
    error('Error. Input must be one of the following: \n 1. GroundSlab \n 2. SuperStrate \n Your input is %s',type) ;
end


end