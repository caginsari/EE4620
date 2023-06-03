function [Vr,Ir] = Residue_GroundSlab(k0,er,h,ksw,z,f,TE_TM_flag)
%% EE4620 Assignment 3:Residue function
% Residue_GroundSlab(k0,er,h,ksw,z,air_slab_flag,TE_TM_flag)
% k0: free space propagation constant
% ksw: is the propagation constant found from the iteration method(Newton-Raphson Method)
% er: is the relative permittivity of the material Air/Slab
% air_slab_flag : the input must be Air, or Slab
% TE_TM_flag : Either TE or TM
ks = k0 *sqrt(er) ;
deltk = k0 ./500 ;
Dkrhop = ksw + deltk./2 ;
Dkrhom = ksw - deltk./2 ;

Dp = Den_GroundSlab(ks,er,h,Dkrhop,f,TE_TM_flag) ;
Dm = Den_GroundSlab(ks,er,h,Dkrhom,f,TE_TM_flag) ;
Dprime = (Dp-Dm) ./ deltk ;

[Vr,Ir] = Residue_GroundSlab_Slab(k0,ksw,er,h,z,TE_TM_flag,Dprime);
[Vr(z > h),Ir(z > h)] = Residue_GroundSlab_Air(k0,er,h,ksw,z(z > h),TE_TM_flag,Dprime);
 
end

% if 1==strcmpi(air_slab_flag,'Air')
%     [Vtmr,Itmr] = Residue_GroundSlab_Air(k0,er,h,ksw,z,TE_TM_flag,Dprime);
% elseif 1==strcmpi(air_slab_flag,'Slab')
%     [Vtmr,Itmr] = Residue_GroundSlab_Slab(k0,ksw,er,h,z,TE_TM_flag,Dprime);
% else
%     error('Must be one of the following inputs \n 1.Air \n 2.Slab \n Your input is %s',air_slab_flag)