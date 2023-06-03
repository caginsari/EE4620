function [S11] = SpectralDomains_Refl_coeff(Z, Z0)

Y = 1./Z ; 

S11= -Y.*Z0 ./(2+Y.*Z0) ;

end
