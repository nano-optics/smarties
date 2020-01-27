function epsAg = epsAg(lambda)
%% epsAg
% Returns the wavelength-dependent relative dielectric function of silver (Ag)
% 
% This function uses the analytical expression given in Eq. (E.1).
% of Principles of surface-enhanced Raman spectroscopy (Elsevier 2009).
% The exp(-i omega t) convention is assumed.
%
% Parameters:
% - lambda: scalar, vector, or matrix
%           wavelengths in NANOMETERS (nm)
% 
% Returns:
% - epsAg: same size as lambda
%          epsilon(lambda) as a complex number.
%
% from Eq. E.1 of Principles of surface-enhanced Raman spectroscopy
% (Elsevier 2009)
%
% Dependency: 
% none

eps_infty = 4.0;
lambda_p = 282.0;
mu_p = 17000.0;
epsAg = eps_infty *(1-1./(lambda_p.^2 *( (1./lambda).^2 + 1i./(mu_p.*lambda))));

