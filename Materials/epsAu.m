function epsAu = epsAu(lambda)
%% epsAu
% Returns the wavelength-dependent relative dielectric function of gold (Au)
% 
% This function uses the analytical expression given in Eq. (E.2)
% of Principles of surface-enhanced Raman spectroscopy (Elsevier 2009).
% The exp(-i omega t) convention is assumed.
%
% Parameters:
% - lambda: scalar, vector, or matrix
%           wavelengths in NANOMETERS (nm)
% 
% Returns:
% - epsAu: same size as lambda
%          epsilon(lambda) as a complex number.
%
% from Eq. E.2 of Principles of surface-enhanced Raman spectroscopy
% (Elsevier 2009)
%
% Dependency: 
% none

eps_infty = 1.54;
lambda_p = 177.5;
mu_p = 14500.0;
A1=1.27;
lambda1=470.0;
mu_p1=1900.0;
A2=1.1;
lambda2=325.0;
mu_p2=1060.0;
phi=-pi/4;
 
epsAu = eps_infty * (1 - 1 ./ (lambda_p.^2 *( (1./lambda).^2 + 1i./(mu_p.*lambda))))...
    + A1 / lambda1 *(exp(1i*phi)./(1/lambda1-1./lambda-1i/mu_p1)+exp(-1i*phi)./(1/lambda1+1./lambda+1i/mu_p1))...
    + A2 / lambda2 *(exp(1i*phi)./(1/lambda2-1./lambda-1i/mu_p2)+exp(-1i*phi)./(1/lambda2+1./lambda+1i/mu_p2));

