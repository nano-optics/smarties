function [stC, stAbcdnm, stEsurf] = slvForSurfaceFieldSpectrum(stParams, stOptions)
  %% slvForSurfaceFieldSpectrum
% Calculates the X-sections, expansion coeffs, and surface fields for many wavelengths
% Input:
%       - stParams: struct
%              The following parameters should be defined in stParams:
%              - a: semi-axis along x,y
%              - c: semi-axis along z
%              - lambda: [L x 1] wavelength (in same unit as a and c)
%              - k1: [L x 1] or [1 x 1]
%                    wavevector in embedding medium (of refractive index nM) (k1=2*pi*nM/lambda)
%              - s: [L x 1] or [1 x 1]
%              - N: number of multipoles for T-matrix
%              - nNbTheta: number of thetas for quadratures
%              - sIncType: string such a 'KxEz' defining the incident field)
%              - or stIncPar: struct defining the incident field (alternative
%              to sIncType)
%              - nNbThetaPst: number of thetas for surface field calculations
%              - epsilon1: [L x 1] or [1 x 1] dielectric constant of embedding medium
%              - epsilon2: [L x 1] or [1 x 1] dielectric constant of spheroid particle
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
% Output:
%       - stC: struct with cross-sections Cext, Csca, Cabs. Orientation-averaged
%               cross-section Cextoa, Cscaoa, Cabsoa are also included. All [L x 1].
%       - stAbcdnm: structure containing the field expansion coefficients
%                anm, bnm, pnm, qnm, (and cnm, dnm if R was calculated),
%                anm, bnm are [1 x P], others are [L x P]
%       - stEsurf: structure containing the wavelength-dependent surface field
%                  and averages as explained in pstSurfaceField.
%
% Dependency: 
% pstMakeStructForField, pstSurfaceField, slvForFixedSpectrum, vshMakeIncidentParameters

N = stParams.N;

% Incident field is defined by sIncType or stIncPar
if isfield(stParams,'stIncPar')
    stIncPar = stParams.stIncPar;
else
    sIncType=stParams.sIncType;
    stIncPar = vshMakeIncidentParameters(sIncType, N);
end
stParams.stIncPar=stIncPar;

stOptions.bGetR = true; % Will need R for surface-fields

% First get all expansion coefficients for all wavelength
[stC, stAbcdnm] = slvForFixedSpectrum(stParams,stOptions);

% This now gets a result structure for postprocessing
stResE=pstMakeStructForField(stAbcdnm, stParams);

nNbThetaPst = stParams.nNbThetaPst; % number of theta for evaluating fields

disp 'Calculating surface-fields for all wavelengths...'
% This gets the surface field for all wavelength
stEsurf=pstSurfaceField(stResE,nNbThetaPst);

end
