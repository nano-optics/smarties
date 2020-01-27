function [stC, stAbcdnm, stEsurf] = slvForSurfaceField(stParams, stOptions)
  %% slvForSurfaceField
% Calculates the cross-sections, expansion coefficients, and surface field for fixed orientation
% Input:
%       - stParams: struct
%              The following parameters should be defined in stParams:
%              - a: semi-axis along x,y
%              - c: semi-axis along z
%              - k1: wavevector in embedding medium (of refractive index nM) (k1=2*pi*nM/lambda)
%              - s: relative refractive index (s=n_Particle / nM)
%              - N: number of multipoles for T-matrix
%              - nNbTheta: number of thetas for quadratures
%              - sIncType: string such a 'KxEz' defining the incident field)
%              - or stIncPar: struct defining the incident field (alternative
%              to sIncType)
%              - nNbThetaPst: number of thetas for surface field calculations
%              - lambda: wavelength (in same unit as a and c)
%              - epsilon1: dielectric constant of embedding medium
%              - epsilon2: dielectric constant of spheroid particle
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
% Output:
%       - stC: struct with cross-sections Cext, Csca, Cabs. Orientation-averaged
%               cross-section Cextoa, Cscaoa, Cabsoa are also included.
%       - stAbcdnm: structure containing the field expansion coefficients
%                anm, bnm, pnm, qnm, (and cnm, dnm if R was calculated),
%                each being a [1 x P] row.
%       - stEsurf: structure containing the surface field and averages
%                  as explained in pstSurfaceField.
%
% Dependency: 
% pstMakeStructForField, pstSurfaceField, slvForFixed, vshMakeIncidentParameters

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
% First get all expansion coefficients
[stC, stAbcdnm] = slvForFixed(stParams,stOptions);

% This now gets a result structure for postprocessing
stResE=pstMakeStructForField(stAbcdnm, stParams);

nNbThetaPst = stParams.nNbThetaPst; % number of theta for evaluating fields

stEsurf=pstSurfaceField(stResE,nNbThetaPst);

end
