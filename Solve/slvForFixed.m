function [stC, stAbcdnm] = slvForFixed(stParams, stOptions, stGeometry)
  %% slvForFixed
% Calculates the expansion coefficients and cross-sections for fixed orientation
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
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
%       - stGeometry (optional): struct with geometry
%
% Output:
%       - stC: struct with cross-sections Cext, Csca, Cabs. Orientation-averaged
%               cross-sections Cextoa, Cscaoa, Cabsoa are also included.
%       - stAbcdnm: structure containing the field expansion coefficients
%                anm, bnm, pnm, qnm, (and cnm, dnm if R was calculated),
%                each being a [1 x P] row.
%
% Dependency: 
% pstGetCrossSections, rvhGetFieldCoefficients, slvForT, vshMakeIncidentParameters

k1 = stParams.k1;
N = stParams.N;

% Incident field is defined by sIncType or stIncPar
if isfield(stParams,'stIncPar')
    stIncPar = stParams.stIncPar;
else
    sIncType=stParams.sIncType;
    stIncPar = vshMakeIncidentParameters(sIncType, N);
end

% First we need to solve for T
if nargin<3
    [stCoa, CstTRa] = slvForT(stParams,stOptions);
else
    [stCoa, CstTRa] = slvForT(stParams,stOptions,stGeometry);
end

% Get the field expansion coefficients from T and R for a given incident
% excitation
stAbcdnm = rvhGetFieldCoefficients(N, CstTRa, stIncPar);

% Calculate the (Ext, Abs, Sca) cross-sections for this excitation
stC = pstGetCrossSections(k1, stAbcdnm);

% Combine stCoa with stC:
stC.Cextoa = stCoa.Cext;
stC.Cscaoa = stCoa.Csca;
stC.Cabsoa = stCoa.Cabs;

end
