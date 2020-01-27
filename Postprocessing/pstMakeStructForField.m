function stRes = pstMakeStructForField(stAbcdnm, nNmax, lambda, epsilon2, epsilon1, stIncPar, a, c)
  %% pstMakeStructForField
% Creates the structure with necessary parameters to calculate surface fields
% 
% pstMakeStructForField(stAbcdnm, nNmax, lambda, epsilon2, epsilon1, stIncPar)
% Create a structure that may be used to evaluate the fields.
% Can also be called as pstMakeStructForField(stAbcdnm,stParams) where
% stParams is a struct with all the relevant fields.
%
% Input:
%       stAbcdnm: [struct] The field expansion coefficients, as from
%                 rvhGetFieldCoefficients. This may be for one or several
%                 wavelengths.
%       nNmax:    [1 x 1] The maximum value of N in the calculations.
%       lambda:   [L x 1] The value(s) of the wavelength used.
%       epsilon2:[L x 1] The value(s) of the dielectric function in the
%                 particle
%       epsilon1: [1 x 1] The dielectric function of the surrounding medium
%       stIncPar: [struct] A structure defining the incident plane wave, as
%                 from vshMakeIncidentParameters
%       a (optional): [1 x 1] Semi-axis of spheroid along x,y
%       c (optional): [1 x 1] Semi-axis of spheroid along z
%
% Output:
%       stRes: A struct with all the same fields as stAbcdnm, as well as
%       all of the other input parameters as fields of that struct
%
% Dependency: 
% vshMakeIncidentParameters

if ~isfield(stAbcdnm,'cnm')
    disp ('Pb in pstMakeStructForField: the structure stAbcdnm should contain cnm and dnm for internal fields');
end

stRes=stAbcdnm;

if isstruct(nNmax)
    stParams = nNmax;
    stRes.nNmax=stParams.N;
    stRes.lambda=stParams.lambda;
    stRes.epsilon2=stParams.epsilon2;
    stRes.epsilon1=stParams.epsilon1;
    stRes.a=stParams.a;
    stRes.c=stParams.c;
    % Incident field is defined by sIncType or stIncPar
    if isfield(stParams,'stIncPar')
        stRes.stIncPar = stParams.stIncPar;
    else
        sIncType=stParams.sIncType;
        stRes.stIncPar = vshMakeIncidentParameters(sIncType, stParams.N);
    end
else
    stRes.nNmax=nNmax;
    stRes.lambda=lambda;
    stRes.epsilon2=epsilon2;
    stRes.epsilon1=epsilon1;
    stRes.stIncPar=stIncPar;
    if nargin>6
        stRes.a=a;
        stRes.c=c;
    end
end
