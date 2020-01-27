function [stC, stAbcdnm] = slvForFixedSpectrum(stParams, stOptions)
  %% slvForFixedSpectrum
% Calculates the cross-sections and expansion coefficients
% for fixed orientation and multiple wavelengths
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
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
% Output:
%       - stC: struct with cross-sections Cext, Csca, Cabs. Orientation-averaged
%               cross-section Cextoa, Cscaoa, Cabsoa are also included. All [L x 1].
%       - stAbcdnm: structure containing the field expansion coefficients
%                anm, bnm, pnm, qnm, (and cnm, dnm if R was calculated),
%                anm, bnm are [1 x P], others are [L x P]
%
% Dependency: 
% pstGetCrossSections, rvhGetAverageCrossSections, rvhGetFieldCoefficients,
% rvhGetSymmetricMat, rvhGetTRfromPQ, rvhTruncateMatrices,
% slvGetOptionsFromStruct, sphCalculatePQ, sphEstimateDelta,
% sphEstimateNB, sphMakeGeometry, vshGetIncidentCoefficients,
% vshMakeIncidentParameters

c = stParams.c;
a = stParams.a;

lambda = stParams.lambda; % [L x 1]
s = stParams.s; % [L x 1]
k1 = stParams.k1; % [L x 1]
L=length(lambda);

% For convenience, k1 and s are stored in a struct
stk1s.k1=k1;
stk1s.s=s;

N = stParams.N;
nNbTheta = stParams.nNbTheta;

[bGetR,Delta,NB,absmvec,bGetSymmetricT,bOutput] = slvGetOptionsFromStruct(stParams,stOptions);

% Make structure describing spheroidal geometry and quadrature points for
% numerical integrations
stGeometry = sphMakeGeometry(nNbTheta, a, c);

if Delta<0 % then need to estimate Delta
    % We here estimate delta using only one lambda, for which the product
    % k1|s| is largest
    [Delta, T2211err]= sphEstimateDelta(stGeometry, stk1s);
    if isnan(Delta)
        warning('ERROR: Delta could not be found. Results are likely to be non-converged. Try choosing Delta manually instead.');
        return;
    end
    disp (['Delta estimated to \Delta=', int2str(Delta),' with relative error in T_{11}^{22,m=1} of ',  num2str(T2211err)]);
end

NQ = N+Delta;% NQ>=N: Maximum multipole order for computing P and Q matrices

% Estimating NB
if NB<=0
    NB=sphEstimateNB(NQ, stGeometry, stk1s);
end
if NB<NQ
    NB=NQ; % NB must be at least NQ
end

% Incident field is defined by sIncType or stIncPar
if isfield(stParams,'stIncPar')
    stIncPar = stParams.stIncPar;
else
    sIncType=stParams.sIncType;
    stIncPar = vshMakeIncidentParameters(sIncType, N);
end


% Now main loop over lambda
fprintf('\nLoop over %d lambda values...\n', L);

% Coefficients of incident wave do not change with wavelength
stIncEabnm=vshGetIncidentCoefficients(N,stIncPar);

% Initialize variables to store results for each wavelength
stAbcdnm.anm = stIncEabnm.anm; % [1 x P] where P=N(N+2)
stAbcdnm.bnm = stIncEabnm.bnm; % [1 x P] where P=N(N+2)
P=N*(N+2); % number of elements in p-index
L=length(lambda);
stAbcdnm.pnm = zeros(L,P);
stAbcdnm.qnm = zeros(L,P);
if bGetR
    stAbcdnm.cnm = zeros(L,P);
    stAbcdnm.dnm = zeros(L,P);
end

Cscaoa = zeros(L,1);
Cextoa = zeros(L,1);
Cabsoa = zeros(L,1);

stParams1.bOutput = bOutput;
for lInd=1:L
    if bOutput
        disp(['lambda = ', num2str(lambda(lInd))]);
    end

    stParams1.k1 = stParams.k1(lInd);
    stParams1.s = stParams.s(lInd);
    % This calculates P and Q
    CstPQa = sphCalculatePQ(NQ, absmvec, stGeometry, stParams1, NB);

    % Get T (and possibly R)
    CstTRa = rvhGetTRfromPQ(CstPQa,bGetR);

    % If needed, discard higher order multipoles
    % (which are affected by the finite size of P and Q)
    if NQ>N
        CstTRa = rvhTruncateMatrices(CstTRa, N);
    end
    % T and R matrices now go up to N multipoles

    % If required, symmetrize the T-matrix
    if bGetSymmetricT
        CstTRa = rvhGetSymmetricMat(CstTRa, {'st4MT'});
    end

    % Calculate the (Ext, Abs, Sca) cross-sections for orientation-averaged excitation
    stCoaOneLambda = rvhGetAverageCrossSections(stParams1.k1, CstTRa);
    Cscaoa(lInd,1) = stCoaOneLambda.Csca;
    Cextoa(lInd,1) = stCoaOneLambda.Cext;
    Cabsoa(lInd,1) = stCoaOneLambda.Cabs;

    % Get the field expansion coefficients from T and R for a given incident
    % excitation (defined earlier in stIncPar)
    stAbcdnmOneLambda = rvhGetFieldCoefficients(N, CstTRa, stIncPar, stIncEabnm);
    stAbcdnm.pnm(lInd,:) = stAbcdnmOneLambda.pnm;
    stAbcdnm.qnm(lInd,:) = stAbcdnmOneLambda.qnm;
    if bGetR % internal fields only if R has been calculated
        stAbcdnm.cnm(lInd,:) = stAbcdnmOneLambda.cnm;
        stAbcdnm.dnm(lInd,:) = stAbcdnmOneLambda.dnm;
    end

end

% Calculate the (Ext, Abs, Sca) cross-sections for all wavelength in one go
stC = pstGetCrossSections(k1, stAbcdnm);

% Combine orientation-averaged cross-sections with stC:
stC.Cextoa = Cextoa;
stC.Cscaoa = Cscaoa;
stC.Cabsoa = Cabsoa;
stC.lambda = lambda;
end
