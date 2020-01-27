function stCoa = slvForTSpectrum(stParams, stOptions)
  %% slvForTSpectrum
% Calculates the T-matrix and orientation-averaged properties for multiple wavelengths
%
% Note that this function does not return the T-matrix to avoid having to
% store all matrices for all wavelengths.
% Input:
%       - stParams: struct
%              The following parameters should be defined in stParams:
%              - a: semi-axis along x,y
%              - c: semi-axis along z
%              - lambda: [L x 1] wavelength (in same unit as a and c)
%              - k1: [L x 1] or [1 x 1]
%                    wavevector in embedding medium (of refractive index nM) (k1=2*pi*nM/lambda)
%              - s: [L x 1] or [1 x 1]
%                    relative refractive index (s=n_Particle / nM)
%              - N: number of multipoles for T-matrix
%              - nNbTheta: number of thetas for quadratures
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
% Output:
%       - stCoa: struct with wavelength-dependent orientation-averaged cross-sections
%                Cext, Csca, Cabs, and lambda; each [L x 1]
%
% Dependency: 
% rvhGetAverageCrossSections, rvhGetSymmetricMat, rvhGetTRfromPQ,
% rvhTruncateMatrices, slvGetOptionsFromStruct, sphCalculatePQ,
% sphEstimateDelta, sphEstimateNB, sphMakeGeometry

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
        disp ('ERROR: Delta could not be found. Results are likely to be non-converged. Try choosing Delta manually instead.');
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

% Now main loop over lambda
tic
disp ' ';
disp (['Loop over ', int2str(L), ' lambdas...']);

stCoa = struct();
stCoa.Csca = zeros(L,1);
stCoa.Cext = zeros(L,1);
stCoa.Cabs = zeros(L,1);

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
    stCoa.Csca(lInd,1) = stCoaOneLambda.Csca;
    stCoa.Cext(lInd,1) = stCoaOneLambda.Cext;
    stCoa.Cabs(lInd,1) = stCoaOneLambda.Cabs;
end
stCoa.lambda=lambda;
end
