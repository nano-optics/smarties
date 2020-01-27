%% ScriptSolveForTSpectrum
% Example script showing how to obtain the orientation-averaged optical
% properties for a spheroid as a function of wavelength.
% Plots the wavelength-dependent spectra for orientation-averaged
% extinction, scattering, and absorption cross-sections.
%%

%% Initialization
%
% Note that you need to run InitPath in the root folder first to add
% required folders to the Matlab path so that functions can be called
% Alternatively, uncomment the following line
%
%   run('..\InitPath');
%
% The following parameters should be defined:
%%
% * a: semi-axis along x,y
% * c: semi-axis along z
% * N: number of multipoles for T-matrix
% * nNbTheta: number of thetas for quadratures
% * lambda: wavelength (in same unit as a and c)
% * k1: wavevector in embedding medium (of refractive index nM) (k1=2*pi*nM/lambda)
% * s: relative refractive index (s=n_Particle / nM)
%
% lambda, k1, and s can here be wavelength-dependent vectors [L x 1]

clear
close all

%% Parameters of the scattering problem
% We define parameters for a gold nanorod in water, modeled as a prolate
% spheroid
%
% <<../fig/schematicp.png>>
%
a = 15; % in nm
c = 45; % in nm, i.e. 30 x 90nm full-axes
lambda = (400:5:900).'; % in nm
epsilon2 = epsAu(lambda);
epsilon1 = 1.33^2; % for water

%% Collect simulation parameters in a structure
k1 = 2*pi./lambda * sqrt(epsilon1);
s = sqrt(epsilon2)/sqrt(epsilon1);

stParams.a=a; stParams.c=c;
stParams.k1=k1; stParams.s=s;
stParams.lambda=lambda;
stParams.epsilon2=epsilon2;
stParams.epsilon1=epsilon1;

% Optional parameters may also be defined as follows:
stOptions.bGetR = false;
stOptions.Delta = 0;
stOptions.NB = 0; % NB will be estimated automatically
stOptions.bGetSymmetricT = false;
stOptions.bOutput = false; % false to suppress messages in lambda-loop

%% Convergence parameters

% Maximum multipole order for T-matrix and series expansions of fields
N = 15;
% Number of points for Gaussian quadratures to compute integrals in P and Q matrices
nNbTheta = 40;

% Those could also be estimated automatically for some desired accuracy as follows
% [N, nNbTheta] = sphEstimateNandNT(stParams, stOptions, 1e-15);

% Add those to the parameters structure
stParams.N=N; stParams.nNbTheta=nNbTheta;


%% T-matrix calculation
% Solve for T (all wavelengths)
tic;
stCoa = slvForTSpectrum(stParams,stOptions);

fprintf('\nT-matrices (N = %d, nNbTheta= %d) ... done in %.2f seconds.\n', N, nNbTheta, toc);

% To test for convergence and accuracy, we choose the wavelength with the largest
% k1|s| and repeat the calculation with N=N+5 and nNbTheta=nNbTheta+5
[~,indWorst]=max(abs(stParams.k1 .* stParams.s));
stParams2 = pstGetParamsStructOneLambda(stParams,lambda(indWorst));
stParams2.N=stParams2.N+5;
stParams2.nNbTheta=stParams2.nNbTheta+5;
fprintf('\nConvergence testing for lambda = %.2f.\n', lambda(indWorst));
tic;
[stCoa2, ~] = slvForT(stParams2,stOptions);
relerrExt = (abs(stCoa.Cext(indWorst)./stCoa2.Cext-1));
relerrSca = (abs(stCoa.Csca(indWorst)./stCoa2.Csca-1));
relerrAbs = (abs(stCoa.Cabs(indWorst)./stCoa2.Cabs-1));
fprintf('\nT-matrix (N = %d, nNbTheta= %d) ... done in %.2f seconds.\n', stParams2.N, stParams2.nNbTheta, toc);

%% Plotting orientation-averaged cross-sections
figure('Name','ScriptSolveForTSpectrum');
plot(lambda,[stCoa.Cext,stCoa.Csca,stCoa.Cabs]);
legend({['<Cext> (err. ', num2str(relerrExt,3),')'], ...
    ['<Csca> (err. ', num2str(relerrSca,3),')'], ...
    ['<Cabs> (err. ', num2str(relerrAbs,3),')']}, ...
    'Location','Best');
title(['a=', num2str(a), ', c=',num2str(c),', N=', int2str(N), ', Nt=', int2str(nNbTheta)]);
xlabel('Wavelength [nm]')
ylabel('Cross-section [nm^2]')
