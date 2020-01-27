function [nNbTheta, err] = sphEstimateNbTheta(stParams, stOptions,maxAcc)
%% sphEstimateNbTheta
% Estimates nNbTheta needed for accurate calculations
%
% sphEstimateNbTheta(stParams, stOptions, maxAcc)
% finds an estimate for the number of quadrature points, nNbTheta, to
% compute the integrals of the T-matrix with relative accuracy maxAcc
% If maxAcc is not specified, then the best achievable accuracy is found.
% The test is based on convergence of the orientation-averaged extinction
% cross-section taking into account only the m=0 and m=1 matrices.
% It returns the estimated converged relative precision.
%
% Input:
%       - stParams:   Structure containing simulation parameters. The
%                   maximum number of multipoles can be specified in
%                   stParams.N. If stParams.nNbTheta is defined, then the
%                   convergence test only checks for values larger than
%                   this.
%       - stOptions: struct with optional parameters, see
%              slvGetOptionsFromStruct for details.
%       - maxAcc: desired relative accuracy (default 1e-20). If this cannot
%       be reached, the best possible accuracy is returned.
%
% Dependency:
% slvForT

%set default parameters
if nargin<3
    maxAcc=1e-20; % by default, the function will look for a plateau
end
minAcc = 1e-3;
maxAcc = min(maxAcc,minAcc);


absmvec = [0,1]; % only m=0 and 1 to be faster

% This works on only one wavelength, so we choose the largest k1 * s
% as representative of the worst case
% Find max and min relative refractive index
[~,ind] = max(abs(stParams.k1 .* stParams.s));

stParam1.a=stParams.a;
stParam1.c=stParams.c;
stParam1.s =stParams.s(ind);
stParam1.k1 =stParams.k1(ind);
stOptions.absmvec=absmvec;
stParam1.N=stParams.N;

thetatable=[5:9, 10:5:50, 60:10:140, 160:20:300,350:50:700,800:100:2000];
if isfield(stParams,'nNbTheta') % if defined, only look for larger values
    thetatable=thetatable(thetatable>=stParams.nNbTheta);
end
h=max(stParams.a,stParams.c)/min(stParams.a,stParams.c);
if h>=10 % for large aspect ratio, set a minimum nNBTheta of 7*h (note this estimate is empirical)
    thetatable=thetatable(thetatable>=7*h);
end

ntmax=length(thetatable);

% Number of points used for linear fit to check if error has reached a plateau
NforConv = 4;
Afit = [ones(NforConv,1), (1:NforConv).'];
% Store errors
Qerr= zeros(ntmax,1);
Q=1e100;

warning('off', 'SMARTIES:missingm'); % suppress warnings in rvhGetAverageCrossSections

for nt = 1:ntmax % Loop over nNbTheta's
    stParam1.nNbTheta=thetatable(nt);
    [stC, ~] = slvForT(stParam1, stOptions);
    Qnew=stC.Cext;

    % Relative error
    Qerr(nt) = abs (Q./Qnew-1);
    if nt>1
%        fprintf('nt=%d, nNbTheta= %d, err= %e\n',nt,thetatable(nt-1),Qerr(nt));
        if Qerr(nt)<maxAcc % then desired convergence was achieved
            nNbTheta = thetatable(nt-1);
            err = Qerr(nt);
            warning('on', 'SMARTIES:missingm'); % reactivate warnings in rvhGetAverageCrossSections
            return;
        end
    end

    % Now check if plateau has been reached
    % Minimum requirement for convergence testing
    if nt > NforConv && max(Qerr((nt-NforConv+1):nt))<minAcc
        % Then test if the last NforConv errors have flattened out by a
        % linear regression to the last NforConv points
        coef=Afit \ log10(abs(Qerr((nt-NforConv+1):nt)));
        slope = coef(2);
        % if slope<0, then still converging
        if slope>=0 % this means less than no longer better over NforConv steps
            nNbTheta = thetatable(nt-NforConv);
            err = mean(Qerr((nt-NforConv+1):nt));
            warning('on', 'SMARTIES:missingm'); % reactivate warnings in rvhGetAverageCrossSections
            return;
        end
    end
    % Prepapre for next step
    Q=Qnew;

end
% disp('Convergence for nNbTheta was not found..')
nNbTheta = NaN;
err=NaN;
warning('on', 'SMARTIES:missingm'); % reactivate warnings in rvhGetAverageCrossSections
