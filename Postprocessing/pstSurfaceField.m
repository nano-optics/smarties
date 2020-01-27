function stEsurf = pstSurfaceField(stRes, stRtfunc, stPinmTaunm)
  %% pstSurfaceField
% Calculates electric field on the surface
% 
% pstSurfaceField(stRes,stRtfunc,stPinmTaunm) calculates the electric field
% (for each m) on the surface defined by stRtfunc, which should be the
% surface of the particle (but may have more points than used for the
% integration). It also returns a collection of averaged field values.
% These field values are for just outside the particle, and include
% averages of the components in the parallel and perpendicular directions.
%
% Input:
%           stRes: A structure containing parameters of the simulation,
%                  with fields lambda [L x 1], epsilon2 [L x 1], epsilon1
%                  [1 x 1], stIncPar, cnm [L x P], dnm [L x P],
%                  nNmax [1 x 1], a [1 x 1], and c [1 x 1]
%           stRtfunc (or nNbPts): A structure containing the geometry, which is valid
%                     for theta on [0;pi], as from rvhGetThetasForAverage
%                     If stRtfunc.sInt='Pts', then no average are
%                     calculated.
%                     If an integer, it is interpreted as nNbPts, and stRtfunc
%                     is recalculated with nNbPts points
%           stPinmTaunm: (optional) The angular functions as from
%                        vshPinmTaunm, to save computing these multiple
%                        times. Else is calculated by this function
%
% Output:
%           stEsurf: A structure containing fields
%               - CErm: Cell {2N+1 x 1} each [L x T] theta-dependent radial (in
%                       spherical coordinates) electric field, with each
%                       cell entry being the contribution from a different
%                       m, from m=-N to m=N
%               - CEtm: Cell {2N+1 x 1}, each [L x T] Same for E_theta
%               - CEfm: Cell {2N+1 x 1}, each [L x T] Same for E_phi
%               - stRtfunc: structure defining the geometry
%               - theta: [1 x T] The values of theta used to calculate the fields
%               - roftheta:    [1 x T] The values of r at which the field is
%                        calculated
%               - lambda:[L x 1] The wavelengths used (from stRes)
%               - stIncPar: Structure containing incident parameters (from
%                           stRes)
%               - nNmax: [1 x 1] Maximum value of N, as from stRes
%               - epsilon: [L x 1] Dielectric function inside the scatterer
%               - S:     [1 x 1] The surface area of the scatterer
%               - MLocPerpAve: [L x 1] The average perpendicular M
%               - MLocParaAve: [L x 1] The average parallel M
%               - MLocAve: [L x 1] The average M=|E|^2
%               - EAve:   [L x 1] The average |E|
%               - F0E4Ave: [L x 1] The average F^0_E^4=|E|^4
%
% Dependency: 
% rvhGetThetasForAverage, sphMakeGeometry, vshEFaverages [private],
% vshEgenThetaAllPhi, vshPinmTaunm, vshSquareCEm [private]

% Test if stRtfunc is structure or needs to be recalculated.
if ~isstruct(stRtfunc)
    nNbPts=stRtfunc;
    stRtfunc=sphMakeGeometry(nNbPts,stRes.a,stRes.c);
    stRtfunc=rvhGetThetasForAverage(stRtfunc); % get thetas over entire range [0,pi]
end

% get theta dependence if not provided
if nargin < 3
    stPinmTaunm=vshPinmTaunm(stRes.nNmax,stRtfunc.theta); % fields are [T x nNmax]
end

rRow = transpose(stRtfunc.r); % [1 x T]
drdtRow=transpose(stRtfunc.drdt); % [1 x T]
thetaRow=transpose(stRtfunc.theta); % [1 x T]

% Get field inside on the surface using the VSH series expansion with
% coefficients cnm and dnm and regular VSH (with 'j')
stEsurf=vshEgenThetaAllPhi(stRes.lambda,stRes.epsilon2,stRes.cnm,stRes.dnm, ...
    rRow,thetaRow,'j',stPinmTaunm);

nNmax=stRes.nNmax;

% Apply boundary condition

% Note that normal is n=n_r e_r + n_t e_t with
% nr=(r./realsqrt(r.^2+drdt.^2)); % [1 x T]
% nt=(-drdt./realsqrt(r.^2+drdt.^2)); % [1 x T]
% tangential unit is h=-n_t e_r + n_r e_t
% nothing is changed along phi
s2m1=((0*stRes.lambda+stRes.epsilon2)./stRes.epsilon1) - 1; % (s^2-1) [L x 1]
nrnt = (-rRow.*drdtRow./(rRow.^2+drdtRow.^2)); % [1 x T]
nr2 = (rRow.^2./(rRow.^2+drdtRow.^2)); % [1 x T]

Arr = 1 + (s2m1 * nr2); % [L x T]
Art = s2m1 * nrnt; % [L x T]
Att = 1 + (s2m1 * (1-nr2)); % [L x T]

for mind=1:(2*nNmax+1)
    % Rewrite outside fields from inside fields using boundary conditions
    ErIn = stEsurf.CErm{1,mind};
    EtIn = stEsurf.CEtm{1,mind};
    stEsurf.CErm{1,mind}=Arr .* ErIn + Art .* EtIn;
    stEsurf.CEtm{1,mind}=Art .* ErIn + Att .* EtIn;
    % CEfm is unchanged
end


% Fill in other fields
stEsurf.stRtfunc=stRtfunc;
stEsurf.theta=stRtfunc.theta;
stEsurf.roftheta=stRtfunc.r;
stEsurf.lambda=stRes.lambda;
stEsurf.stIncPar=stRes.stIncPar;
stEsurf.nNmax = nNmax;
stEsurf.epsilon = stRes.epsilon2;

if ~strcmpi(stRtfunc.sInt,'Pts')
% Calculate surface-averages
    stEsurf=vshEFaverages(stEsurf,stRtfunc);
end

end


function stEsurf=vshEFaverages(stEsurf,stRtfunc)
% Calculate average enhancement factors on surfaces
% vshEFaverages(stEsurf,stRtfunc) calculates
% surface-averaged quantities of relevance from the field calculated
% either from vshEgenThetaAllPhi.
%
% Input:
%        stEsurf: contains surface fields and other infos
%        stRtfunc: structure defining the geometry, which must be the same
%                  as that used to calculate stEsurf
%
% Output:
%   stEsurf: same structure stEsurf with the additional fields:
% - MLocParaAve: [LR x 1] wavelength- or r-dependent average MLocPara
% - MLocPerpAve: [LR x 1] wavelength- or r-dependent average MLocPerp
% - MLocAve: [LR x 1] wavelength- or r-dependent average MLoc
% - F0E4Ave: [LR x 1] wavelength- or r-dependent average SERS EF F^0_{E4}
% - Eave: [LR x 1] wavelength- or r-dependent average |E|
%

r = stRtfunc.r; % [T x 1]
drdt = stRtfunc.drdt; % [T x 1]
dS = 2*pi*r.*realsqrt(r.^2+(drdt).^2) ...
    .*stRtfunc.wTheta; % [T x 1]
S=sum(dS,1); % [1 x 1]
stEsurf.S=S;
dSovS = dS / S; % [T x 1]

dSdivOvS = (dS./(r.^2+(drdt).^2)) / S; % [T x 1]

rRow = transpose(r); % [1 x T]
drdtRow = transpose(drdt);

nNmax=stEsurf.nNmax;

stEsurf.MLocPerpAve=0;
stEsurf.MLocAve=0;
stEsurf.EAve=0;
for m =-nNmax:nNmax

    Em2theta=abs(stEsurf.CErm{1,m+nNmax+1}).^2 + ...
    abs(stEsurf.CEtm{1,m+nNmax+1}).^2+ ...
    abs(stEsurf.CEfm{1,m+nNmax+1}).^2; % [L x T]

    % Do surface-integrals as sums using matrix-vector products
    stEsurf.MLocAve = stEsurf.MLocAve + Em2theta * dSovS; % [L x 1]
    stEsurf.EAve = stEsurf.EAve + sqrt(Em2theta) * dSovS; % [L x 1]

    stEsurf.MLocPerpAve = stEsurf.MLocPerpAve + ...
         (abs(   bsxfun(@times,stEsurf.CErm{1,m+nNmax+1},rRow) ...
        - bsxfun(@times,stEsurf.CEtm{1,m+nNmax+1},drdtRow)    ).^2) ...
        * dSdivOvS; % [L x 1]
end;

stEsurf.MLocParaAve=stEsurf.MLocAve-stEsurf.MLocPerpAve;

% average SERS EF =<|E|^4> (E4 approximation)

CErm2=vshSquareCEm(stEsurf.CErm);
CEtm2=vshSquareCEm(stEsurf.CEtm);
CEfm2=vshSquareCEm(stEsurf.CEfm);

stEsurf.F0E4Ave = (abs(CErm2{1,1}+CEtm2{1,1}+CEfm2{1,1}).^2) * dSovS;


for qq=1:(2*nNmax)
    stEsurf.F0E4Ave = stEsurf.F0E4Ave + ...
        (2*abs(CErm2{1,qq+1}+CEtm2{1,qq+1}+CEfm2{1,qq+1}).^2) * dSovS;
end

end

function CEm2=vshSquareCEm(CEm)
% Calculates the sum of the squares of the field components
% vshSquareCEm(CEm) Calculates the values of the modulus squared of
% whatever CEm represents (eg E or M) for different components in m.
% That is, if some quantity X is defined as
%       X = \sum^m=N_m=-N X_m e^{im*phi}
% then
%       |X|^2 = \sum^q=2N_q=0 Y_q e^{iq*phi}
% where
%       Y_m = X_m X*_-m
% This method returns the partical sums of Y_m from m=0 to m=2N as a cell.
% Typically this is used to calculate M=|E|^2 or F=|M|^2=|E|^4.
%
% Input:
%        CEm: [1 x M] Cell where M=2N+1, for m from -N to N. This contains
%        the m components of the quantity being considered.
%
% Output:
%        CEm2:[1 x M] Cell, where M=2N+1, for q from 0 to 2N. The partial
%        sums of CEm{m}*conj(CEm{-m})
%

% CEm is cell 1 x M for m=-N:N
% CEm{m+N+1} is for m
% CEm2 is cell 1 x M for q=0:2N
% CEm2{q+1} is for q

nNbM=length(CEm);
nNmax=(nNbM-1)/2;

CEm2=cell(1,nNbM);

for qq=0:(2*nNmax)
    Sq=0;
    for tt=(qq+1):nNbM
        Sq=Sq+CEm{1,tt}.*conj(CEm{1,tt-qq});
    end
    CEm2{1,qq+1}=Sq;
end
end
