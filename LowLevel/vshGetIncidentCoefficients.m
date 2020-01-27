function stIncEabnm = vshGetIncidentCoefficients(nMax, stIncPar)
%% vshGetIncidentCoefficients
% Calculates expansion coefficients of an incident plane wave
% 
% vshGetIncidentCoefficients(nMax,stIncPar) calculates the expansion
% coefficients for a general linearly polarized incident plane wave
% (assuming E0=1) for all m up to multipole nNmax.
% The expressions are given in Eqs. C.56-C.59) of [Mishchenko 2002].
%
% Input:
%       - nNmax:    Maximum multipole order
%       - stIncPar: A structure with three angles defining the incident
%                   wave. thetap and phip define the wavevector direction
%                   while alphap defines the field polarization
%
% Output: a structure with fields anm and bnm, each [1 x P] with the
%         expansion coefficients for all 1<=n<nNmax , -n<=m<=n
%         The indices (n,m) are condensed in a linear index (the p-index)
%         with the convention p = n(n+1) + m.
%         Note that P = nNmax (nNmax+2)
%
% Dependency: 
% vshPinmTaunm

alphap=stIncPar.alphap;
thetap=stIncPar.thetap;
phip=stIncPar.phip;

nPmax=nMax*(nMax+2); % for p-index

% First get dbarnm (see user guide)
nvec=1:nMax; % N
factn = 1i.^nvec .* realsqrt(4*pi*(2*nvec+1)./(nvec.*(nvec+1))); % [1 x N] row
mvec=-nMax:nMax;
factm= - (-1).^mvec .*exp(-1i*mvec*phip); % [1 x M] row
dbarnm = zeros(1,nPmax); % [1 x P] row
for n=1:nMax % Loop on n
    m=-n:n; % all m at a time [1 x M]
    ind=n*(n+1)+m; % [1 x M]
    dbarnm(1,ind)=factn(n) * factm(m+nMax+1); % [1 x M]
end

% now get E.B and E.C part
stPTp = vshPinmTaunm(nMax,thetap);

minusECnmstar = cos(alphap)*1i*stPTp.pinm + sin(alphap).*stPTp.taunm;
iEBnmstar = 1i*cos(alphap)*stPTp.taunm + sin(alphap)*stPTp.pinm; % [1 x P] matrix

% multiply to get final result (Eqs. C.57 and C.58 of [Mishchenko 2002])
stIncEabnm.anm = dbarnm .* minusECnmstar;
stIncEabnm.bnm = dbarnm .* iEBnmstar;
