function stEforPhi = vshEthetaForPhi(stEsurf, phi0)
%% vshEthetaForPhi
% Calculate E at a given phi for all theta from stEsurf struct
% 
% vshEthetaForPhi(stEsurf,phi0) calculates the E field for all lambda's and
% theta's in a plane defined by phi = phi0
%
% Input:
%       stEsurf: structure containing field
%       phi0:    [1 x 1] The value of phi at which the field is to be
%                calculated
%
% Output:
%        stEforPhi: structure with 3 fields, each [L x T]
%           - Er: The radial field component
%           - Et: The theta field component
%           - Ef: The phi field component
%           and with additional fields:
%           - theta [1 x T]
%           - phi0 [1 x 1]
%
% Dependency: 
% none

nNmax = (length(stEsurf.CErm)-1)/2;
theta = stEsurf.theta;

nNbLambda=size(stEsurf.CErm{1},1);
nNbTheta = length(theta);

stEforPhi.Er = zeros(nNbLambda,nNbTheta);
stEforPhi.Et = zeros(nNbLambda,nNbTheta);
stEforPhi.Ef = zeros(nNbLambda,nNbTheta);
stEforPhi.theta = theta;
stEforPhi.phi0 = phi0;

for m =-nNmax:nNmax

    expphase = exp(1i*m*phi0);

    stEforPhi.Er = stEsurf.CErm{1,m+nNmax+1}*expphase + stEforPhi.Er;
    stEforPhi.Et = stEsurf.CEtm{1,m+nNmax+1}*expphase + stEforPhi.Et;
    stEforPhi.Ef = stEsurf.CEfm{1,m+nNmax+1}*expphase + stEforPhi.Ef;

end;
end
