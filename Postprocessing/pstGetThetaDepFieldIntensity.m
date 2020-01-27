function M = pstGetThetaDepFieldIntensity(stEsurf, phi0, lambda0)
  %% pstGetThetaDepFieldIntensity
% Calculates surface field intensity M=|E|^2(theta) for some phi
% 
% pstPlotThetaDepFieldIntensity(stEsurf,phi0) plots surface field intensity
% M(theta) for the phi-values in phi0
%
% Input:    - stEsurf: structure with surface field
%           - phi0: vector of F phi-values
%           - lambda0 (optional): lambda at which plot is done, by default
%                                 it uses lambda(1)
% Output:   - M: [F x T] is |E(theta)|^2 for each phi-value in phi0
%
% Dependency: 
% vshEthetaForPhi

if nargin<3
    indLambda = 1;
else
    indLambda = find(stEsurf.lambda == lambda0, 1);
    if isempty(indLambda)
        indLambda = 1;
    end
end

% Calculate M=|E|^2 on the surface for all theta at phi=phi0
M=zeros(length(phi0),length(stEsurf.theta)); % [F x T]
CsLeg=cell(length(phi0),1); % For legend

if length(stEsurf.lambda)==1 % only one lambda, easiest case
    for nf=1:length(phi0)
        stEphi=vshEthetaForPhi(stEsurf,phi0(nf));
        M(nf,:)=abs(stEphi.Er).^2+abs(stEphi.Et).^2+abs(stEphi.Ef).^2; % calculate |E|^2 on the surface
        CsLeg{nf}=['phi = ', num2str(phi0(nf))];
    end
else
    % Need to extract results for one lambda
    stEsurf1.theta = stEsurf.theta;
    N=stEsurf.nNmax;
    stEsurf1.CErm=cell(1,2*N+1);
    stEsurf1.CEtm=cell(1,2*N+1);
    stEsurf1.CEfm=cell(1,2*N+1);
    for m=-N:N
        stEsurf1.CErm{1,m+N+1} = stEsurf.CErm{1,m+N+1}(indLambda,:);
        stEsurf1.CEtm{1,m+N+1} = stEsurf.CEtm{1,m+N+1}(indLambda,:);
        stEsurf1.CEfm{1,m+N+1} = stEsurf.CEfm{1,m+N+1}(indLambda,:);
    end
    for nf=1:length(phi0)
        stEphi=vshEthetaForPhi(stEsurf1,phi0(nf));
        M(nf,:)=abs(stEphi.Er).^2+abs(stEphi.Et).^2+abs(stEphi.Ef).^2; % calculate |E|^2 on the surface
    end
end

end
