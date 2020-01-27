function stCrossSection = pstGetCrossSections(k1, stAbcdnm)
%% pstGetCrossSections
% Calculates cross sections from expansion coefficients
% 
% pstGetCrossSections(k1, stAbcdnm) calculates absorption, scattering and
% extinction cross sections from the expansion coefficients (i.e. for a
% fixed orientation. This may be called for multiple wavelengths.
%
% Input:
%		 - k1			incident wavevector [L x 1]
%		 - stAbcdnm   structure with expansion coefficients including:
%                   - anm, bnm [L x P] or [1 x P] (incident field)
%                   - pnm, qnm [L x P] (scattered field)
%
% Output:
%       stCrossSection: Structure with three fields
%           - Cext: [L x 1] wavelength-dependent extinction cross-section
%           - Csca: [L x 1] wavelength-dependent scattering cross-section
%           - Cabs: [L x 1] wavelength-dependent absorption cross-section
%
% Dependency: 
% none

pnm = stAbcdnm.pnm; % [L x P]
qnm = stAbcdnm.qnm; % [L x P]
anm = stAbcdnm.anm; % [1 x P] (if lambda-independent) or [L x P]
bnm = stAbcdnm.bnm; % [1 x P] or [L x P]

% Eq. 5.18b of [Mishchenko 2002]
CscaSum=sum(abs(pnm).^2,2) + sum(abs(qnm).^2,2); % [L x 1]

% Eq. 5.18a of [Mishchenko 2002]
CextSum = sum(bsxfun(@times,conj(pnm),anm),2) + sum(bsxfun(@times,conj(qnm),bnm),2);

stCrossSection.Csca= 1./((k1).^2) .* CscaSum; % [L x 1]

stCrossSection.Cext = -1./((k1).^2) .* real(CextSum); % [L x 1]

stCrossSection.Cabs = stCrossSection.Cext - stCrossSection.Csca; % [L x 1]

end
