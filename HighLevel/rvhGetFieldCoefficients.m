function stAbcdnm = rvhGetFieldCoefficients(nNmax, CstTRa, stIncPar, stIncEabnm)
%% rvhGetFieldCoefficients
% Calculates the field coefficients from the T-(R-)matrix for a fixed orientation
% 
% rvhGetFieldCoefficients(nNmax, CstTRa, stIncPar, stIncEabnm) calculates the field
% expansion coefficients from the T/R-matrices for a given incident plane wave.
% This calculate the coefficients for one wavelength only.
% If R is not defined in CstTRa, then the internal fields are not computed.
% If stIncEabnm is given, the expansion coefficients for the incident wave
% are not recalculated.
% This method is valid for scatterers with rvh symmetry.
%
% If CstTRa contains m values that are not required, then they are ignored.
%
% Input:
%	nNmax:		The maximum order to use in the expansions (should be less than or
%               equal to the T-matrix size)
%	CstTR:		A cell containing the T and R matrix for each m [1 x M]. This
%               should contain all m present in stIncPar.absmvec, and
%               should make use of the rvh symmetry
%   stIncPar:	A structure containing information about the incident
%				field, as from vshMakeIncidentParameters
%   stIncEabnm (optional):	A structure containing the expansion
%                           coefficients anm and bnm of the incident wave
%                           These are recalculated if not provided
%
% Output:
%	stAbcdnm:	A structures with fields anm, bnm, cnm, dnm, pnm,
%               qnm (all [1 x P]), where P=nNmax(nNmax+2).
%               cnm and dnm are only included if R is defined in CstTRa.
%
% Dependency: 
% vshGetIncidentCoefficients

% Coeff of incident wave
if nargin<4
    stIncEabnm=vshGetIncidentCoefficients(nNmax,stIncPar);
end

if size(CstTRa,1)~=1
    disp('CstTRa should be {1 x M} for rvhGetFieldCoefficients, i.e. one wavelength at a time');
end
bGetR = false;
if length(CstTRa{1}.CsMatList)>1
    if strcmp(CstTRa{1}.CsMatList{2},'st4MR') % R-matrix is included
        bGetR = true;
    end
end

% Define incident wave parameters
absmvec=stIncPar.absmvec(stIncPar.absmvec<=nNmax); % truncate to relevant m's

% number of m-values in T, not all of which might be needed for calculations
numM = size(CstTRa,2);

M = length(absmvec);

% dimension of (n,m) coupled indices
nPmax = nNmax*(nNmax+2);


anm = stIncEabnm.anm;   %[1 x P]
bnm = stIncEabnm.bnm;

acol = transpose(anm);   % [P x 1]
bcol = transpose(bnm);   % [P x 1]

pcol = zeros(nPmax, 1);  % [P x 1]
qcol = zeros(nPmax, 1);

if bGetR
    ccol = zeros(nPmax, 1);  % [P x 1]
    dcol = zeros(nPmax, 1);
end

mindforT = zeros(1, M);

% Find the corresponding m-indices in CstTRa (in case they are not in
% standard order)
for mind11 = 1:M
    for mind=1:numM
        if CstTRa{1,mind}.st4MTeo.m==absmvec(mind11)
            mindforT(mind11)=mind;
            break
        end
    end
end

for mind11 = 1:M
    mind = mindforT(mind11);
    m=absmvec(mind11);

    nmin=max(m,1);
    evenodd=mod(nmin,2);
    nvece = (nmin+evenodd):2:nNmax; % indices for even n [1 x Ne]
    nveco = (nmin+1-evenodd):2:nNmax; % indices for odd n [1 x No]

    pvece = nvece.*(nvece+1)+m; %[1 x Ne]
    pveco = nveco.*(nveco+1)+m; %[1 x No]
    nvecTe=1:length(nvece); %(1:Ne)
    nvecTo=1:length(nveco); %(1:No)

    pcol(pvece, 1) = CstTRa{1,mind}.st4MTeo.M11(nvecTe,nvecTe) * acol(pvece,1) ...
        + CstTRa{1,mind}.st4MTeo.M12(nvecTe,nvecTo) * bcol(pveco,1);
    qcol(pveco, 1) = CstTRa{1,mind}.st4MTeo.M21(nvecTo,nvecTe) * acol(pvece,1) ...
        + CstTRa{1,mind}.st4MTeo.M22(nvecTo,nvecTo) * bcol(pveco,1);

    pcol(pveco, 1) = CstTRa{1,mind}.st4MToe.M11(nvecTo,nvecTo) * acol(pveco,1) ...
        + CstTRa{1,mind}.st4MToe.M12(nvecTo,nvecTe) * bcol(pvece,1);
    qcol(pvece, 1) = CstTRa{1,mind}.st4MToe.M21(nvecTe,nvecTo) * acol(pveco,1) ...
        + CstTRa{1,mind}.st4MToe.M22(nvecTe,nvecTe) * bcol(pvece,1);

    if bGetR
        ccol(pvece, 1) = CstTRa{1,mind}.st4MReo.M11(nvecTe,nvecTe) * acol(pvece,1) ...
            + CstTRa{1,mind}.st4MReo.M12(nvecTe,nvecTo) * bcol(pveco,1);
        dcol(pveco, 1) = CstTRa{1,mind}.st4MReo.M21(nvecTo,nvecTe) * acol(pvece,1) ...
            + CstTRa{1,mind}.st4MReo.M22(nvecTo,nvecTo) * bcol(pveco,1);

        ccol(pveco, 1) = CstTRa{1,mind}.st4MRoe.M11(nvecTo,nvecTo) * acol(pveco,1) ...
            + CstTRa{1,mind}.st4MRoe.M12(nvecTo,nvecTe) * bcol(pvece,1);
        dcol(pvece, 1) = CstTRa{1,mind}.st4MRoe.M21(nvecTe,nvecTo) * acol(pveco,1) ...
            + CstTRa{1,mind}.st4MRoe.M22(nvecTe,nvecTe) * bcol(pvece,1);
    end

    if (m~=0) % need to do m<0 also
        pvecne = nvece.*(nvece+1)-m; %[1 x Ne]
        pvecno = nveco.*(nveco+1)-m; %[1 x No]
        % For negative m, using Eq. 5.37 of [Mishchenko 2002]
        pcol(pvecne, 1) = CstTRa{1,mind}.st4MTeo.M11(nvecTe,nvecTe) * acol(pvecne,1) ...
            - CstTRa{1,mind}.st4MTeo.M12(nvecTe,nvecTo) * bcol(pvecno,1);
        qcol(pvecno, 1) = - CstTRa{1,mind}.st4MTeo.M21(nvecTo,nvecTe) * acol(pvecne,1) ...
            + CstTRa{1,mind}.st4MTeo.M22(nvecTo,nvecTo) * bcol(pvecno,1);

        pcol(pvecno, 1) = CstTRa{1,mind}.st4MToe.M11(nvecTo,nvecTo) * acol(pvecno,1) ...
            - CstTRa{1,mind}.st4MToe.M12(nvecTo,nvecTe) * bcol(pvecne,1);
        qcol(pvecne, 1) = - CstTRa{1,mind}.st4MToe.M21(nvecTe,nvecTo) * acol(pvecno,1) ...
            + CstTRa{1,mind}.st4MToe.M22(nvecTe,nvecTe) * bcol(pvecne,1);

        if bGetR
            ccol(pvecne, 1) = CstTRa{1,mind}.st4MReo.M11(nvecTe,nvecTe) * acol(pvecne,1) ...
                - CstTRa{1,mind}.st4MReo.M12(nvecTe,nvecTo) * bcol(pvecno,1);
            dcol(pvecno, 1) = - CstTRa{1,mind}.st4MReo.M21(nvecTo,nvecTe) * acol(pvecne,1) ...
                + CstTRa{1,mind}.st4MReo.M22(nvecTo,nvecTo) * bcol(pvecno,1);

            ccol(pvecno, 1) = CstTRa{1,mind}.st4MRoe.M11(nvecTo,nvecTo) * acol(pvecno,1) ...
                - CstTRa{1,mind}.st4MRoe.M12(nvecTo,nvecTe) * bcol(pvecne,1);
            dcol(pvecne, 1) = - CstTRa{1,mind}.st4MRoe.M21(nvecTe,nvecTo) * acol(pvecno,1) ...
                + CstTRa{1,mind}.st4MRoe.M22(nvecTe,nvecTe) * bcol(pvecne,1);
        end
    end
end
stAbcdnm.pnm = transpose(pcol); % [1 x P] matrix
stAbcdnm.qnm = transpose(qcol); % [1 x P] matrix
stAbcdnm.anm = anm; % [1 x P]
stAbcdnm.bnm = bnm; % [1 x P]
if bGetR
    stAbcdnm.cnm = transpose(ccol); % [1 x P] matrix
    stAbcdnm.dnm = transpose(dcol); % [1 x P] matrix
end

end
