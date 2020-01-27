function stEsca = pstGetNearField(stRes, stRtfunc, stRprime, stPinmTaunm)
%% pstGetNearField
% Calculates scattered field from integral of surface fields
%
% stE = pstGetNearField(stRes, stRtfunc, stRprime) calculates the scattered field
% at points defined in stRprime from the surface-field using a
% surface-integral expression (Eq. 5.168 of Mishchenko 2002). The double integral
% is carried out using the same number of phi's as theta's on the surface
% to avoid calculating more quadrature points.
%
% Input:
%           stRes: A structure containing parameters of the simulation,
%                  with fields lambda [L x 1], epsilon2 [L x 1], epsilon1
%                  [1 x 1], stIncPar, cnm [L x P], dnm [L x P],
%                  nNmax [1 x 1], a [1 x 1], and c [1 x 1]
%           stRtfunc (or nNbPts): A structure containing the geometry, which is valid
%                     for theta on [0;pi], as from rvhGetThetasForAverage
%                     If stRtfunc.sInt='Pts', then no averages are
%                     calculated.
%                     If an integer, it is interpreted as nNbPts, and stRtfunc
%                     is recalculated with nNbPts points
%           stRprime: struct with fields r, theta, phi [1 x RP]
%                     points where field is to be calculated:
%           stPinmTaunm: (optional) The angular functions as from
%                        vshPinmTaunm, to save computing these multiple
%                        times. Else is calculated by this function
%
% Output:
%           stE: structure with the fields at requested points
%
% Dependency:
% crossComp [private], dotComp [private], pstGetResStructOneLambda,
% pstGetTangentialFields [private],rvhGetThetasForAverage, sph2cart2 [private],
% sphMakeGeometry, vshEgenThetaAllPhi, vshEthetaForPhi, vshPinmTaunm

if length(stRes.lambda)>1 % Only works for one lambda, so take lambda(1)
    disp 'pstGetNearField only works for one lambda, so take lambda(1)';
    stRes = pstGetResStructOneLambda(stRes,stRes.lambda(1));
end

% Test if stRtfunc is structure or needs to be recalculated.
if ~isstruct(stRtfunc)
    nNbPts=stRtfunc;
    stRtfunc=sphMakeGeometry(nNbPts,stRes.a,stRes.c);
    stRtfunc=rvhGetThetasForAverage(stRtfunc); % get thetas over entire range [0,pi]
end

% get theta dependence if not provided
if nargin < 4
    stPinmTaunm=vshPinmTaunm(stRes.nNmax,stRtfunc.theta); % fields are [T x nNmax]
end

epsilon2 = stRes.epsilon2;
epsilon1 = stRes.epsilon1;
lambda = stRes.lambda;

k1= 2*pi*sqrt(epsilon1)/lambda; % incident wavenumber
k2 = 2*pi*sqrt(epsilon2)/lambda; % internal wavenumber

[stNcrE, stNcrH] = pstGetTangentialFields(stRes, stRtfunc, stPinmTaunm);

% We use the same quadrature for phi but need to extend to the [-pi;pi] range
stPhi.phi = cos(stRtfunc.theta).*pi; % transform to be from -pi to pi, not 0 to pi; [T x 1]
stPhi.wPhi= stRtfunc.wTheta * pi; % [T x 1]
stPhi.nNbPhi=stRtfunc.nNbTheta;

Ex = zeros(1, length(stRprime.r));% [1 x RP]
Ey = zeros(1, length(stRprime.r));% [1 x RP]
Ez = zeros(1, length(stRprime.r));% [1 x RP]

cost=cos(stRtfunc.theta); % [T x 1]
sint=sin(stRtfunc.theta); % [T x 1]

% Observation points in cartesian coordinates
[rpx, rpy, rpz] = sph2cart2(stRprime.r, stRprime.theta, stRprime.phi); % all [1 x RP]

theta = stRtfunc.theta; % [T x 1]
wTheta = stRtfunc.wTheta;% [T x 1]
r = stRtfunc.r;% [T x 1]
drdt = stRtfunc.drdt;% [T x 1]


% Loop over phi
for ff=1:stPhi.nNbPhi
    phi = stPhi.phi(ff);
    wPhi = stPhi.wPhi(ff);
    cosf=cos(phi);
    sinf=sin(phi);

    % Calculate tangential fields at phi
    stNcrEforPhi = vshEthetaForPhi(stNcrE, phi);
    stNcrHOverH0forPhi = vshEthetaForPhi(stNcrH, phi);

    % get surface currents/dipoles in cartesian coords
    % Magnetic dipoles are m = n x E
    mx = stNcrEforPhi.Er.' .*sint*cosf ... % [T x 1]
         + stNcrEforPhi.Et.' .*cost*cosf - stNcrEforPhi.Ef.' * sinf;
    my = stNcrEforPhi.Er.' .*sint*sinf + ...
        stNcrEforPhi.Et.' .*cost*sinf + stNcrEforPhi.Ef.' *cosf;
    mz = stNcrEforPhi.Er.' .*cost - stNcrEforPhi.Et.' .*sint;
    % Electric dipoles are p = n x H/H0
    px = stNcrHOverH0forPhi.Er.' .*sint*cosf ... % [T x 1]
         + stNcrHOverH0forPhi.Et.' .*cost*cosf - stNcrHOverH0forPhi.Ef.' * sinf;
    py = stNcrHOverH0forPhi.Er.' .*sint*sinf + ...
        stNcrHOverH0forPhi.Et.' .*cost*sinf + stNcrHOverH0forPhi.Ef.' *cosf;
    pz = stNcrHOverH0forPhi.Er.' .*cost - stNcrHOverH0forPhi.Et.' .*sint;

    % source points on surface in cartesian
    [rx, ry, rz] = sph2cart2(r, theta, phi); % all [T x 1]

    % vector from points on surface to observation point
    Rx = bsxfun(@minus, rpx, rx);% [T x RP]
    Ry = bsxfun(@minus, rpy, ry);% [T x RP]
    Rz = bsxfun(@minus, rpz, rz);% [T x RP]
    R2 = abs(Rx).^2+abs(Ry).^2+abs(Rz).^2;% [T x RP]
    R = realsqrt(R2);% [T x RP]

    % R.p/abs(R)^2
    RdotpovR2 = dotComp(Rx, Ry, Rz, px, py, pz)./R2; % [T x RP]

    % (R.p)R/abs(R)^2
    eRdotpeRx = Rx.*RdotpovR2;
    eRdotpeRy = Ry.*RdotpovR2;
    eRdotpeRz = Rz.*RdotpovR2;

    % p - (R.p)R/abs(R)^2
    pminusRpRx = bsxfun(@minus, px, eRdotpeRx); % [T x RP]
    pminusRpRy = bsxfun(@minus, py, eRdotpeRy); % [T x RP]
    pminusRpRz = bsxfun(@minus, pz, eRdotpeRz); % [T x RP]

    % 3(R.p)R/abs(R)^2-p
    ThreeRpRminuspx = bsxfun(@minus, 3*eRdotpeRx, px); % [T x RP]
    ThreeRpRminuspy = bsxfun(@minus, 3*eRdotpeRy, py);
    ThreeRpRminuspz = bsxfun(@minus, 3*eRdotpeRz, pz);

    % R/abs(R) x m
    [eRcrossmx, eRcrossmy, eRcrossmz] = crossComp(Rx./R, Ry./R, Rz./R, mx, my, mz); % [T x RP]

    multiplier = (1./(k1*R).^2 - 1i./(k1*R)); % [T x RP]
    pGx=pminusRpRx + ThreeRpRminuspx.*multiplier; % [T x RP]
    pGy=pminusRpRy + ThreeRpRminuspy.*multiplier; % [T x RP]
    pGz=pminusRpRz + ThreeRpRminuspz.*multiplier; % [T x RP]

    multiplier2 = (1i*k1-1./(R)); % [T x RP]
    mdGx = eRcrossmx.*multiplier2; % [T x RP]
    mdGy = eRcrossmy.*multiplier2; % [T x RP]
    mdGz = eRcrossmz.*multiplier2; % [T x RP]

    eikrOver4piR = exp(1i*k1*R)./(4*pi*R);% [T x RP]
    integrandx = eikrOver4piR.*(k2*pGx + mdGx);%[T x RP]
    integrandy = eikrOver4piR.*(k2*pGy + mdGy);%[T x RP]
    integrandz = eikrOver4piR.*(k2*pGz + mdGz);%[T x RP]

    dS = r.*wTheta.*wPhi.*sqrt(r.^2+drdt.^2); % [T x 1]

    % Do sums as [1 x T] * [T x RP] products
    Ex = Ex + dS.' * integrandx; % [1 x RP]
    Ey = Ey + dS.' * integrandy; % [1 x RP]
    Ez = Ez + dS.' * integrandz; % [1 x RP]
end

E = realsqrt(abs(Ex).^2+abs(Ey).^2+abs(Ez).^2); % [1 x RP]
stEsca.E = E; % [1 x RP]
stEsca.Ex = Ex; % [1 x RP]
stEsca.Ey = Ey; % [1 x RP]
stEsca.Ez = Ez; % [1 x RP]
stEsca.x = rpx; % [1 x RP]
stEsca.y = rpy; % [1 x RP]
stEsca.z = rpz; % [1 x RP]

end



function [x, y, z] = sph2cart2(r, t, f)
x = r.*sin(t).*cos(f);
y = r.*sin(t).*sin(f);
z = r.*cos(t);
end

% function [r, t, f] = cart2sph2(x, y, z)
% r = sqrt(x.^2+y.^2+z.^2);
% t = acos(z./r);
% f = atan2(y, x);
% end

function [cx, cy, cz] = crossComp(x1, y1, z1, x2, y2, z2)
cx = bsxfun(@times, y1, z2) - bsxfun(@times, y2, z1);
cy = bsxfun(@times, z1, x2) - bsxfun(@times, z2, x1);
cz = bsxfun(@times, x1, y2) - bsxfun(@times, x2, y1);
end

function [dotProd] = dotComp(x1, y1, z1, x2, y2, z2)
% dotProd = x1.*conj(x2) + y1.*conj(y2) + z1.*conj(z2);
dotProd = bsxfun(@times, x1, x2) + bsxfun(@times, y1, y2) + bsxfun(@times, z1, z2);
end


function [stNcrE, stNcrH] = pstGetTangentialFields(stRes,stRtfunc,stPinmTaunm)
%% pstGetTangentialFields
% Calculate tangential electric and magnetic fields on the surface
%
% auxGetTangentialFields(stRes,stRtfunc,stPinmTaunm) calculates the
% tangential electric and magnetic fields returned as equivalent surface currents
% n x E and n x H, where n is the exterior surface normal.
%
% Input:
%           stRes: A structure containing parameters of the simulation,
%                  with fields lambda [1 x 1], epsilon2 [1 x 1], epsilon1
%                  [1 x 1], stIncPar, cnm [1 x P], dnm [1 x P],
%                  nNmax [1 x 1]
%           stRtfunc: A structure containing the geometry, which is valid
%                     for theta on [0;pi], as from rvhGetThetasForAverage
%           stPinmTaunm: (optional) The angular functions as from
%                        vshPinmTaunm, to save computing these multiple
%                        times. Else is calculated by this function
%
% Output:
%           stNcrE: A structure containing n x E in fields CErm, CEtm, CEfm
%           stNcrH: A structure containing n x E in fields CErm, CEtm, CEfm
%

rRow = transpose(stRtfunc.r); % [1 x T]
drdtRow=transpose(stRtfunc.drdt); % [1 x T]
thetaRow=transpose(stRtfunc.theta); % [1 x T]

% Get E field inside on the surface using the VSH series expansion with
% coefficients cnm and dnm and regular VSH (with 'j')
stNcrE=vshEgenThetaAllPhi(stRes.lambda,stRes.epsilon2,stRes.cnm,stRes.dnm, ...
    rRow,thetaRow,'j',stPinmTaunm);

nNmax=stRes.nNmax;

% Swap coefficients cnm and dnm to get H/H0 where H0 = k2/(i*omega*mu_0)
stNcrH=vshEgenThetaAllPhi(stRes.lambda,stRes.epsilon2,stRes.dnm, stRes.cnm,...
    rRow,thetaRow,'j',stPinmTaunm);

% Note that normal is n=n_r e_r + n_t e_t with
nr=rRow./realsqrt(rRow.^2+drdtRow.^2); % [1 x T]
nt=-drdtRow./realsqrt(rRow.^2+drdtRow.^2); % [1 x T]
%nf=0;
% tangential unit vectors are t1=-n_t e_r + n_r e_t and t2=e_phi

for mind=1:(2*nNmax+1)
    % Replaces E by n x E and H by n x H
    % n x E:
    Vr=stNcrE.CErm{1,mind};
    Vt=stNcrE.CEtm{1,mind};
    stNcrE.CErm{1,mind}=nt .* stNcrE.CEfm{1,mind};
    stNcrE.CEtm{1,mind}=- nr .* stNcrE.CEfm{1,mind};
    stNcrE.CEfm{1,mind}=nr.*Vt - nt.*Vr;
    % n x H:
    Vr=stNcrH.CErm{1,mind};
    Vt=stNcrH.CEtm{1,mind};
    stNcrH.CErm{1,mind}=nt .* stNcrH.CEfm{1,mind};
    stNcrH.CEtm{1,mind}=- nr .* stNcrH.CEfm{1,mind};
    stNcrH.CEfm{1,mind}=nr.*Vt - nt.*Vr;
end

stNcrE.theta = thetaRow;
stNcrH.theta = thetaRow;

end
