function stEAllPhi = vshEgenThetaAllPhi(lambda, epsilon, pnm, qnm, rt, theta, sBessel, stPinmTaunm)
%% vshEgenThetaAllPhi
% Calculates the field on a surface r(theta) using series expansions
% 
% vshEgenThetaAllPhi(lambda,epsilon,pnm,qnm,rt,theta,sBessel,stPinmTaunm)
% Low-level function to calculate VSH expansions for r(theta) and
% many lambda. The VSHs used can either be of type 1 (sBessel='j') or
% type 3 (sBessel='h1').
% Use rt(theta)=Inf (with 'h1') to obtain far-field properties.
% The fields Ecr, Ect, Esf given in the results are discussed in the
% supplementary information.
%
% Input:
% - lambda: [L x 1] wavelengths
% - epsilon: [1 x 1] or [L x 1] dielectric function where field is evaluated
% - pnm, qnm: 2 matrices [L x P] where P = nNmax(nNmax+2) containing the
%              coefficients p_{n,m} and q_{n,m} for the field expansion of
%              M^(i)_{n,m} and N^(i)_{n,m}, respectively
% - rt:     possibly row vector [1 x T] possibly zero (if sBessel='j')
%           spherical coordinate r (in nm) of points, for each
%           corresponding theta. If values are Inf, then these correspond
%           to far-field values.
% - theta:  possibly row vector [1 x T]
%           with spherical coordinate theta of points
% - sBessel: string defining the Bessel function to be used
%           sBessel='j' for regular VSH or sBessel='h1' for irregular
% - stPinTaun: structure (optional)
%              with functions of theta pi_n and tau_n
%              if omitted, then the functions are computed from scrath.
%              It is faster to pass this structure as argument if these
%              functions have already been calculated
%
% Output:
%           stEAllPhi, structure with 3 fields representing the three
%           components E_r, E_t=E_theta, E_f=E_\phi of the field.
%           Each is a cell {2N+1 x 1}, one element A_m for each m=-N .. N
%           A_m is a matrix [L x T] for the wavelength and theta dependence
%           such that A = sum_{m=-N..N} A_m exp(i m phi)
%           Explicitly, the fields are named
%        - CErm: {2nNmax+1 x [L x T]} wavelength-and-theta-dependent Er
%        - CEtm: {2nNmax+1 x [L x T]} wavelength-and-theta-dependent Etheta
%        - CEfm: {2nNmax+1 x [L x T]} wavelength-and-theta-dependent Ephi
%
% Dependency: 
% vshGetZnAll [private], vshPinmTaunm

nPmax=size(pnm,2);
nNmax= round(-1+realsqrt(1+nPmax));
nNbLambda=length(lambda);
if (size(theta,1)~=1)
    disp 'vshEgenThetaAllPhi error: theta must be row vector.';
end
if (length(rt)~=length(theta))
    if ((rt(1)~=0) && (rt(1)~=Inf))
        disp 'vshEgenThetaAllPhi error: theta and rt must be same size row vectors.';
    end
end
nNbTheta = length(theta);

n=1:nNmax; %[1 x nNmax]
% need n-dep coeff for series mu_n *sqrt(n *(n+1)) and mu_n/sqrt(n *(n+1))
muntimes=realsqrt((2*n+1).*n.*(n+1)/(4*pi));% mu_n [1 x nNmax] for Et and Ef
mundivdgen=muntimes./(n.*(n+1));

stEAllPhi.theta=theta;
stEAllPhi.roftheta=rt;
stEAllPhi.CErm = cell(1,2*nNmax+1);
stEAllPhi.CEfm = cell(1,2*nNmax+1);
stEAllPhi.CEtm = cell(1,2*nNmax+1);

if rt(1)==0 % then all rt(ii) should be zeros
    for m=[-nNmax:-2, 2:nNmax]
        stEAllPhi.CErm{m+nNmax+1} = zeros(nNbLambda,nNbTheta);
        stEAllPhi.CEtm{m+nNmax+1} = zeros(nNbLambda,nNbTheta);
        stEAllPhi.CEfm{m+nNmax+1} = zeros(nNbLambda,nNbTheta);
    end;
    % special case where r0=0
    coef1=1/realsqrt(6*pi);
    coef2=coef1/realsqrt(2); % [1 x 1]
    % Results are all [L x T] matrices
    % m=0
    stEAllPhi.CErm{1,1+nNmax}=(coef1 * qnm(:,2)) * cos(theta); % [[Lx1] x [1xT]] = [L x T]
    stEAllPhi.CEtm{1,1+nNmax}=(-coef1 * qnm(:,2)) * sin(theta); % [L x T]
    stEAllPhi.CEfm{1,1+nNmax}=zeros(nNbLambda,nNbTheta); % [L x T]
    % m=1
    stEAllPhi.CErm{1,2+nNmax}=(-coef2 * qnm(:,3)) * sin(theta); % [L x T]
    stEAllPhi.CEtm{1,2+nNmax}=(-coef2 * qnm(:,3)) * cos(theta); % [L x T]
    stEAllPhi.CEfm{1,2+nNmax}=(-1i*coef2 * qnm(:,3)) * ones(1,nNbTheta); % [L x T]
    % m=-1
    stEAllPhi.CErm{1,nNmax}=(coef2 * qnm(:,1)) * sin(theta); % [L x T]
    stEAllPhi.CEtm{1,nNmax}=(coef2 * qnm(:,1)) * cos(theta); % [L x T]
    stEAllPhi.CEfm{1,nNmax}=(-1i*coef2 * qnm(:,1)) * ones(1,nNbTheta); % [L x T]
    disp 'r0=0 in vshEgenThetaAllPhi';
    return
end

% r~=0 from here

% get Zn(rho) for radial dependence and derived functions
if ~(rt(1)==Inf)
    % Matrix product of [L x 1] by [1 x T]
    kr = (2*pi*sqrt(epsilon)./lambda) * rt; % [L x T]
    % Column vector [LT x 1] for all [L x T] arguments
    rhocol = (reshape(transpose(kr),nNbLambda*nNbTheta,1));
    stZnAllcol=vshGetZnAll(nNmax,rhocol,sBessel); % fields are [LT x nNmax]
else % Special case for radiation profile
    stZnAllcol.Z1=ones(nNbLambda*nNbTheta,nNmax);
    stZnAllcol.Z2=ones(nNbLambda*nNbTheta,nNmax);
    stZnAllcol.Z0=ones(nNbLambda*nNbTheta,nNmax);
end

% get theta dependence if not provided
if nargin < 8
    stPinmTaunm=vshPinmTaunm(nNmax,transpose(theta)); % fields are [T x nNmax]
end


% loop over lambda to avoid using 3-dimension matrices:
% At a fixed lambda, the sum over n for all theta can be carried out
% using a matrix product of the theta-and-n-dependent [T x nNmax] matrix
% by a n-dependent column vector [nNmax x 1]
% the result, a [T x 1] matrix is then transposed to a [1 x T] line

for m = -nNmax:nNmax

    if(m==0)
        nvec = 1:nNmax;
        pvec = nvec.*(nvec +1);
        pinm = zeros(nNbTheta,nNmax);
        dnm = stPinmTaunm.pn0;
    else
        nvec = abs(m):nNmax; % [1 x N2]
        pvec = nvec.*(nvec+1) + m;    % [1 x N2]
        pinm = stPinmTaunm.pinm(:,pvec);   % [T x N2];
        dnm = bsxfun(@times,pinm,transpose(sin(theta))/m); % [T x N2]
    end;
    taunm = stPinmTaunm.taunm(:,pvec); % [T x N2]

    Ersum=zeros(nNbLambda,nNbTheta);
    Etsum=zeros(nNbLambda,nNbTheta);
    Efsum=zeros(nNbLambda,nNbTheta);

    if (rt(1)==Inf)
        % for far-field radiation profile
        qnmForZ1=zeros(nNbLambda,1); % [L x 1]
        ipnmForZ0=pnm(:,pvec); % [L x N2]
        qnmForZ2=qnm(:,pvec); % [L x N2]
        n=1:nNmax;
        mundivd=mundivdgen.*((-1i).^(n+1));
    else
        qnmForZ1=qnm(:,pvec); % [L x N2]
        ipnmForZ0=1i*pnm(:,pvec); % [L x N2]
        qnmForZ2=qnm(:,pvec); % [L x N2]
        mundivd=mundivdgen;
    end

    for ll=1:nNbLambda
        indInRhocol=(1:nNbTheta)+(ll-1)*nNbTheta;
        % for Er, vecN=d_{n,1} * Z_n^1(rho) * mu_n * n *(n+1)
        vecNdep=transpose(qnmForZ1(ll,:).*muntimes(1,nvec)); % [N2 x 1]
        % Ersum=sum_n(pi_n(t) * vecN_n)
        % Do the sum as a matrix product
        Ersum(ll,:)=transpose((dnm.*stZnAllcol.Z1(indInRhocol,nvec)) * vecNdep); % [1 x T]
        % for Et and Ef
        vecNdep=transpose(ipnmForZ0(ll,:).*mundivd(1,nvec)); % [N2 x 1]
        vecNdep2=transpose(qnmForZ2(ll,:).*mundivd(1,nvec)); % [N2 x 1]
        % Do the sums as matrix products
        tmp1=(pinm.*stZnAllcol.Z0(indInRhocol,nvec)) * vecNdep;
        tmp2=(taunm.*stZnAllcol.Z2(indInRhocol,nvec)) * vecNdep2;
        Etsum(ll,:)=transpose(tmp1+tmp2); % [1 x T]

        tmp1=(taunm.*stZnAllcol.Z0(indInRhocol,nvec)) * vecNdep;
        tmp2=(pinm.*stZnAllcol.Z2(indInRhocol,nvec)) * vecNdep2;
        Efsum(ll,:)=transpose(tmp1+tmp2); % [1 x T]
    end

    stEAllPhi.CErm{1,m + nNmax+1} = (-1)^m*Ersum;
    stEAllPhi.CEtm{1,m + nNmax+1} = (-1)^m*Etsum;
    stEAllPhi.CEfm{1,m + nNmax+1} = 1i*(-1)^m*Efsum;
end

end


function stZnAll=vshGetZnAll(nNmax, rho,sBessel)
% computes the three Zn(rho) auxiliary functions for the radial
% dependence of VSHs for n=1 to nNmax
% can be used for both regular VSHs (based on j(rho)) or
% irregular VSHs (based on h1(rho)).
%
% Parameters:
% - nNmax: scalar integer
%          number of n in series
% - rho:   column vector [R x 1] (no zero components allowed, even for regular VSH, for speed optimization)
%          arguments of the VSHs
% - sBessel: string defining the Bessel function to be used
%           sBessel='j' for or sBessel='h1'
%
%
% Returns: stZnAll structure with 3 fields
%          containing matrices [R x nNmax]
% - stZnAll.Z0 is Z_n^0(rho)=z_n(rho)
% - stZnAll.Z1 is Z_n^1(rho)=z_n(rho)/rho
% - stZnAll.Z2 is Z_n^2(rho)=[rho*z_n(rho)]'/rho
%

if ~isempty(find(rho==0,1))
    disp 'Warning: rho=0 arguments not allowed in vshZnAll...'
end

n=1:nNmax;
nm1=0:nNmax;
nu=nm1+0.5;

f=zeros(length(rho),nNmax+1);

for rhoind=1:length(rho)
    f(rhoind,:)=besselj(nu,rho(rhoind));
    if ~isempty(find(f(rhoind,:)==0,1)) % beyond double precision
        disp(['Warning: Bessel (j) calculation went beyond precision in ', mfilename]);
        disp(['x=', num2str(rho(rhoind)), ' Nmax=', int2str(max(nm1))]);
    end
end


if strcmpi(sBessel,'h1')
    for rhoind=1:length(rho)
        y=bessely(nu,rho(rhoind));
        if find(isinf(y)==true,1) % beyond double precision
           disp(['Warning: Bessel (y) calculation went beyond precision in ', mfilename]);
            disp(['x=', num2str(rho(rhoind)), ' Nmax=', int2str(max(nm1))]);
        end
        f(rhoind,:)=f(rhoind,:)+1i*y;
    end
else
    if ~strcmpi(sBessel,'j')
        disp 'Error in vshGetZnAll: wrong sBessel string, assuming j...'
    end
end

% f is matrix [R x nNmax+1] of cylindrical Bessel
% Z_{n+0.5}(rho), n=0..nNmax

f=bsxfun(@times,f,sqrt((pi/2)./rho)); % [R x nNmax+1]
% f is now matrix of spherical Bessel
% z_n(rho), n=0..nNmax or equivalently z_{n-1}(rho), n=1..nNmax+1

stZnAll.Z0=f(:,2:(nNmax+1)); % [R x nNmax]
stZnAll.Z1=bsxfun(@times,stZnAll.Z0, 1./ rho); % [R x nNmax]

% Computes: Z2_n=z_{n-1} - n Z1_n
stZnAll.Z2=f(:,n) -bsxfun(@times,stZnAll.Z1,n);

end
