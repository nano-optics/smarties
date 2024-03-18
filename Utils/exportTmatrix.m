function [T, q, qp] = exportTmatrix( stT, complete, out, format )
%% exportTmatrix
% Reshaping to long format and exporting T-matrix entries to a text file
%
% PARAMETERS:
% - stT: structure containing T-matrix elements, as returned by the program
% - complete: [logical] compute negative m's
% - out: optional output filename
% - format: string format for text output
%
% RETURNS: T-matrix elements and indices consolidated in a single matrix
% The output format consists of 8 columns:
% s sp m mp n np Tr Ti
% whereby
% * m  1st m-index
% * mp 2nd m-index (identical, due to rotational symmetry)
% * n  1st n-index
% * np 2nd n-index
% * s  1st block index (electric/magnetic)
% * sp 2nd block index (electric/magnetic)
% * Tr real(T_sspmmpnnp)
% * Ti imag(T_sspmmpnnp)
%
% Dependency:
% none

if(nargin < 2)
    complete = true;
end
if(nargin < 3)
    out = [];
end
if(nargin < 4)
%     format = '%d %d %d %d %d %d % d % d %.15g %.15g\n';
    format = '%d %d %d %d % d % d %.15g %.15g\n';
end

mMax = length(stT);

%% extract T-matrix values for each m and list with corresponding indices
% loop over positive m-values [following storage as independent cells]
tmp = cell(1,mMax); % temporarily store in a cell array
for i_m =  1:mMax

    % get all values for current m
    [M, nvec] = combine_oeeo(stT{i_m});
    % create array of indices
% seems to be incorrect order (transposed)
%     [np,sp, n,s] = ndgrid(nvec,1:2, nvec,1:2);
    [n,s, np,sp] = ndgrid(nvec,1:2, nvec,1:2);
    % since m-cell storage was
    %
    %  M11(n,np) | M12(n,np)
    %  --------- + ---------
    %  M21(n,np) | M22(n,np)
    %
    % to grab values by columns, want to vary n, then s, then np, then sp
    % want to vary
    % reshape as vectors
    vect = M(:);
    vecn = n(:); vecnp = np(:);
    vecs = s(:); vecsp = sp(:);
    % m = mp is constant
    vecm = 0*vecs + (i_m - 1);
    vecmp = vecm; % axisym
    % strip analytical zeros due to eo-oe symmetry
    ids = ~isnan(vect);
    tmp{i_m} = [vecs(ids) vecsp(ids) vecn(ids) vecnp(ids) vecm(ids) vecmp(ids) real(vect(ids)) imag(vect(ids))];

end

% combine all values for positive ms
% now have an array of the form
% vecs vecsp vecm vecmp vecn vecnp treal timag
T = cell2mat(tmp.');

% strip remaining analytic zeroes for m=0
ids = all(T(:,[7 8]) ~= 0, 2);
T = T(ids,:);

%% add values for negative m's if requested
if complete

    ids = T(:,5) ~= 0; % don't duplicate m=0 cases

    vecm = -T(ids,5);
    vecmp = vecm; % axisym
    % duplicate n and s indices (for m!=0)
    vecs =  T(ids,1);
    vecsp = T(ids,2);
    vecn =  T(ids,3);
    vecnp = T(ids,4);

    % T_{-m} = T_{m} if s=sp, -T_{m} otherwise
    sgn = (-1).^(T(ids,1) + T(ids,2));
    vectr = sgn .* T(ids,7);
    vecti = sgn .* T(ids,8);

    % update T to include negative m's
    T2 = [vecs vecsp vecn vecnp vecm vecmp vectr vecti];
    T = [T ; T2];

end


%% reorder rows by increasing s sp n np m mp (slow to fast)
% i.e. start with T11 block
% start with n=1, np=1
% vary m=mp=0:n
% then n=1, np=2, ...
% Note: careful with sortrows when T contains complex numbers,
% it doesn't treat negative m's like we'd want...
% so now taking real part and imag separately...

% T = sortrows(T, [6 5 4 3 2 1]);
% T = sortrows(T, 1:6);

T = sortrows(T, 1:6);
% [vecs vecsp vecn vecnp vecm vecmp vectr vecti];
% vecp = p_index(T(:,3), T(:,5)); % n, m -> p 
% vecpp = p_index(T(:,4), T(:,6));
% T = [vecp vecpp T];

% p = p_index(T(:,3), T(:,5));
% pp = p_index(T(:,4), T(:,6));

Lmax = max(T(:,3));
 
q = lms_index(T(:,3), T(:,5),T(:,1),Lmax);
qp = lms_index(T(:,4), T(:,6),T(:,2),Lmax);


% write to a file
if(~isempty(out))
    if strcmp(out, 'stdout')
        fileID = 1;
    else
        fileID = fopen(out, 'w');
    end
    fprintf(fileID, '%d elements of T-matrix\n', size(T, 1));
    fprintf(fileID, 's sp n np m mp Tr Ti \n');
    fprintf(fileID, format, T.');
    if ~strcmp(out, 'stdout')
        fclose(fileID);
    end
end



end


function [M, nvec] = combine_oeeo(stMa, fieldname)
%% modified from rvhGetFullMatrix to have NaNs
% Returns full matrix from a struct stored in block-rvh form
%
% Input:
%           - st4Ma: struct of matrix in block-rvh form
%                   with a CsMatList field and at least two fields (ending in "eo" and "oe")
%                   each is a struct with fields M11,M12,M21,M22,m,ind1,ind2
% Output:
%           - M: the full square matrix of size [2Nm x 2Nm] where Nm=N+1-m
%                (or N if m=0)
%                Note that Nm=length(ind1)+length(ind2)
%           - nvec: [Nm x 1] the n (or k) - values each block of the matrix
%                   corresponds to
%
% Dependency:
% none

if(nargin < 2)
    fieldname = 'st4MT';
end


% Get oe struct
% interpolation of field name
st4M = stMa.([fieldname 'oe']);
% st4M = stMa.st4MToe;

ind1=st4M.ind1;
ind2=st4M.ind2;
N1=length(ind1);
N2=length(ind2);

Nm=N1+N2;
M=NaN(2*Nm); % fill with NAs to keep track of analytical zeros
m=st4M.m;

M(ind1,ind1) = st4M.M11;
M(ind1,Nm+ind2) = st4M.M12;
M(Nm+ind2,ind1) = st4M.M21;
M(Nm+ind2,Nm+ind2) = st4M.M22;
nvec = ( (max(m,1)): (Nm + max(m,1) - 1)) .';

% Get eo struct
st4M = stMa.([fieldname 'eo']);
% st4M = stMa.st4MTeo;

% Complete the matrix
% (note that in principle, ind1 is same as ind2 before and vice versa)
ind1=st4M.ind1;
ind2=st4M.ind2;
N1=length(ind1);
N2=length(ind2);
Nm=N1+N2;

M(ind1,ind1) = st4M.M11;
M(ind1,Nm+ind2) = st4M.M12;
M(Nm+ind2,ind1) = st4M.M21;
M(Nm+ind2,Nm+ind2) = st4M.M22;

end


function [p] = p_index(l,m)
% (n,m) -> p = n * (n+1) + m
p = l .* (l + 1) + m;
end

function [q] = lms_index(l,m,s,lmax)

% (l,m) -> p = l * (l+1) + m
p = p_index(l,m);
pmax = lmax*(lmax+1)+lmax;
% 𝑙(𝑞,𝑝)=(𝑞−1)𝑝𝑚𝑎𝑥+𝑝
q = (s-1) * pmax + p;
end
