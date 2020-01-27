function [M, nvec] = rvhGetFullMatrix(stMa, sMatName)
  %% rvhGetFullMatrix
% Returns full matrix from a struct stored in block-rvh form
% 
% Input:
%           - st4Ma: struct of matrix in block-rvh form
%                   with a CsMatList field and at least two fields (ending in "eo" and "oe")
%                   each is a struct with fields M11,M12,M21,M22,m,ind1,ind2
%           - sMatName (optional): string with name of matrix ,i.e. "st4MP"
%                       this name must be in the CsMatList field of st4Ma
%                       If unspecified sMAtName is first element of
%                       st4Ma.CsMatList
% Output:
%           - M: the full square matrix of size [2Nm x 2Nm] where Nm=N+1-m
%                (or N if m=0)
%                Note that Nm=length(ind1)+length(ind2)
%           - nvec: [Nm x 1] the n (or k) - values each block of the matrix
%                   corresponds to
%
% Dependency: 
% none

if nargin<2
    sMatName = stMa.CsMatList{1};
end

% Get oe struct
st4M = stMa.([sMatName, 'oe']);

ind1=st4M.ind1;
ind2=st4M.ind2;
N1=length(ind1);
N2=length(ind2);

Nm=N1+N2;
M=zeros(2*Nm);
m=st4M.m;

M(ind1,ind1) = st4M.M11;
M(ind1,Nm+ind2) = st4M.M12;
M(Nm+ind2,ind1) = st4M.M21;
M(Nm+ind2,Nm+ind2) = st4M.M22;
nvec = ( (max(m,1)): (Nm + max(m,1) - 1)) .';

% Get eo struct
% The line below was wrong in Version 1.0, it has now been corrected
st4M = stMa.([sMatName, 'eo']);

% Complete the matrix (note that in principle, ind1 is same as ind2 before and vice
% versa)
ind1=st4M.ind1;
ind2=st4M.ind2;
N1=length(ind1);
N2=length(ind2);
Nm=N1+N2;

M(ind1,ind1) = st4M.M11;
M(ind1,Nm+ind2) = st4M.M12;
M(Nm+ind2,ind1) = st4M.M21;
M(Nm+ind2,Nm+ind2) = st4M.M22;
