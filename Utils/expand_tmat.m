function [T, u, up] = expand_tmat(stT, qmax)
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

[T, q, qp] = exportTmatrix( stT, true, [], [] );
% convert indices to (u,up) for tmat.h5 convention
[u] = treams_indexing(q, qmax);
[up] = treams_indexing(qp, qmax);

end
