function [ST] = sparseTmatrix(stT)
%% sparseTmatrix
% Reshaping T-matrix as complete sparse matrix
%
% The rows and columns follow the convention of Nieminen et al,
% from fastest to slowest: m [-n:n], n [1:Nmax], s [M, N]
%
% Dependency:
% none

mMax = length(stT);
Nmax = mMax - 1;
block_size = Nmax*(Nmax+2);

% grab the values, add negative m's, and reshape to long format
T = exportTmatrix(stT, true); 

% vecs  = T(:,3);
% vecsp = T(:,4);
% vecn  = T(:,5);
% vecnp = T(:,6);
% vecm  = T(:,7);
% vecmp = T(:,8);
% values = T(:,9) + 1i * T(:,10);
vecs  = T(:,1);
vecsp = T(:,2);
vecn  = T(:,3);
vecnp = T(:,4);
vecm  = T(:,5);
vecmp = T(:,6);
values = T(:,7) + 1i * T(:,8);

% now create the sparse matrix

% combined index for the whole matrix (2x2 blocks)
vecp =  (vecs  - 1) * block_size + p_index(vecn, vecm);
vecpp = (vecsp - 1) * block_size + p_index(vecnp,vecmp);
%[p_index(vecn, vecm) vecp p_index(vecnp,vecmp) vecpp]

% fprintf('%2d %2d %g\n',  [vecp vecpp values].');
% note: if wanting to visualise, be aware that Matlab removes 0s in sparse
% a trick is to add 1 to check them
%ST = sparse(vecp,vecpp,1+values,2*block_size,2*block_size);
ST = sparse(vecp, vecpp, values, 2*block_size, 2*block_size);


end


function [out1] = p_index(in1,in2)
% (n,m) -> p = n * (n+1) + m
out1 = in1 .* (in1 + 1) + in2;
end


