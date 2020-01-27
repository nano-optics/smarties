function [] = storeGLquadrature()
%% storeGLquadrature
% Calculates and stores points and weights for integral quadrature.
% This calculates Nt in steps of 5 from 50 to 505, then in steps of 100 from
% 600 to 2000, as well as 5 above each of those values (605, 705 etc).
%
%	Output:
%	A structure with fields
% - Nt: vector containing the number of nodes
% - values: cell array containing matrices of [theta, wTheta]
%
% Dependency:
% updateGLquadrature

Nt1 = 50:5:505;
% from 600 to 2000
% generate pairs N, N+5 by steps of 100
% 50 55 100 105 etc. til 2000 2005
Nt2 = 600:100:2000;
Nt3 = reshape([Nt2; Nt2+5], 1, 2*size(Nt2, 2));
Nt = [Nt1, Nt3];

updateGLquadrature(Nt, false);%write new file without keeping added values
end
