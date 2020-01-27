function quadTable = updateGLquadrature(Nt, keepOld)
%% updateGLquadrature
% Calculates and stores points and weights for integral quadrature. This
% calculates for Nt not already stored in Utils/quadTable.mat, and returns
% the full structure of values. If the file doesn't exist, then it
% populates it with the values Nt. It also saves the values back into
% Utils/quadTable.mat.
%
%   Input:
% - Nt: vector containing the nodes to be in quadTable
% - keepOld: boolean, default true to add to existing file, false to overwrite.
%
%	Output:
%	A structure with fields 
% - Nt: vector containing the number of nodes
% - values: cell array containing matrices of [theta, wTheta]
%
% Dependency: 
% auxInitLegendreQuad

if nargin < 2
    keepOld = true;
end

quadTableFile = [fileparts(mfilename('fullpath')), '/quadTable.mat'];
    
if(keepOld && exist(quadTableFile, 'file') == 2)
    load(quadTableFile)
    assert(exist('quadTable', 'var')==1, 'The file quadTable.mat should contain a structure quadTable')
    %we only want to add in entries that aren't already there
    Nt = setdiff(Nt, quadTable.Nt);
    quadTable.Nt = [quadTable.Nt, Nt];%append new entries to old ones
    indexOffset = length(quadTable.values);
    oldValues = quadTable.values;
    nbNt = length(Nt);
    quadTable.values = cell(1, nbNt);
    %copy old values to new cell
    for(ii = 1:length(oldValues))
        quadTable.values{ii} = oldValues{ii};
    end
else
    quadTable = struct;
    quadTable.Nt = Nt;
    nbNt = length(Nt);
    quadTable.values = cell(1,nbNt);
    indexOffset = 0;
end


for(iNt = 1:nbNt)
    [xi,wi]=auxInitLegendreQuad(Nt(iNt));
    quadTable.values{iNt+indexOffset} = [acos(xi), wi];
end

save(quadTableFile, 'quadTable')
end
