function InitPath()
% Initialization function to add all relevant directories to Matlab path
% Must be run once after starting Matlab and before running other
% scripts/functions

% get current directory
sTmp=what;
% add it and all subdirectory to matlab path
addpath(genpath(sTmp.path));

sMsg=sprintf('Directory : %s and all subdirectories added to Matlab path',sTmp.path);
disp(sMsg);
