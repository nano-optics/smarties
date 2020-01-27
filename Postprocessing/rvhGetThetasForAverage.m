function stRtfunc = rvhGetThetasForAverage(stRtfunc)
  %% rvhGetThetasForAverage
% Modifies geometry for postprocessing by extending theta range to [0;pi]
% 
% rvhGetThetasForAverage(stRtfunc) modifies the geometry for particles with
% rvh symmetry, so that the range of theta is [0;pi]
% The input must have fields r, theta, drdt, wTheta, nNbTheta, as it will
% if it is created by sphMakeGeometry
%
% Input:
%           stRtfunc:  A geometry which only contains thetas from [0;pi/2],
%           as made by sphMakeGeometry.
%
% Output:
%           stRtfunc2: A geometry, including quadrature nodes and weights,
%           with twice as many theta's on [0;pi]
%
% Dependency: 
% none

% double up thetas to take into account of pi/2<theta<pi
stRtfunc.theta=[stRtfunc.theta;pi-stRtfunc.theta(end:-1:1)]; % [2T x 1]
stRtfunc.wTheta = [stRtfunc.wTheta/2;stRtfunc.wTheta(end:-1:1)/2];
stRtfunc.nNbTheta = 2*stRtfunc.nNbTheta;

stRtfunc.r = [stRtfunc.r;stRtfunc.r(end:-1:1)]; %[2T x 1]
stRtfunc.drdt = [stRtfunc.drdt;-stRtfunc.drdt(end:-1:1)];   %[T x 1]
