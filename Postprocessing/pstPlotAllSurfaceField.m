function pstPlotAllSurfaceField(nNbPts, stResE, lambda0, showMesh, log)
  %% pstPlotAllSurfaceField
% Makes a 3D plot of surface field intensity M=|E|^2(theta,phi)
% 
% pstPlotAllSurfaceField(nNbPts,stResE,lambda0) makes a 3D plot of
% the surface field intensity M=|E|^2(theta,phi). Note that this requires to
% compute the surface-field from scratch.
%
% Input:    - nNbPts: number of theta and phi. Total number of points for
%             plot is [ nNbPts+1 x nNbPts+1 ]
%           - stResE: structure with results from T-matrix and expansion
%             coefficients and a and c
%           - lambda0 (optional): lambda at which plot is done, by default
%                                 it uses lambda(1)
%           - showMesh (optional): [boolean] display mesh lines
%           - log (optional): [boolean] display the log_10 of intensities
%
% Dependency: 
% isOctave, pstGetResStructOneLambda, pstSurfaceField, sphMakeGeometry,
% vshEthetaForPhi, vshMakeIncidentParameters

if nargin<5
    log = false;
end

if nargin<4
    showMesh = false;
end

if nargin<3
    lambda0 = stResE.lambda(1);
end

N = stResE.nNmax;

% Incident field is defined by sIncType or stIncPar
if ~isfield(stResE,'stIncPar')
    sIncType=stResE.sIncType;
    stResE.stIncPar = vshMakeIncidentParameters(sIncType, N);
end


a=stResE.a;
c=stResE.c;
% For this plot, we need specific points given by the following function
[x,y,z] = ellipsoid(0,0,0,a,a,c,nNbPts);
% We can deduce the theta's and phi's from x,y,z
% use atan2 to avoid infinities problems
thetamat=pi/2-atan2(z,sqrt((x).^2+(y).^2)); % [E x E] where E=nNbPts+1
phimat=atan2(y,x); % [E x E]

theta=thetamat(:,1); % [E x 1] all columns of thetamat are identical
% Need to recalculate geometry for those theta's
stRtfunc2=sphMakeGeometry(0,a,c,theta);

% Get surface fields for those theta's (note that averaging will be wrong
% since we are not using quadrature nodes)
if length(stResE.lambda)==1 % only one lambda
    stEsurf2=pstSurfaceField(stResE,stRtfunc2);
else
    stRes1 = pstGetResStructOneLambda(stResE,lambda0);
    stEsurf2=pstSurfaceField(stRes1,stRtfunc2);
end

phi=phimat(2,:); % [1 x E] all rows of phimat except first and last are identical
M=zeros(length(theta),length(phi)); % [E x E]
% loop over phi
for nf=1:length(phi)
    stEforPhi=vshEthetaForPhi(stEsurf2,phi(nf));
    M(:,nf)=transpose(abs(stEforPhi.Er).^2+abs(stEforPhi.Et).^2+abs(stEforPhi.Ef).^2);
end

if(log)
figure('Name','Surface field intensity - log10(M)');
h=surf(x,y,z,log10(M),'EdgeColor','none','LineStyle','none');
else
figure('Name','Surface field intensity - M - linear scale');
h=surf(x,y,z,M,'EdgeColor','none','LineStyle','none');
end
%hpl1=gca;
%set(hpl1,'Outerposition',[0,0.5,1,0.5]);
axis tight; axis equal; axis off;
view([-37.5+180 30]);
set(gca,'Visible','off');
set(h,'LineStyle','none');

hold on;
colT=[0.05 0.05 0.05];
if showMesh % user wants to visualise mesh
  nn=1:20:(nNbPts+1);
  h2=mesh(x(nn,nn),y(nn,nn),z(nn,nn),0+0*z(nn,nn));
  set(h2,'EdgeColor',colT,'LineWidth',1.2);
end
line([-2*a,2*a],[0,0],[0,0],'LineStyle','-','Color',colT);
line([0,0],[-2*a,2*a],[0,0],'LineStyle','-','Color',colT);
hcb=colorbar('eastOutside','FontSize',12,'FontWeight','bold','YColor',colT);

if(log)
title(hcb, 'log10(|E/E0|^2)')
else
title(hcb, 'M = |E/E0|^2')
end

if ~isOctave()
  axis vis3d
end
rotate3d on
