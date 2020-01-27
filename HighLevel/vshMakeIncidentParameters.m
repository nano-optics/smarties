function stIncPar = vshMakeIncidentParameters(sIncType, nMax, thetap, phip, alphap)
%% vshMakeIncidentParameters
% Create struct with parameters of the incident electric field
% 
% vshMakeIncidentParameters(sIncType, nMax, thetap, phip, alphap) Creates a
%		 structure containing the parameters of the incident electric
%		 field for plane wave excitation.
%        The parameters can be specified using EITHER a string
%		 describing it, OR by specifying the angles of the wavevector, and
%		 the electric field.
%
%	Input:
%		sIncType  One of 'KxEz', 'KxEy', 'KyEx', 'KyEz', 'KzEx', 'KzEy'
%          where the first x,y,z refers to the wavevector direction,
%         and the second x,y,z (necessarily different to the first one)
%         is the direction of the linearly-polarised electric field. 
%         If the string does not fit this format, then the angles
%         thetap, phip, alpha, must be specified.
%
%		nMax	  The maximum order used in the expansions, required to
%				  determine the values of m required in the calculation.
%
%		thetap	  (optional) This, along with phip, defines the direction
%				  of the incident wavevector k.
%
%		phip	  (optional) This, along with thetap, defines the direction
%				  of the incident wavevector k.
%
%		alphap	  (optional) The defines the orientation of the electric
%				  field , in the plane orthogonal to wavevector k.
%	Output:
%       stIncPar: structure with thetap, phip, and alphap
%
% Dependency: 
% none

switch(sIncType)
	case 'KzEx'
		stIncPar.type='KzEx';
		stIncPar.thetap=0;
		stIncPar.phip=0;
		stIncPar.alphap=0;
		absmvec=1:1;
	case 'KxEz'
		stIncPar.type='KxEz';
		stIncPar.thetap=pi/2;
		stIncPar.phip=0;
		stIncPar.alphap=pi;
		absmvec=0:nMax;
	case 'KxEy'
		stIncPar.type='KxEy';
		stIncPar.thetap=pi/2;
		stIncPar.phip=0;
		stIncPar.alphap=pi/2;
		absmvec=0:nMax;
	case 'KzEy'
		stIncPar.type='KzEy';
		stIncPar.thetap=0;
		stIncPar.phip=0;
		stIncPar.alphap=pi/2;
		absmvec=1:1;
	case 'KyEz'
		stIncPar.type='KyEz';
		stIncPar.thetap=pi/2;
		stIncPar.phip=pi/2;
		stIncPar.alphap=pi;
		absmvec=0:nMax;
	case 'KyEx'
		stIncPar.type='KyEx';
		stIncPar.thetap=pi/2;
		stIncPar.phip=pi/2;
		stIncPar.alphap=-pi/2;
		absmvec=0:nMax;
	otherwise
		if (nargin~=5)
			disp 'vshMakeIncidentParameters error: sIncType not recognized. thetap, phip, alphap missing..., using KzEx as default';
			stIncPar.type='KzEx';
			stIncPar.thetap=0;
			stIncPar.phip=0;
			stIncPar.alphap=0;
			absmvec=1:1;
		else
			stIncPar.type='general';
			stIncPar.thetap=thetap;
			stIncPar.phip=phip;
			stIncPar.alphap=alphap;
		    absmvec=0:nMax;
		end
end

stIncPar.absmvec = absmvec;
