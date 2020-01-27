function stSM = pstScatteringMatrixOA(CstTRa,lambda,Csca,nNbTheta)
%% pstScatteringMatrixOA
% Calculates scattering matrix for random orientation
% 
% stSM = pstScatteringMatrixOA(CstTRa,lambda,Csca,nNbTheta) calculates the
% quantities related to the scattering matrix of an ensemble of
% randomly-oriented scatterers.
% This function and its subroutines have been directly translated (with
% permission) from the Fortran routines written by Mishchenko and co-workers
% and publicly available at
% http://www.giss.nasa.gov/staff/mmishchenko/t_matrix.html
% Many of the variable names have been kept, along with the loop structure,
% which will make it very inefficient to run in MATLAB, but can be used for
% comparison and accuracy checks of the Fortran codes.
% The method consists in computing first a number of expansion coefficients
% alpha_n and beta_n up to some maximum n<=LMAX, from which the scattering
% matrix elements F11 etc... can be computed accurately at any angle in a second stage.
%
% Input:
%           CstTRa: A cell of structure with the T-matrix (for all m)
%           lambda: wavelength (in same unit as Csca)
%           Csca: orientation-averaged scattering cross-section
%                 Note that Csca and lambda are only needed to get the
%                 correct scaling factor for the scattering matrix
%           nNbTheta (optional, default is 2):
%                  number of angles for which the scattering
%                  matrix is computed (the angles are linearly spaced
%                  between 0 and pi). This parameter is the same as NPNA in
%                  the Fortran codes.
%
% Output:
%           stSM: structure with the following fields (which are similar to
%           the results returned by the Fortran codes
%               - LMAX [scalar]: maximum multipole order for expansion
%                   (called smax in Eqs. 4.75-4.80 of Mishchenko 2002).
%               - ALF1n,ALF2n,ALF3n,ALF4n,BET1n,BET2n: [N+1 x 1]
%                   expansion coefficients for the scattering matrix as
%                   defined in Eqs. 4.75-4.80 and obtained
%                   from Eqs. 4.109-4.114 of Mishchenko 2002.
%                   Note that those start at n=0 and are zero for n>LMAX
%               - AsymP [scalar]: Asymmetry parameter (<cos theta>)
%               - theta, thetadeg [nNBTheta x 1]: angles in Radian and
%                   Degree for the scattering matrix calculation
%               - F11, F22, F33, F44, F12, F34: scattering matrix elements
%                   as defined in Mishchenko 2002.
%               - AllAB [N+1 x 7]: An array with all the expansion
%                   coefficients as given in the output of the Fortran codes
%               - AllF [nNbTheta x 7]: An array with all the angles and
%                   scattering matrix elements as given in the output of the
%                   Fortran codes
%
% Dependency: 
% none

if nargin<4
    nNbTheta=2; % only 2 angles for scattering matrix F
end

%C *****  CALCULATION OF THE ARRAYS B1 AND B2  *****

K1=1;
K2=0;
K3=0;
K4=1;
K5=1;
K6=2;

N = length(CstTRa)-1; % assumes that all m have been calculated, i.e. absmvec=0:N

nvec=(1:N);
kvec=(1:N).';

% SSI(n) is 2n+1
% SSJ(n) is sqrt(2n+1)
FFkn = sqrt(bsxfun(@times,1./(2*kvec+1),(2*nvec+1)));

% Factor for A:
FAkn = bsxfun(@times, 1i.^kvec./sqrt(2*kvec+1), 1i.^(-nvec));

% Factorials array (as logs)
logFact=zeros(4*(N+1),1);
for ii=3:(4*(N+1))
	logFact(ii)=logFact(ii-1)+0.5*log(ii-1);
end

% Sign array (-1)^(n+1)
SSIGN=ones(4*(N+1),1);
SSIGN(2:2:end)=-1;


T1mk=zeros(2*N+1,N);
T2mk=zeros(2*N+1,N);
B1n1mn = zeros(N+1,2*N+1,N);
B2n1mn = zeros(N+1,2*N+1,N);

% Get T in normal form
CT11=cell(N+1,1);
CT12=cell(N+1,1);
CT21=cell(N+1,1);
CT22=cell(N+1,1);
for m=0:N
    nminm1=max(1,m)-1;
    CT11{m+1}=zeros(N);
    CT11{m+1}(nminm1+CstTRa{m+1}.st4MTeo.ind1,nminm1+CstTRa{m+1}.st4MTeo.ind1)=CstTRa{m+1}.st4MTeo.M11;
    CT11{m+1}(nminm1+CstTRa{m+1}.st4MToe.ind1,nminm1+CstTRa{m+1}.st4MToe.ind1)=CstTRa{m+1}.st4MToe.M11;
    CT12{m+1}=zeros(N);
    CT12{m+1}(nminm1+CstTRa{m+1}.st4MTeo.ind1,nminm1+CstTRa{m+1}.st4MTeo.ind2)=CstTRa{m+1}.st4MTeo.M12;
    CT12{m+1}(nminm1+CstTRa{m+1}.st4MToe.ind1,nminm1+CstTRa{m+1}.st4MToe.ind2)=CstTRa{m+1}.st4MToe.M12;
    CT21{m+1}=zeros(N);
    CT21{m+1}(nminm1+CstTRa{m+1}.st4MTeo.ind2,nminm1+CstTRa{m+1}.st4MTeo.ind1)=CstTRa{m+1}.st4MTeo.M21;
    CT21{m+1}(nminm1+CstTRa{m+1}.st4MToe.ind2,nminm1+CstTRa{m+1}.st4MToe.ind1)=CstTRa{m+1}.st4MToe.M21;
    CT22{m+1}=zeros(N);
    CT22{m+1}(nminm1+CstTRa{m+1}.st4MTeo.ind2,nminm1+CstTRa{m+1}.st4MTeo.ind2)=CstTRa{m+1}.st4MTeo.M22;
    CT22{m+1}(nminm1+CstTRa{m+1}.st4MToe.ind2,nminm1+CstTRa{m+1}.st4MToe.ind2)=CstTRa{m+1}.st4MToe.M22;
end

% Calculation of B1 and B2 mnk


for n=1:N % Do 100
    
    % C *****  CALCULATION OF THE ARRAYS T1 AND T2  *****
    for nn=1:N
        mmax = min(n,nn);
        for m=0:mmax
            mInd = m+1+N;
            % This is to make it eo-oe compatible
            TT1=CT11{m+1}(n,nn);
            TT2=CT12{m+1}(n,nn);
            TT3=CT21{m+1}(n,nn);
            TT4=CT22{m+1}(n,nn);
            T1=TT1+TT2;
            T2=TT3+TT4;
            % T1mnk is [2N+1 x N]
            T1mk(mInd,nn) = T1+T2;
            T2mk(mInd,nn) = T1-T2;
            if m>0 % do m<0
                T1=TT1-TT2;
                T2=TT3-TT4;
                mIndn = N + 1 - m;
                T1mk(mIndn,nn) = T1-T2;
                T2mk(mIndn,nn) = T1+T2;
            end
        end
    end
    % C  *****  END OF THE CALCULATION OF THE ARRAYS T1 AND T2  *****
 
    
% uncomment the following to return intermediate results
%     stSM.CT1R{n}=real(T1mk);
%     stSM.CT1I{n}=imag(T1mk);
%     stSM.CT2R{n}=real(T2mk);
%     stSM.CT2I{n}=imag(T2mk);
% 
%     stSM.CA1R{n}=zeros(N,N+1+n);
%     stSM.CA1I{n}=zeros(N,N+1+n);
%     stSM.CA2R{n}=zeros(N,N+1+n);
%     stSM.CA2I{n}=zeros(N,N+1+n);
    
    nn1max=N+1+n; % n1 in B goes from |m-1| to ???
    for n1Ind=1:nn1max % loop over nn1=n1+1 in B
        n1=n1Ind-1;
        
        % C  *****  CALCULATION OF THE ARRAYS A1 AND A2  *****
        
        G1 = CCG(n,n1,N,K1,K2,logFact,SSIGN);
        % G1 cg(n,all m,n1,m1=0, all nn, mm=m) is [2N+1 x N+1]
        nnmax=min(N,n1+n);
        nnmin=max(1,abs(n-n1));
        kn=n+n1Ind;
        A1k = zeros(N,1);
        A2k = zeros(N,1);
        for nn=nnmin:nnmax % loop over n' for A
            nnInd=nn+1;
            SIG=SSIGN(kn+nn); % (-1)^(kn+nn+1)
            mIndMax=min(n,nn)+ N + 1;
            AA1=0;
            AA2=0;
            % do sum on m
            for mInd=(N+1):mIndMax % loop on m from 0 to min(n,nn)
                m=mInd-(N+1);
                SSS=G1(mInd,nnInd); % cg(n,m,n1,0,nn,m)
                R1=T1mk(mInd,nn);
                R2=T2mk(mInd,nn);
                if m>0 % Do same for negative m
                    mIndn =N+1-m;
                    R1=R1+T1mk(mIndn,nn)*SIG;
                    R2=R2+T2mk(mIndn,nn)*SIG;
                end
                AA1=AA1+SSS*R1;
                AA2=AA2+SSS*R2;
            end
            A1k(nn)=AA1*FAkn(nn,n);
            A2k(nn)=AA2*FAkn(nn,n);
        end
        
        % C  *****  END OF THE CALCULATION OF THE ARRAYS A1 AND A2 ****
%         stSM.CA1R{n}(:,n1Ind)=real(A1k);
%         stSM.CA1I{n}(:,n1Ind)=imag(A1k);
%         stSM.CA2R{n}(:,n1Ind)=real(A2k);
%         stSM.CA2I{n}(:,n1Ind)=imag(A2k);
%         stSM.FA{n}=[real(FAkn(:,n)),imag(FAkn(:,n))];
        
        % K3=0, K4=1
        G2=CCG(n,n1,N,K3,K4,logFact,SSIGN);
        % G2 cg(n,all m,n1,m1=1-m, all nn, mm=1) is [2N+1 x N+1]
        mmax=min(n1+1,n);
        mmin=max(-n1+1,-n);
        for m=mmin:mmax
            mInd=m+N+1;
            BB1=0;
            BB2=0;
            for nn=nnmin:nnmax
                nnInd=nn+1;
                SSS=G2(mInd,nnInd);
                BB1=BB1+SSS*A1k(nn);
                BB2=BB2+SSS*A2k(nn);
            end
            B1n1mn(n1Ind,mInd,n)=BB1;
            B2n1mn(n1Ind,mInd,n)=BB2;
        end
    end
end

%    C  *****  END OF THE CALCULATION OF THE ARRAYS B1 AND B2 ****

D1mkn = zeros(2*N+1,N,N);
D2mkn = zeros(2*N+1,N,N);
D3mkn = zeros(2*N+1,N,N);
D4mkn = zeros(2*N+1,N,N);
D5mkn = zeros(2*N+1,N,N);

% C  *****  CALCULATION OF THE ARRAYS D1,D2,D3,D4, AND D5  *****

for n=1:N
    for nn=1:N
        mInd=min(n,nn);
        mIndMax=N+1+mInd;
        mIndMin=N+1-mInd;
        n1IndMax=mIndMax;
        for mInd=mIndMin:mIndMax
            m=mInd-(N+1);
            n1IndMin=abs(m-1)+1;
            DD1=0;
            DD2=0;
            % loop over n1 for sums of BB in D
            for n1Ind=n1IndMin:n1IndMax
                n1=n1Ind-1;
                XX=2*n1+1;
               DD1=DD1+XX*real(B1n1mn(n1Ind,mInd,n)*conj(B1n1mn(n1Ind,mInd,nn)));
               DD2=DD2+XX*real(B2n1mn(n1Ind,mInd,n)*conj(B2n1mn(n1Ind,mInd,nn)));
            end
            D1mkn(mInd,nn,n)=DD1; % Note the real part is not present in Eqs. 5.132-5.136 of the book
            D2mkn(mInd,nn,n)=DD2; % but it would give the same final results
        end
        mmax=min(n,nn+2);
        mmin=max(-n,-nn+2);
        mIndMax=N+1+mmax;
        mIndMin=N+1+mmin;
        for mInd=mIndMin:mIndMax
            m=mInd-(N+1);
            n1IndMin=abs(m-1)+1;
            DD3=0;
            DD4=0;
            DD5=0;
            mInd2=-m+2+N+1;
            for n1Ind=n1IndMin:n1IndMax
                n1=n1Ind-1;
                XX=2*n1+1;
                % Same comment as above regarding the real part
                DD3=DD3+XX*real(B1n1mn(n1Ind,mInd,n)*conj(B1n1mn(n1Ind,mInd2,nn)));
                DD4=DD4+XX*real(B2n1mn(n1Ind,mInd,n)*conj(B2n1mn(n1Ind,mInd2,nn)));
                DD5=DD5+XX*(B2n1mn(n1Ind,mInd,n)*conj(B1n1mn(n1Ind,mInd2,nn)));
            end
            D3mkn(mInd,nn,n)=DD3;
            D4mkn(mInd,nn,n)=DD4;
            D5mkn(mInd,nn,n)=DD5;
        end
    end
end
% C  *****  END OF THE CALCULATION OF THE D-ARRAYS *****

stSM.ALF1n=zeros(2*N+1,1);
stSM.ALF2n=zeros(2*N+1,1);
stSM.ALF3n=zeros(2*N+1,1);
stSM.ALF4n=zeros(2*N+1,1);
stSM.BET1n=zeros(2*N+1,1);
stSM.BET2n=zeros(2*N+1,1);

% C  *****  CALCULATION OF THE EXPANSION COEFFICIENTS *****

DK=lambda^2/(4*Csca*pi);
for lInd=1:(2*N+1) % loop on s (l here)
    ll=lInd-1;
    G1L=0; % Real
    G2L=0; % Real
    G3L=0; % Real
    G4L=0; % Real
    G5L=0; % Complex
    SL=(2*ll+1)*DK;
    for n=1:N % sum from n=1
        nnmin=max(1,abs(n-ll));
        nnmax=min(N,n+ll);
        if nnmin <= nnmax % only to this if terms in sum
            G1 = CCG(n,ll,N,K1,K2,logFact,SSIGN); % K1=1, K2=0
            if ll>=2
                G2 = CCG(n,ll,N,K5,K6,logFact,SSIGN); % K5=1, K6=2
            end
            nll=n+ll;
            for nn=nnmin:nnmax % sum over i here nn
                nnInd=nn+1;
                mmax=min(n,nn);
                mIndMin=N+1-mmax;
                mIndMax=N+1+mmax;
                SI=SSIGN(nll+nnInd); % (-1)^(n+ll+nn)
                DM1=0;
                DM2=0;
                for mInd=mIndMin:mIndMax % sum over m
                    m=mInd-(N+1);
                    if m>=0
                        SSS1=G1(mInd,nnInd);
                    else
                        mIndn = N+1-m;
                        SSS1=G1(mIndn,nnInd)*SI;
                    end
                    DM1=DM1+SSS1*D1mkn(mInd,nn,n);
                    DM2=DM2+SSS1*D2mkn(mInd,nn,n);
                end
                FFN=FFkn(nn,n); % Prefactors
                SSS=G1(N+2,nnInd)*FFN; % G1(m=1, nn)
                G1L=G1L+SSS*DM1;
                G2L=G2L+SSS*DM2*SI;
                if ll>=2
                    DM3=0;
                    DM4=0;
                    DM5=0;
                    mmax=min(n,nn+2);
                    mmin=max(-n,-nn+2);
                    mIndMax=N+1+mmax;
                    mIndMin=N+1+mmin;
                    for mInd=mIndMin:mIndMax
                        m=mInd-(N+1);
                        mIndn=N+1-m;
                        SSS1=G2(mIndn,nnInd);
                        DM3=DM3+SSS1*D3mkn(mInd,nn,n);
                        DM4=DM4+SSS1*D4mkn(mInd,nn,n);
                        DM5=DM5+SSS1*D5mkn(mInd,nn,n);
                    end
                    G5L=G5L-SSS*DM5;
                    SSS=G2(N,nnInd)*FFN;
                    G3L=G3L+SSS*DM3;
                    G4L=G4L+SSS*DM4*SI;
                end % end if l<=2
            end % Loop on nn
        end % end if nnmin <= nnmax
    end % Loop on n
    G1L=G1L*SL;
    G2L=G2L*SL;
    G3L=G3L*SL;
    G4L=G4L*SL;
    G5L=G5L*SL;
    stSM.ALF1n(lInd)=G1L+G2L;
    stSM.ALF2n(lInd)=G3L+G4L;
    stSM.ALF3n(lInd)=G3L-G4L;
    stSM.ALF4n(lInd)=G1L-G2L;
    stSM.BET1n(lInd)=real(G5L)*2;
    stSM.BET2n(lInd)=imag(G5L)*2;
    % We here keep more orders than in the Fortran codes for maximum precision
     LMAX=ll; 
     if abs(G1L)<1e-15
         break
     end
end
stSM.LMAX=LMAX;
stSM.AsymP=stSM.ALF1n(2)/3; % asymmetry parameter (<cos theta>)

% Uncomment to return intermediate results
% stSM.D1mkn=D1mkn;
% stSM.D2mkn=D2mkn;
% stSM.D3mkn=D3mkn;
% stSM.D4mkn=D4mkn;
% stSM.D5mkn=D5mkn;
% stSM.B1n1mn=B1n1mn;
% stSM.B2n1mn=B2n1mn;



% C****************************************************************
%  
% C    CALCULATION OF THE SCATTERING MATRIX FOR GIVEN EXPANSION
% C    COEFFICIENTS
%  
% C    A1,...,B2 - EXPANSION COEFFICIENTS
% C    LMAX - NUMBER OF COEFFICIENTS MINUS 1
% C    N - NUMBER OF SCATTERING ANGLES
% C        THE CORRESPONDING SCATTERING ANGLES ARE GIVEN BY
% C        180*(I-1)/(N-1) (DEGREES), WHERE I NUMBERS THE ANGLES

theta=linspace(0,pi,nNbTheta).';
thetadeg=theta*180/pi;
stSM.theta=theta;
stSM.thetadeg=thetadeg;
stSM.F11=zeros(nNbTheta,1);
stSM.F22=zeros(nNbTheta,1);
stSM.F33=zeros(nNbTheta,1);
stSM.F44=zeros(nNbTheta,1);
stSM.F12=zeros(nNbTheta,1);
stSM.F34=zeros(nNbTheta,1);


DN=1/(nNbTheta-1);
DA=DN*pi; % theta
DB=180*DN; % thetadeg
%L1MAX=LMAX+1;
TB=-DB;
TAA=-DA;
%  1003 FORMAT(' ',5X,'<',8X,'F11',8X,'F22',8X,'F33',
%      & 8X,'F44',8X,'F12',8X,'F34')
D6=sqrt(6)/4;
for I1=1:nNbTheta % Loop over scattering angles theta
    TAA=TAA+DA; % theta(I1)
    TB=TB+DB; % thetadeg(I1)
    U=cos(TAA); % u= cos(theta)
    F11=0;
    F2=0;
    F3=0;
    F44=0;
    F12=0;
    F34=0;
    P1=0;
    P2=0;
    P3=0;
    P4=0;
    PP1=1;
    PP2=(1+U)^2/4;
    PP3=(1-U)^2/4;
    PP4=D6*(U^2-1);
    for L=0:LMAX
        L1=L+1;
        F11=F11+stSM.ALF1n(L1)*PP1;
        F44=F44+stSM.ALF4n(L1)*PP1;
        if L~=LMAX
            PL1=2*L+1;
            P=(PL1*U*PP1-L*P1)/L1;
            P1=PP1;
            PP1=P;
        end
        if L>=2
            F2=F2+(stSM.ALF2n(L1)+stSM.ALF3n(L1))*PP2;
            F3=F3+(stSM.ALF2n(L1)-stSM.ALF3n(L1))*PP3;
            F12=F12+stSM.BET1n(L1)*PP4;
            F34=F34+stSM.BET2n(L1)*PP4;
            if L~=LMAX
                PL2=(L*L1)*U;
                PL3=(L1*(L^2-4));
                PL4=1/(L*(L1^2-4));
                P=(PL1*(PL2-4)*PP2-PL3*P2)*PL4;
                P2=PP2;
                PP2=P;
                P=(PL1*(PL2+4)*PP3-PL3*P3)*PL4;
                P3=PP3;
                PP3=P;
                P=(PL1*U*PP4-sqrt(L^2-4)*P4)/sqrt(L1^2-4);
                P4=PP4;
                PP4=P;
            end
        end
    end % next L
    F22=(F2+F3)/2;
    F33=(F2-F3)/2;
    stSM.F11(I1)=F11;
    stSM.F22(I1)=F22;
    stSM.F33(I1)=F33;
    stSM.F44(I1)=F44;
    stSM.F12(I1)=F12;
    stSM.F34(I1)=F34;
end % next I1

    stSM.AllAB = [(1:(2*N+1)).',stSM.ALF1n,stSM.ALF2n,stSM.ALF3n,stSM.ALF4n,stSM.BET1n,stSM.BET2n];
    stSM.AllF = [thetadeg,stSM.F11,stSM.F22,stSM.F33,stSM.F44,stSM.F12,stSM.F34];
end







function GG = CCG(n,n1,NMAX,K1,K2,logFact,SSIGN)
% C******************************************************************
% C
% C   CALCULATION OF CLEBSCH-GORDAN COEFFICIENTS
% C   (N,M:N1,M1/NN,MM)
% C   FOR GIVEN N AND N1. M1=MM-M, INDEX MM IS FOUND FROM M AS
% C   MM=M*K1+K2
% C
% C   INPUT PARAMETERS :  N,N1,NMAX,K1,K2
% C                               N.LE.NMAX
% C                               N.GE.1
% C                               N1.GE.0
% C                               N1.LE.N+NMAX
% C   OUTPUT PARAMETERS : GG(M+NPN6,NN+1) - ARRAY OF THE CORRESPONDING
% C                                       COEFFICIENTS
% C                               /M/.LE.N
% C                               /M1/=/M*(K1-1)+K2/.LE.N1
% C                               NN.LE.MIN(N+N1,NMAX)
% C                               NN.GE.MAX(/MM/,/N-N1/)
% C   IF K1=1 AND K2=0, THEN 0.LE.M.LE.N

%      REAL*8 GG(N,N+1),CD(0:2N),CU(0:2N)
GG = zeros(2*NMAX+1,NMAX+1); % for m=-NMAX:NMAX, nn=0:NMAX
CD = zeros(2*NMAX+1); % for n1=0:2NMAX
CU = zeros(2*NMAX+1); % for n1=0:2NMAX

NPN6 = NMAX+1;
if not(0<=n1  && n1<=NMAX+n && n>=1 && n<=NMAX)
    disp ('ERROR IN CCG')
    GG=0;
    return
end
NNF=min(n+n1,NMAX);
mInit=NPN6-n;
mFin=NPN6+n;
if (K1==1 && K2==0)
    mInit=NPN6;
end
% Loop over m
for mInd=mInit:mFin
    m=mInd-NPN6;
    mm=m*K1+K2;
    m1=mm-m;
    % only calculate if needed
    if abs(m1)<=n1
        NNL=max(abs(mm),abs(n-n1));
        % only calculate if needed
        if NNL<=NNF
            NNU=n+n1;
            NNM=floor((NNU+NNL)*0.5);
            if NNU==NNL
                NNM=NNL;
            end
            C=CCGIN(n,n1,m,mm,logFact,SSIGN);
            CU(NNL+1)=C;
            if NNL~=NNF
                C2=0;
                C1=C;
                % Loop on nn
                for nn=NNL+1:min(NNM,NNF)
                    A=double((nn+mm)*(nn-mm)*(n1-n+nn));
                    A=A*((n-n1+nn)*(n+n1-nn+1)*(n+n1+nn+1));
                    A=(4*nn^2)/A;
                    A=A*((2*nn+1)*(2*nn-1));
                    A=sqrt(A);
                    if nn==1
                        B=0.5*(m-m1);
                        D=0;
                    else
                        B=(2*nn*(nn-1));
                        B=((2*m-mm)*nn*(nn-1)-mm*n*(n+1)+mm*n1*(n1+1))/B;
                        D=(4*(nn-1)*(nn-1));
                        D=D*((2*nn-3)*(2*nn-1));
                        D=((nn-mm-1)*(nn+mm-1)*(n1-n+nn-1))/D;
                        D=D*((n-n1+nn-1)*(n+n1-nn+2)*(n+n1+nn));
                        D=sqrt(D);
                    end % end if nn==1
                    C=A*(B*C1-D*C2);
                    C2=C1;
                    C1=C;
                    CU(nn+1)=C;
                end % Loop on nn
                if NNF>NNM
%                    C = CGDIRECT(n,m,n1,m1,NNU,mm);
                    C = CGDIRECT(n,m,n1,m1,logFact);
                    CD(NNU+1)=C;
                    if NNU ~= (NNM+1)
                        C2=0;
                        C1=C;
                        % Decreasing loop on nn
                        for nn=(NNU-1):-1:(NNM+1)
                            A=double((nn-mm+1)*(nn+mm+1)*(n1-n+nn+1));
                            A=A*((n-n1+nn+1)*(n+n1-nn)*(n+n1+nn+2));
                            A=(4*(nn+1)^2)/A;
                            A=A*((2*nn+1)*(2*nn+3));
                            A=sqrt(A);
                            B=double(2*(nn+2)*(nn+1));
                            B=((2*m-mm)*(nn+2)*(nn+1)-mm*n*(n+1)+mm*n1*(n1+1))/B;
                            D=double(4*(nn+2)^2);
                            D=D*((2*nn+5)*(2*nn+3));
                            D=((nn+mm+2)*(nn-mm+2)*(n1-n+nn+2))/D;
                            D=D*((n-n1+nn+2)*(n+n1-nn-1)*(n+n1+nn+3));
                            D=sqrt(D);
                            C=A*(B*C1-D*C2);
                            C2=C1;
                            C1=C;
                            CD(nn+1)=C;
                        end % Decreasing loop on nn
                    end % end if NNU ~= NNM+1
                end % end if NNF>NNM
            end % end if NNL~=NNF
            
            for nn=NNL:NNF
                if nn<=NNM
                    GG(mInd,nn+1)=CU(nn+1);
                else
                    GG(mInd,nn+1)=CD(nn+1);
                end
            end % Loop on nn
        end % End if
    end % End if
end % End Loop on mInd
end


% *********************************************************************

function C = CGDIRECT (n,m,n1,m1,F)
% Note that nn and mm are unused parameters in the Fortran codes
C=F(2*n+1)+F(2*n1+1)+F(n+n1+m+m1+1)+F(n+n1-m-m1+1);
C=C-F(2*(n+n1)+1)-F(n+m+1)-F(n-m+1)-F(n1+m1+1)-F(n1-m1+1);
C=exp(C);
end


% C*********************************************************************
% C
% C   CALCULATION OF THE CLEBCSH-GORDAN COEFFICIENTS
% C   G=(N,M:N1,MM-M/NN,MM)
% C   FOR GIVEN N,n1,M,MM, WHERE NN=MAX(/MM/,/N-N1/)
% C                               /M/.LE.N
% C                               /MM-M/.LE.N1
% C                               /MM/.LE.N+N1

function G=CCGIN(n,n1,m,mm,F,SSIGN)
m1=mm-m;
if not(n>=abs(m) &&  n1>=abs(m1) && abs(mm)<=(n+n1))
    disp(' ERROR IN SUBROUTINE CCGIN');
    G=0;
    return
end
if (abs(mm)<=abs(n-n1))
    if (n1>n)
        K=n;
        n=n1;
        n1=K;
        K=m;
        m=m1;
        m1=K;
    end
    N2=n*2;
%    M2=m*2; % unused?
    N12=n1*2;
%    M12=m1*2; % unused?
    G=SSIGN(n1+m1+1) * exp(F(n+m+1)+F(n-m+1)+F(N12+1)+F(N2-N12+2)-F(N2+2) ...
        -F(n1+m1+1)-F(n1-m1+1)-F(n-n1+mm+1)-F(n-n1-mm+1));
else
    A=1;
    if mm<0
        mm=-mm;
        m=-m;
        m1=-m1;
        A=SSIGN(mm+n+n1+1);
    end
    G=A*SSIGN(n+m+1)*exp(F(2*mm+2)+F(n+n1-mm+1)+F(n+m+1)+F(n1+m1+1) ...
        -F(n+n1+mm+2)-F(n-n1+mm+1)-F(-n+n1+mm+1)-F(n-m+1)-F(n1-m1+1));
end % if (abs(MM)<=abs(n-n1))

end

