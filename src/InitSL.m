% InitFM    EDIFICE: Set fluid mechanics initial condition in EDIFICE
%
% [CTX] = InitFM(CTX)
%
%   Function initializes the fluid mechanics solution variables for the
%   Stokes flow problem consisting of velocities U, W, and pressure P
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20190418  Tobias Keller
%   modified  20200227   Tobias Keller


function  [CTX]  =  InitSL(CTX)


%***  prepare structures and variables
SL  =  CTX.SL;
FE  =  CTX.FE;
BC  =  CTX.BC;


%*****  initialize velocity field  ****************************************

if strcmp(BC.Type,'ConstStr')
    BC.USides(1,:)  =  BC.BGStrainr.*FE.W;
    BC.USides(2,:)  =  0;
end

dv               =  BC.USides(2) - BC.USides(1);
U                =  BC.USides(1) + repmat(((1:FE.nxU)-1)*dv/(FE.nxU-1),FE.nzU,1);
W                =  -repmat(((1:FE.nzU)'-1)*dv/FE.aspect/(FE.nzU-1),1,FE.nxU);
     
BC.UTopBot       =  zeros(2,FE.nxU);
BC.WTopBot       =  zeros(2,FE.nxU);
BC.USides        =  zeros(2,FE.nzU);
BC.WSides        =  zeros(2,FE.nzU);

BC.UTopBot(1,:)  =  U(1  ,:);
BC.UTopBot(2,:)  =  U(end,:);
BC.WTopBot(1,:)  =  W(1  ,:);
BC.WTopBot(2,:)  =  W(end,:);

BC.USides(1,:)   =  U(:,1  );
BC.USides(2,:)   =  U(:,end);
BC.WSides(1,:)   =  W(:,1  );
BC.WSides(2,:)   =  W(:,end);

SL.U             =  1e-6.*U(:);
SL.W             =  1e-6.*W(:);
    

%*****  initialize dynamic, lithostatic, total pressure fields  ***********

SL.Pref   =  SL.RhoRef.*CTX.PHYS.grav.*FE.CoordP(:,2) + CTX.BC.SurfPres;
SL.P      =  (mean(CTX.PROP.Rho(CTX.MP.Mat)) - SL.RhoRef).*CTX.PHYS.grav.*FE.CoordP(:,2);
SL.Pt     =  SL.P + SL.Pref;


%*** hand back variables and structures
SL.U  =  SL.U(:);
SL.W  =  SL.W(:);
SL.P  =  SL.P(:);

CTX.SL  =  SL;
CTX.BC  =  BC;


end

