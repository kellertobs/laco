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


function  [CTX]  =  InitFM(CTX)


%***  prepare structures and variables
SL  =  CTX.SL;
FE  =  CTX.FE;
BC  =  CTX.BC;


%*****  initialize velocity field  ****************************************

SL.U        =  zeros(FE.NU,1);
SL.W        =  zeros(FE.NU,1);

BC.UTopBot  =  zeros(2,FE.nxU);
BC.WTopBot  =  zeros(2,FE.nxU);
BC.USides   =  zeros(2,FE.nzU);
BC.WSides   =  zeros(2,FE.nzU);
    

%*****  initialize dynamic, lithostatic, total pressure fields  ***********

SL.P      =  zeros(FE.NP,1);
SL.Plith  =  CTX.PROP.Rho(1).*CTX.PHYS.grav.*FE.CoordP(:,2) + CTX.BC.SurfPres;
SL.Pt     =  SL.P + SL.Plith;


%*** hand back variables and structures
SL.U  =  SL.U(:);
SL.W  =  SL.W(:);
SL.P  =  SL.P(:);

CTX.SL  =  SL;
CTX.BC  =  BC;

end

