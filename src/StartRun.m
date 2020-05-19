
% StartRun  Start simulation run in EREBUS
%
%   [CTX] = StartRun(CTX)
%   Prepares data structures and sets initial conditions from input
%   options. Requires as input an application context 'CTX' produced with a
%   EREBUS parameter script
%
%   created   20161115  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20190418  Tobias Keller
%   modified  20200227   Tobias Keller


function  [CTX] = StartRun(CTX)


%*****  initialise timing / counting parameters  **********************

CTX.TIME.total  =  0;
CTX.TIME.istep  =  1;

CTX.IO.frame    =  0;

CTX.SL.it       =  0;
CTX.SL.fnorm    =  1e6;


%*****  initialise finite element grid  *******************************

CTX.FE  =  InitFE(CTX.FE);
CTX     =  InitTopo(CTX);


%*****  initialise material distribution and properties  **************

CTX      =  InitMat(CTX);
CTX.MPo  =  CTX.MP;


%*****  initialise solution variables  ********************************

CTX      =  InitSL(CTX);
CTX.SLo  =  CTX.SL;


%*****  initialise random perturbation field  *************************

CTX.MP.pert  =  InitPert(CTX.FE,CTX.INIT.PertSmooth,CTX.INIT.PertSymmetric);


%*****  update material properties and aux. fields  *******************

CTX  =  UpdateMaterialPoints(CTX);


%*****  Initialise solution vector and scaling matrix  ********************

CTX.SL.S  =  zeros(2*CTX.FE.NU+CTX.FE.NP,1);


%*****  plot initial condition and save to file  **********************

if strcmp(CTX.IO.LivePlot,'ON'); LivePlotting(CTX); end
SaveToFile(CTX);
CTX.IO.frame  =  1;

    
end

