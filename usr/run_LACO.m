%*****  Run reference simulation in EPSL #####  ***************************

clear variables;
    
% choose parameter space to be tested
% A = [0,-1,-2,-4,1,2,4];
A = -1;
[TT,SS] = meshgrid(A,A);
TT = TT(:);  % parameter vector for tectonic stress number
SS = SS(:);  % parameter vector for volume source number

for i = 1:length(TT)
    
    % read in default parameter choices and runtime options
    par_LACO_DEFAULT
    
    % choose to plot results live during simulation,
    % set interval of time steps to plot and store output
    CTX.IO.LivePlot   = 'ON';
    CTX.IO.nwrite     =  50;
    
    % choose number of elements per grid axis
    N                 =  800;

    % set run identifier and create directory to store output files (see README
    % file for full list of runIDs for all preset parameter tests
    if TT(i)>0; Tstr  = 'Tp'; else; Tstr = 'Tm'; end
    if SS(i)>0; Sstr  = 'Sp'; else; Sstr = 'Sm'; end
    CTX.IO.RunID      =  ['LACO_',Tstr,num2str(abs(TT(i))),'_',Sstr,num2str(abs(SS(i))),'_N',num2str(N)];
    if ~isfolder(['../out/',CTX.IO.RunID]); mkdir(['../out/',CTX.IO.RunID]); end
    
    % set tectonic stress number [nondim.]
    % sets tectoncic background stress relative to surface yield strength
    T                 =  TT(i);
    
    % set volume source number [nondim.]
    % sets stress due to volume source relative to yield strength at depth
    S                 =  SS(i);
    
    % get characteristic yield strength at reservoir depth
    Y                 = (CTX.PROP.Coh(2)+CTX.PROP.Frict(2)*CTX.PHYS.grav*CTX.PROP.Rho(2)*CTX.INIT.SrcZLoc/2);
    
    % set tectonic background shear rate [1/s] according to tectonic stress number, T
    CTX.BC.BGStrainr  =  T*Y/(CTX.PROP.Eta(2));
    
    % set in-/deflation rate for magma reservoir [1/s] according to source number, S
    CTX.INIT.SrcAmpl  =  S*Y/(CTX.PROP.Eta(2)*CTX.INIT.SrcREta^0.5);
    
    % set numerical mesh size
    CTX.FE.nx         =  N;
    CTX.FE.nz         =  N;
    
    % run simulation code
    run('../src/LACO.m')
    
end
