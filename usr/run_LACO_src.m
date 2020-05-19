%*****  Run reference simulation in EPSL #####  ***************************

clear variables;
    
% choose parameter space to be tested
A = [1,2,4];
[WW,HH] = meshgrid(A,A);
WW = WW(:).*250;  % parameter vector for tectonic stress number
HH = HH(:).*250;  % parameter vector for volume source number

for i = 1:length(WW)
    
    % read in default parameter choices and runtime options
    par_LACO_DEFAULT
    
    % choose to plot results live during simulation,
    % set interval of time steps to plot and store output
    CTX.IO.LivePlot   = 'ON';
    CTX.IO.nwrite     =  50;
    
    % choose number of elements per grid axis
    N                 =  400;

    % set run identifier and create directory to store output files (see README
    % file for full list of runIDs for all preset parameter tests
    CTX.IO.RunID      =  ['LACO_src_','W',num2str(abs(WW(i))),'_','H',num2str(abs(HH(i))),'_N',num2str(N)];
    if ~isfolder(['../out/',CTX.IO.RunID]); mkdir(['../out/',CTX.IO.RunID]); end
    
    % set half-width of magma reservoir volume source in [m]
    W                 =  WW(i);
    
    % set half-height of magma reservoir volume source in [m]
    H                 =  HH(i);
    
    % set horizontal semi-axis of elliptical volume source [m]
    CTX.INIT.SrcWidth       =  W;                     
    CTX.INIT.SphereRadiusX  =  [CTX.INIT.SrcWidth];

    % set vertical semi-axis of elliptical volume source [m]
    CTX.INIT.SrcHeight      =  H;          
    CTX.INIT.SphereRadiusZ  =  [CTX.INIT.SrcHeight];
    
    % get characteristic yield strength at reservoir depth
    Y                 = (CTX.PROP.Coh(2)+CTX.PROP.Frict(2)*CTX.PHYS.grav*CTX.PROP.Rho(2)*CTX.INIT.SrcZLoc/2);
    
    % set tectonic background shear rate [1/s] according to tectonic stress number, T = -1
    CTX.BC.BGStrainr  =  -1*Y/(CTX.PROP.Eta(2));
    
    % set in-/deflation rate for magma reservoir [1/s] according to source number, S = -1
    CTX.INIT.SrcAmpl  =  -1*Y/(CTX.PROP.Eta(2)*CTX.INIT.SrcREta^0.5);
    
    % set numerical mesh size
    CTX.FE.nx         =  N;
    CTX.FE.nz         =  N;
    
    % run simulation code
    run('../src/LACO.m')
    
end
