%*****  Run reference simulation in EPSL #####  ***************************

clear variables;
    
% choose parameter space to be tested
dt = [1,2,4,8] .* 5e2;

for i = 1:length(dt)
    
    % read in default parameter choices and runtime options
    par_LACO_DEFAULT
    
    % choose to plot results live during simulation,
    % set interval of time steps to plot and store output
    CTX.IO.LivePlot   = 'ON';
    CTX.IO.nwrite     =  floor(CTX.TIME.end/CTX.TIME.spyr/dt(i));
    
    % choose number of elements per grid axis
    N                 =  400;

    % set run identifier and create directory to store output files (see README
    % file for full list of runIDs for all preset parameter tests
    CTX.IO.RunID      =  ['LACO_','dt',num2str(abs(dt(i))),'_N',num2str(N)];
    if ~isfolder(['../out/',CTX.IO.RunID]); mkdir(['../out/',CTX.IO.RunID]); end
    
    % set half-width of magma reservoir volume source in [m]
    CTX.TIME.step     =  dt(i).*CTX.TIME.spyr;     % time step size [s]

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
