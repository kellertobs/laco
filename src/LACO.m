
% *************************************************************************
% **************************    EDIFICE    ********************************
% *************************************************************************
% 
% EDIFICE   Finite-element Stokes flow solver with visco-elastic/brittle-
%           plastic rheology, free surface, and subsurface
%           volume source
%
%   Main routine, uses options and paramters provided in the application 
%   context 'CTX' produced at runtime.
%
%   created   20170427   Tobias Keller
%   modified  20190418   Tobias Keller
%   modified  20200227   Tobias Keller

tic

fprintf('\n\n\n')
fprintf('**************************************************************\n')
fprintf('********************    start EDIFICE    *********************\n')
fprintf('**************************************************************\n\n\n')

%*****  check if parameter structure has been prepared

if ~exist('CTX','var')
    error('No parameter structure available to run EDIFICE');
end


%**************************************************************************
%*****   Initialise simulation run   **************************************
%**************************************************************************

CTX  =  StartRun(CTX);



%**************************************************************************
%*****  Time loop  ********************************************************
%**************************************************************************

fprintf('*****  start time-loop\n\n\n')

while CTX.TIME.total < CTX.TIME.end
        
    %***  update and report time information
    
    CTX.TIME.total = CTX.TIME.total + CTX.TIME.step;
    ReportTiming(CTX);

    
%**************************************************************************
%*****   Solve nonlinear system of equations   ****************************
%**************************************************************************

    disp(   '    Solve nonlinear equations');
    tic;
    
    CTX   =  NonlinearSolver(CTX);
    
    solvetime = toc;
    fprintf(1,'    solve time   =');
    fprintf(1,'  %f\n\n',solvetime);
    

%**************************************************************************
%*****   Manage I/O and live plotting   ***********************************
%**************************************************************************
    
    PrintDiagnostics(CTX);
         
    if  mod(CTX.TIME.istep,CTX.IO.nwrite) == 0
        
        
        %*****  live plotting some output variables  **********************
        
        LivePlotting(CTX);

        
        %*****  save and print output to file  ****************************
        
        SaveToFile(CTX);
    
        PlotOutput(CTX,CTX.IO.frame,1,CTX.IO.frame,'EtaVP','log',[19,23],'srf_rr','flip','print');

        CTX.IO.frame  =  CTX.IO.frame + 1;

    end
    
    CTX.TIME.istep  =  CTX.TIME.istep + 1;

    fprintf(1,'\n\n\n')
    
end


%*****   write endtime diagnostics   **************************************

finaltime = toc;

fprintf('    computational time = ')
fprintf(1,'%f\n\n',finaltime)


