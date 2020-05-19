
%**************************************************************************
%**********  Set Parameters for EDIFICE  **********************************
%**************************************************************************

%*****  Run options  ******************************************************

CTX.IO.SrcDir        =  '../src';                % source directory
CTX.IO.DataDir       =  '../out';                % output directory


%*****  Solver options  ***************************************************

CTX.SL.maxits        =  10;                      % max non-linear iterations
CTX.SL.atol          =  1.e-9;                   % absolute convergence tolerance (precond. residual norm)
CTX.SL.rtol          =  1.e-3;                   % relative convergence tolerance

CTX.SL.theta_it      =  0.5;                     % iterative relaxation parameter
CTX.SL.theta_dt      =  0.5;                     % time-step weighting parameter (0.5 = Cranck-Nicolson; 0 = Forward Euler; 1 = Backward Euler)

CTX.SL.Smooth        =  0.25;                    % Smoothing factor to regularise viscosity
CTX.SL.StabFact      =  1.e-24;                  % Stabilisation factor for P-coefficient matrix
CTX.SL.RhoRef        =  0;                       % reference density (= 0 for full P, != 0 for dynamic P)


%*****  finite-element mesh and Lagrangian particles options  *************

CTX.FE.ElType        = 'Q1P0';                   % finite-element type (Q1P0, Q1Q1, Q2Q1)
CTX.FE.nblock        =  6400;                    % block size for vectorised matrix assembly 

CTX.FE.nx            =  400;                     % number of elements in x-direction on FE mesh
CTX.FE.nz            =  400;                     % number of elements in z-direction on FE mesh

CTX.FE.D             =  8e3;                     % Depth of model domain [m]
CTX.FE.W             =  CTX.FE.D ...             % Width of model domain [m]
                     *  CTX.FE.nx/CTX.FE.nz;

CTX.FE.LagrMesh      =  'ON';                    % switch for Lagrangian mesh (ON/OFF), or Lagrangian free surface (SRF)


%*****  Time stepping options  ********************************************

CTX.TIME.spyr        =  3600*24*365.25;          % seconds per year

CTX.TIME.step        =  1e3 * CTX.TIME.spyr;     % time step size [s]
CTX.TIME.end         =  5e4 * CTX.TIME.spyr;     % stopping time for simulation run [s]


%*****  I/O and live plotting options  ************************************

CTX.IO.nwrite        =  10;                      % time steps between output frames
CTX.IO.LivePlot      = 'ON';                     % switch ON for live output plotting
CTX.IO.PlotStyle     = 'srf_rr';                 % 'img' = image(), 'srf' = surface(), 
                                                 % '..._lr' = reflected left,
                                                 % '..._rr' = reflected right

%*****  Physical parameters  **********************************************

CTX.PHYS.grav        =  9.81;                    % gravity [m/s2]
CTX.PHYS.RConst      =  8.314;                   % universal gas constant [J/K/mol]


%*****  set initial conditions for topography  ****************************

CTX.INIT.TopoMode    =  'peak';                  % initial topography mode 'flat', 'peak', 'slope'
CTX.INIT.TopoHeight  =  1.5e3;                   % set height for initial topography
CTX.INIT.TopoWidth   =  2.5e3;                   % set width for initial topography
CTX.INIT.TopoXLoc    =  CTX.FE.W;                % set x-location for initial topography


%*****  set initial condition for subsurface volume source  ***************

CTX.INIT.SrcAmpl     =  0e-13;                   % amplitude of volume source [1/s]
CTX.INIT.SrcREta     =  1e-4;                    % viscosity reduction factor of volume source [1]

CTX.INIT.SrcWidth    =  1e3;                     % set horizontal semiaxis of elliptical volume source [m]
CTX.INIT.SrcHeight   =  5e2^2/CTX.INIT.SrcWidth; % set vertical semiaxis of elliptical volume source [m]
CTX.INIT.SrcXLoc     =  CTX.FE.W;                % set x-location of volume source
CTX.INIT.SrcZLoc     =  5.e3;                    % set z-location of volume source


%*****  initial conditions for material types  ****************************

CTX.INIT.MatMode       =  'layer';               % initial materials mode ('layer','const')
CTX.INIT.MatThick      =  CTX.INIT.SrcZLoc;      % initial layer thickness [m]
CTX.INIT.Mat           =  [1,2];                 % material types for layers
CTX.INIT.MatLayers     =  2;                     % number of layers

CTX.INIT.AddBlock      =  0;                     % add block of material
CTX.INIT.BlockWidth    =  [0];
CTX.INIT.BlockHeight   =  [0];
CTX.INIT.BlockZLoc     =  [0];
CTX.INIT.BlockXLoc     =  [0];
CTX.INIT.BlockMat      =  [0];

CTX.INIT.AddSphere     =  1;                     % add sphere/ellipse of material
CTX.INIT.SphereRadiusX =  [CTX.INIT.SrcWidth];
CTX.INIT.SphereRadiusZ =  [CTX.INIT.SrcHeight];
CTX.INIT.SphereZLoc    =  [CTX.INIT.SrcZLoc];
CTX.INIT.SphereXLoc    =  [CTX.FE.W];
CTX.INIT.SphereMat     =  [3];

CTX.INIT.AddFault      =  0;                     % add inherited fault zone
CTX.INIT.FaultWidth    =  [0];
CTX.INIT.FaultAngle    =  [0];
CTX.INIT.FaultZLoc     =  [0];
CTX.INIT.FaultXLoc     =  [0];
CTX.INIT.FaultMat      =  [0];


%*****  Options for initial random perturbations  *************************

CTX.INIT.PertSmooth     =  20;                   % set smoothness of random noise
CTX.INIT.PertSymmetric  =  'OFF';                % switch to symmetric noise distribution
CTX.INIT.PertFrict      =  0.03;                 % amplitude of friction angle perturbation
CTX.INIT.PertCoh        =  0;                    % amplitude of cohesion perturbation


%*****  Set material properties of material types  ************************

CTX.PROP.n      =  [      1;       2;       3];  % material ID numbers
CTX.PROP.Rho    =  [   2500;    2500;    2500];  % material densities [kg/m3]
CTX.PROP.Eta    =  [  1e+22;   1e+23;   1e+22];  % material viscosities [Pas]
CTX.PROP.G      =  [ 1.e+11;  1.e+11;  1.e+11];  % material shear moduli [Pa]
CTX.PROP.Coh    =  [   1e+7;    1e+7;    1e+7];  % material cohesions [Pa]
CTX.PROP.Frict  =  [    0.3;     0.3;     0.3];  % material friction angles [deg]


%*****  Options for visco-elastic/brittle-plasti rheology  ****************

CTX.RHEO.ConstViscosity  =  'OFF';               % switch for constant viscosity mode
CTX.RHEO.Plasticity      =  'ON';                % switch for plastic failure mode
CTX.RHEO.Elasticity      =  'ON';                % switch for elastic mode

CTX.RHEO.Strainr0        =  1.e-15;              % reference strain rate [1/s]
CTX.RHEO.PowerLawExp     =  1;

CTX.RHEO.Dmg0            =  0.1;                 % reference damage strain [1]
CTX.RHEO.DmgDFrict       =  1.0;                 % damage-reduction factor for friction angle [1]
CTX.RHEO.DmgDCoh         =  1.0;                 % damage-reduction factor for cohesion [1]
CTX.RHEO.DmgHealing      =  0e-15;               % damage healing rate [1/s]

CTX.RHEO.MaxEta          =  1.e+30;              % maximum cutoff viscosity [Pas]
CTX.RHEO.MinEta          =  1.e+20;              % minimum cutoff viscosity [Pas]
CTX.RHEO.MinGdt          =  1.e+16;              % minimum cutoff elastic strength [Pas]


%***** Boundary conditions for velocity, pressure  ************************

CTX.BC.FreeSurface     =  'ON';                  % switch for free surface top boundary

CTX.BC.Type            =  'ConstStr';            % specify constant strainrate or constant velocity boundaries
CTX.BC.BGStrainr       =   0e-15;                % specify boundary background strainrate [1/s]

CTX.BC.TopBot          =  'fs';                  % top/bot boundaries: free slip 'fs'; no slip 'ns'; free '00', bottom flux 'bt'
CTX.BC.Sides           =  'fs';                  % side boundaries   : free slip 'fs'; no slip 'ns'; simple shear 'ss'; free '00'

CTX.BC.SurfPres        =  1e+5;                  % surface pressure [Pa] (1 atm = 1e+5 Pa)






 
