% InitMat    EDIFICE: Set initial distribution of material types
%
% [CTX] = InitMat(CTX)
%
%   Function initializes the distribution of material types with distinct
%   properties throughout domain
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  [CTX]  =  InitMat(CTX)


%***  prepare structures and variables
FE    =  CTX.FE;
INIT  =  CTX.INIT;
PROP  =  CTX.PROP;
COORD =  CTX.FE.CoordIP;
        
        
%*****  initialize material type on material points  **********************
MP.Mat  =  zeros(FE.NIP,1);


%***  constant composition
if     strcmp(INIT.MatMode(1:5),'const')
    
    MP.Mat(:)  =  INIT.Mat(1);

    
%***  two-layer model (w/ initial perturbation)
elseif strcmp(INIT.MatMode(1:5),'layer')
     
    ind          =  COORD(:,2) >= INIT.MatThick;
    MP.Mat(ind)  =  INIT.Mat(1);
    ind          =  COORD(:,2) <  INIT.MatThick;
    MP.Mat(ind)  =  INIT.Mat(2);
    
%***  three-layer model upper/lower crust over lithosphere   
elseif strcmp(INIT.MatMode(1:5),'crust')
    
    ind          =  COORD(:,2) >= INIT.Moho;
    MP.Mat(ind)  =  INIT.Mat(1);
    ind          =  COORD(:,2) <  INIT.Moho & COORD(:,2) <= INIT.Conrad;
    MP.Mat(ind)  =  INIT.Mat(2);
    ind          =  COORD(:,2) <  INIT.Conrad;
    MP.Mat(ind)  =  INIT.Mat(3);
    
%***  multi-layer stack over detachment layer
elseif strcmp(INIT.MatMode(1:5),'multi')

    h0  =  INIT.MatThick;
    n   =  INIT.MatLayers;
    dh  =  h0/n;
    
    ind          =  COORD(:,2) >= h0;
    MP.Mat(ind)  =  INIT.Mat(1);
    
    for i = 1:2:n
        ind          =  COORD(:,2) < h0-(i-1)*dh & COORD(:,2) >= h0- i   *dh;
        MP.Mat(ind)  =  INIT.Mat(2);
        ind          =  COORD(:,2) < h0- i   *dh & COORD(:,2) >= h0-(i+1)*dh;
        MP.Mat(ind)  =  INIT.Mat(3);
    end
    
end


%***  add block of specified composition
if INIT.AddBlock > 0
    for n = 1:INIT.AddBlock
        ind          =  abs(COORD(:,1)-INIT.BlockXLoc(n)) <= INIT.BlockWidth( n)/2 & ...
                        abs(COORD(:,2)-INIT.BlockZLoc(n)) <= INIT.BlockHeight(n)/2 ;
        MP.Mat(ind)  =  INIT.BlockMat(n);
    end
end
    

%***  add sphere of specified composition
if INIT.AddSphere > 0
    for n = 1:INIT.AddSphere
        ind           =  (COORD(:,1)-INIT.SphereXLoc(n)).^2./INIT.SphereRadiusX(n)^2  ...
                      +  (COORD(:,2)-INIT.SphereZLoc(n)).^2./INIT.SphereRadiusZ(n)^2 <= 1;
        MP.Mat(ind)  =  INIT.SphereMat(n);
    end
end


%***  add pre-existing fault of specified composition
if INIT.AddFault > 0
    for n = 1:INIT.AddFault
        ind  =       COORD(:,2) >= INIT.FaultDepth(n) & ...
                abs((COORD(:,1)-INIT.FaultPos(n))                         ...
             +  tand(INIT.FaultAngle(n))*((GRID.EL.D-COORD(:,2))          ...
             -  0.5.*(GRID.EL.D-INIT.FaultDepth(n)))) <= 0.5*INIT.FaultWidth(n)/cosd(INIT.FaultAngle(n));
        MP.Mat(ind)  =  INIT.FaultMat(n);
    end
end


%***  initialise subsurface volume source field
a  =  exp(-(COORD(:,1)-INIT.SrcXLoc).^2./(INIT.SrcWidth ).^2) ...
   .* exp(-(COORD(:,2)-INIT.SrcZLoc).^2./(INIT.SrcHeight).^2);
MP.SrcdVdt  =  a .* INIT.SrcAmpl;
MP.SrcREta  =  INIT.SrcREta .^ a;


%*****  initialize material properties and auxiliary fields  **************
MP.Plith      =  mean(CTX.PROP.Rho(MP.Mat)).*CTX.PHYS.grav.*FE.CoordIP(:,2) + CTX.BC.SurfPres;
MP.Eta        =  PROP.Eta(MP.Mat);
MP.EtaVP      =  PROP.Eta(MP.Mat);
MP.EtaVEP     =  PROP.Eta(MP.Mat);
MP.Chi        =  zeros(FE.NIP,1);
MP.Gdt        =  PROP.G(MP.Mat).*CTX.TIME.step;
MP.Edot       =  zeros(FE.NIP,3);
MP.Edot(:,1)  =  ones(FE.NIP,1).*CTX.RHEO.Strainr0;
MP.Edot(:,2)  = -ones(FE.NIP,1).*CTX.RHEO.Strainr0;
MP.Edot(:,3)  =  1e-24;
MP.Tau        =  zeros(FE.NIP,3);
MP.Taur       =  MP.Tau;
MP.EII        =  zeros(FE.NIP,4);
MP.EII(:,1)   =  SecondInv(MP.Edot);
MP.EII(:,2)   =  SecondInv(MP.Edot);
MP.TII        =  SecondInv(MP.Tau);
MP.TIIr       =  SecondInv(MP.Taur);
MP.YieldStr   =  PROP.Coh(MP.Mat);
MP.Dmg        =  zeros(FE.NIP,1);
MP.Rho        =  PROP.Rho(MP.Mat);

if ~isfield(CTX.SL,'RhoRef')
    if strcmp(CTX.BC.FreeSurface,'OFF')
        CTX.SL.RhoRef  =  mean(CTX.MP.Rho);
    else
        CTX.SL.RhoRef  =  0;
    end
end

%***  return structures
CTX.MP  =  MP;


end

