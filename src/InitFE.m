
% InitFE    Create 2D finite element mesh for use in EDIFICE
%
% [FE] = InitFE(FE)
%
%   sets  parameters related to the FE mesh and creates coordinate arrays 
%   FE.Q1Coord (#Q1nodes,2), FE.Q2Coord (#Q2nodes,2) and 
%   PAR.FE.ElCoord(#elements,2), containing coordinates of all nodes and 
%   elements on the FE mesh
%
%   created   20140729  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20190418  Tobias Keller
%   modified  20200227   Tobias Keller


function [FE] = InitFE(FE)

%***** Setup geometry of elements and nodes  ******************************

FE.aspect  =  FE.nx/FE.nz;
FE.nxEl    =  FE.nx;  FE = rmfield(FE,'nx');
FE.nzEl    =  FE.nz;  FE = rmfield(FE,'nz');

FE.nxQ1    =  FE.nxEl*1+1;
FE.nzQ1    =  FE.nzEl*1+1;

FE.nxQ2    =  FE.nxEl*2+1;
FE.nzQ2    =  FE.nzEl*2+1;

FE.nxIP    =  FE.nxEl*3;
FE.nzIP    =  FE.nzEl*3;

FE.NEl     =  FE.nxEl*FE.nzEl;
FE.NQ1     =  FE.nxQ1*FE.nzQ1;
FE.NQ2     =  FE.nxQ2*FE.nzQ2;
FE.NIP     =  FE.nxIP*FE.nzIP;

FE.Q1pEl   =  4;
FE.Q2pEl   =  9;
FE.IPpEl   =  9;


%*****  Initialize int points and shape functions on local element  *******

[FE.Si,FE.Wi]         =  IntCoords(FE);
[FE.NiP0,FE.dNdSiP0]  =  ShapeFuncts(FE.Si,'P0');
[FE.NiQ1,FE.dNdSiQ1]  =  ShapeFuncts(FE.Si,'Q1');
[FE.NiQ2,FE.dNdSiQ2]  =  ShapeFuncts(FE.Si,'Q2');


%*****  Initialize DoFs on specified element type  ************************

switch FE.ElType
    
    case 'Q2Q1'
        FE.nxU  =  FE.nxQ2;
        FE.nzU  =  FE.nzQ2;
        FE.nxP  =  FE.nxQ1;
        FE.nzP  =  FE.nzQ1;
        FE.NU   =  FE.NQ2;
        FE.NP   =  FE.NQ1;
        FE.UpEl =  FE.Q2pEl;
        FE.PpEl =  FE.Q1pEl;
        
        FE.NiU     =  FE.NiQ2;
        FE.dNdSiU  =  FE.dNdSiQ2;
        FE.NiP     =  FE.NiQ1;
        FE.dNdSiP  =  FE.dNdSiQ1;
        
    case 'Q1Q1'
        FE.nxU  =  FE.nxQ1;
        FE.nzU  =  FE.nzQ1;
        FE.nxP  =  FE.nxQ1;
        FE.nzP  =  FE.nzQ1;
        FE.NU   =  FE.NQ1;
        FE.NP   =  FE.NQ1;
        FE.UpEl =  FE.Q1pEl;
        FE.PpEl =  FE.Q1pEl;
        
        FE.NiU     =  FE.NiQ1;
        FE.dNdSiU  =  FE.dNdSiQ1;
        FE.NiP     =  FE.NiQ1;
        FE.dNdSiP  =  FE.dNdSiQ1;
        
    case 'Q1P0'
        FE.nxU  =  FE.nxQ1;
        FE.nzU  =  FE.nzQ1;
        FE.nxP  =  FE.nxEl;
        FE.nzP  =  FE.nzEl;
        FE.NU   =  FE.NQ1;
        FE.NP   =  FE.NEl;
        FE.UpEl =  FE.Q1pEl;
        FE.PpEl =  1;
        
        FE.NiU     =  FE.NiQ1;
        FE.dNdSiU  =  FE.dNdSiQ1;
        FE.NiP     =  FE.NiP0;
        FE.dNdSiP  =  FE.dNdSiP0;
end

FE.DOFU    =  (1:2:2*FE.NU ).';
FE.DOFW    =  (2:2:2*FE.NU ).';
FE.DOFP    =  (1:1:  FE.NP ).'+2*FE.NU;


%*****  Initialize global mapping arrays  *********************************

FE  =  FEMapping(FE);


%*****  Initialize global Q2 coordinates  *********************************

if ~isfield(FE,'W'); FE.W = FE.D*FE.nxQ1./FE.nzQ1; end   % get width from numerical aspect ratio 
                                                         % unless specified otherwise in input parameters
    
hx     =  FE.W/FE.nxEl/2;
hz     =  FE.D/FE.nzEl/2;
x      =  0:hx:FE.W;
z      =  0:hz:FE.D;
[X,Z]  =  meshgrid(x,z);

Coord         =  zeros(FE.nzQ2,FE.nxQ2,2);
Coord(:,:,1)  =  X(:,:);
Coord(:,:,2)  =  Z(:,:);

FE.CoordQ2    =  reshape(Coord,FE.NQ2,2);
    
[FE]     =  UpdateFE(FE);


end




