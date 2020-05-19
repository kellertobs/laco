
% InitTop    EDIFICE: Initialise topography
%
% [CTX] = InitTopo(CTX)
%
%   Sets initial topography of FE mesh and adjusts coordinate arrays
%
%   created   20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function [CTX] = InitTopo(CTX)

FE    =  CTX.FE;
INIT  =  CTX.INIT;

X     =  FE.CoordQ2(FE.MapQ2(1,:),1);


%***** Set initial topography from input options  *************************

switch CTX.INIT.TopoMode
    
    case 'flat'
        Topo  =  FE.CoordQ2(FE.MapQ2(1,:),2);

    case 'slope'
        Topo  =  INIT.TopoHeight - X.*INIT.TopoHeight./FE.W;
        
    case 'peak'
        Topo  =  INIT.TopoHeight - INIT.TopoHeight .* exp(-(X-INIT.TopoXLoc).^2./INIT.TopoWidth^2);
        
    otherwise
        error('Unknown choice for TopoMode. Try again with "flat", "slope", or "peak"');
end

FE.CoordQ2(FE.MapQ2(1,:),2)  =  Topo;
FE                           =  UpdateFE(FE);


%*****  Adjust interior of FE mesh to topography  *************************

FE  =  RemeshFE(FE);

CTX.FE  =  FE;

end




