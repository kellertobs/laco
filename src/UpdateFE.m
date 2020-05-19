
% UpdateFE    Update properties of the finite element mesh
%
% [FE]  =  UpdateFE(FE)
%
%   Updates finite element mesh parameters after changes to the geometry 
%   due to lagrangian mesh deformation or remeshing.
%
%   created   20140806  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  [FE]  =  UpdateFE(FE)

%*****  Update global Element coordinates  ********************************

FE.CoordEl  =  PQ2El(FE.CoordQ2,FE);


%*****  Update global Q1 node coordinates  ********************************

FE.CoordQ1  =  PQ2Q1(FE.CoordQ2,FE);


%*****  Update global integration point coordinates  **********************

FE.CoordIP  =  PQ2IP(FE.CoordQ2,FE);


%*****  Update surface coordinates  ***************************************

FE.SurfEl  =  FE.CoordEl(FE.MapEl(1,:),2);
FE.SurfQ1  =  FE.CoordQ1(FE.MapQ1(1,:),2);
FE.SurfQ2  =  FE.CoordQ2(FE.MapQ2(1,:),2);


%*****  Select appropriate coordinates for FM DoFs  ***********************

switch FE.ElType
    case 'Q2Q1'
        FE.CoordU  =  FE.CoordQ2;
        FE.CoordP  =  FE.CoordQ1;
    case 'Q1Q1'
        FE.CoordU  =  FE.CoordQ1;
        FE.CoordP  =  FE.CoordQ1;
    case 'Q1P0'
        FE.CoordU  =  FE.CoordQ1;
        FE.CoordP  =  FE.CoordEl;
end


%*****  Update element deformation state  *********************************

FE.MaxDef  =  max( abs(FE.CoordQ1(FE.El2Q1(:,3),1)-FE.CoordQ1(FE.El2Q1(:,1),1)) ...
                ./ abs(FE.CoordQ1(FE.El2Q1(:,2),2)-FE.CoordQ1(FE.El2Q1(:,1),2)) );
FE.MinDef  =  min( abs(FE.CoordQ1(FE.El2Q1(:,3),1)-FE.CoordQ1(FE.El2Q1(:,1),1)) ...
                ./ abs(FE.CoordQ1(FE.El2Q1(:,2),2)-FE.CoordQ1(FE.El2Q1(:,1),2)) );
FE.AvgDef  =  mean(abs(FE.CoordQ1(FE.El2Q1(:,3),1)-FE.CoordQ1(FE.El2Q1(:,1),1)) ...
                ./ abs(FE.CoordQ1(FE.El2Q1(:,2),2)-FE.CoordQ1(FE.El2Q1(:,1),2)) );

                   
%*****  Update element volume  ********************************************

COORDX    =  reshape(FE.CoordQ2(FE.El2Q2,1),FE.NEl,FE.Q2pEl);
COORDZ    =  reshape(FE.CoordQ2(FE.El2Q2,2),FE.NEl,FE.Q2pEl);
FE.ElVol  =  zeros(FE.NEl,1);

for ip = 1:FE.IPpEl
    dNdSi     =  FE.dNdSiQ2(:,:,ip);
    Jx        =  COORDX*dNdSi';
    Jy        =  COORDZ*dNdSi';
    detJ      =  Jx(:,1).*Jy(:,2) - Jx(:,2).*Jy(:,1);
    FE.ElVol  =  FE.ElVol + detJ.*FE.Wi(ip);
end


%*****  Update mesh dimensions  *******************************************

FE.D       =  mean(FE.CoordQ1(FE.MapQ1(end,:),2));
FE.W       =  mean(FE.CoordQ1(FE.MapQ1(:,end),1))-mean(FE.CoordQ1(FE.MapQ1(:,1),1));

FE.hxQ1    =  FE.W./(FE.nxQ1-1);
FE.hzQ1    =  FE.D./(FE.nzQ1-1);
FE.hxQ2    =  FE.W./(FE.nxQ2-1);
FE.hzQ2    =  FE.D./(FE.nzQ2-1);

FE.aspect  =  FE.W/FE.D;

end