
% RemeshFE    EDIFICE: Remesh distorted FE mesh
%
% [FE] = RemeshFE(FE)
%
%   Remeshes finite-element mesh after deformation in full or surface
%   Lagrangian mesh mode. Surface topography is preserved on a regularly
%   spaced mesh in horizontal direction. Distortion from matching
%   topography is faded out with depth into the interior of the mesh.
%
%   created   20170427  Tobias Keller
%   modified  20200227   Tobias Keller


function  [FE]  =  RemeshFE(FE)


Xo  =  FE.CoordQ2;


%*****  Get regular reference mesh of appropriate dimensions  *************

if ~strcmp(FE.LagrMesh,'ON')
    
    %***  Re-initialise mesh of original dimensions
    hx            =  FE.W/(FE.nxQ2-1);
    hz            =  FE.D/(FE.nzQ2-1);
    x             =  0:hx:FE.W;
    z             =  0:hz:FE.D;
    
    [X,Z]         =  meshgrid(x,z);
    
    Coord         =  zeros(FE.nzQ2,FE.nxQ2,2);
    Coord(:,:,1)  =  X(:,:);
    Coord(:,:,2)  =  Z(:,:);
    
else
    
    %***  Initialise new box mesh according to mean boundary deformation
    W     =  Xo(FE.MapQ2(:,FE.nxQ2),1)-Xo(FE.MapQ2(:,1),1);
    D     =  Xo(FE.MapQ2(FE.nzQ2,:),2)-Xo(FE.MapQ2(1,:),2);
    hx    =  mean(W)/(FE.nxQ2-1);
    hz    =  mean(D)/(FE.nzQ2-1);
    minX  =  mean(Xo(FE.MapQ2(:,1),1));
    minZ  =  mean(Xo(FE.MapQ2(1,:),2));
    
    Coord(:,:,1)  =  repmat((minX:hx:mean(W)+minX) ,FE.nzQ2,1);
    Coord(:,:,2)  =  repmat((minZ:hz:mean(D)+minZ)',1,FE.nxQ2);
    
end


%*****  Interpolate deformed boundary topography on regularly spaced mesh

%***  top boundary
Xo(FE.MapQ2(1      ,:),2)  =  spline(Xo(FE.MapQ2(1      ,:),1),Xo(FE.MapQ2(1      ,:),2),Coord(1      ,:,1));

%*** bottom boundary
Xo(FE.MapQ2(FE.nzQ2,:),2)  =  spline(Xo(FE.MapQ2(FE.nzQ2,:),1),Xo(FE.MapQ2(FE.nzQ2,:),2),Coord(FE.nzQ2,:,1));

%***  left boundary
Xo(FE.MapQ2(:,      1),1)  =  spline(Xo(FE.MapQ2(:,      1),2),Xo(FE.MapQ2(:,      1),1),Coord(:,1      ,2));

%***  right boundary
Xo(FE.MapQ2(:,FE.nxQ2),1)  =  spline(Xo(FE.MapQ2(:,FE.nxQ2),2),Xo(FE.MapQ2(:,FE.nxQ2),1),Coord(:,FE.nxQ2,2));


%*****  Re-apply boundaries deformation if required  **********************

apply_topx      =  0;
apply_topz      =  0;
apply_botx      =  0;
apply_botz      =  0;
apply_leftx     =  0;
apply_leftz     =  0;
apply_rightx    =  0;
apply_rightz    =  0;
n               =  2;

if     strcmp(FE.LagrMesh,'SRF')
    apply_topz      =  1;
elseif strcmp(FE.LagrMesh,'ON')
    apply_topz      =  1;
    apply_botz      =  1;
    apply_leftx     =  1;
    apply_rightx    =  1;
end

% lower boundary topography
if apply_topx == 1
    diff         = Coord(1,:,1) - Xo(FE.MapQ2(1,:),1)';
    Coord(:,:,1) = Coord(:,:,1) - repmat(diff,FE.nzQ2,1) .* repmat(linspace(1,0,FE.nzQ2).',1,FE.nxQ2).^n;
end

if apply_topz == 1
    diff         = Coord(1,:,2) - Xo(FE.MapQ2(1,:),2)';
    Coord(:,:,2) = Coord(:,:,2) - repmat(diff,FE.nzQ2,1) .* repmat(linspace(1,0,FE.nzQ2).',1,FE.nxQ2).^n;
end

% upper boundary topography
if apply_botx == 1
    diff         = Coord(FE.nzQ2,:,1) - Xo(FE.MapQ2(FE.nzQ2,:),1)';
    Coord(:,:,1) = Coord(:,:,1) - repmat(diff,FE.nzQ2,1) .* repmat(linspace(1,0,FE.nzQ2).',1,FE.nxQ2).^n;
end

if apply_botz == 1
    diff         = Coord(FE.nzQ2,:,2) - Xo(FE.MapQ2(FE.nzQ2,:),2)';
    Coord(:,:,2) = Coord(:,:,2) - repmat(diff,FE.nzQ2,1) .* repmat(linspace(1,0,FE.nzQ2).',1,FE.nxQ2).^n;
end

% left boundary topography
if apply_leftx == 1
    diff         = Coord(:,1,1) - Xo(FE.MapQ2(:,1),1);
    Coord(:,:,1) = Coord(:,:,1) - repmat(diff,1,FE.nxQ2) .* min(1,max(0,((FE.W-Coord(:,:,1))./FE.W))).^n;
end

if apply_leftz == 1
    diff         = Coord(:,1,2) - Xo(FE.MapQ2(:,1),2);
    Coord(:,:,2) = Coord(:,:,2) - repmat(diff,1,FE.nxQ2) .* min(1,max(0,((FE.W-Coord(:,:,1))./FE.W))).^n;
end

% right boundary topography
if apply_rightx == 1
    diff         = Coord(:,FE.nxQ2,1) - Xo(FE.MapQ2(:,FE.nxQ2),1);
    Coord(:,:,1) = Coord(:,:,1) - repmat(diff,1,FE.nxQ2) .*  min(1,max(0,(Coord(:,:,1)./FE.W))).^n;
end

if apply_rightz == 1
    diff         = Coord(:,FE.nxQ2,2) - Xo(FE.MapQ2(:,FE.nxQ2),2);
    Coord(:,:,2) = Coord(:,:,2) - repmat(diff,1,FE.nxQ2) .*  min(1,max(0,(Coord(:,:,1)./FE.W))).^n;
end


%*****  write updated coordinates into FE  ********************************

FE.CoordQ2  =  reshape(Coord,FE.NQ2,2);
FE          =  UpdateFE(FE);

end
