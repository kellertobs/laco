% AdvectFE    EDIFICE: Advection of FE mesh

% CTX = AdvectionFE(CTX) 

%   Function advects the nodes of the finite element mesh
%   with the solid velocity field stored in CTX.VP.[U,W].
%
%   modified  20170427  Tobias Keller
%   modified  20200227  Tobias Keller


function  [CTX] = AdvectFE(CTX)

FE    =  CTX.FE;
FEo   =  CTX.FEo;
SL    =  CTX.SL;
SLo   =  CTX.SLo;
thdt  =  CTX.SL.theta_dt;


%*****  advect Q2 coordinates with solid velocity field  ******************

U    =  thdt.*SL.U + (1-thdt).*SLo.U;
W    =  thdt.*SL.W + (1-thdt).*SLo.W;

if FE.NU == FE.NQ1
    U  =  PQ1Q2(U,FE);
    W  =  PQ1Q2(W,FE);
end

FE.CoordQ2(:,1)  =  FEo.CoordQ2(:,1) + U .* CTX.TIME.step;
FE.CoordQ2(:,2)  =  FEo.CoordQ2(:,2) + W .* CTX.TIME.step;

FE  =  UpdateFE(FE);

if (FE.MaxDef/FE.MinDef > 100 || max(FE.ElVol)/min(FE.ElVol) > 100)
    error('!!! FE mesh too deformed to continue !!!')
end

CTX.FE = FE;

end

