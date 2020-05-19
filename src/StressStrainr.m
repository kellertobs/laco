% StressStrainr    EDIFICE: Update stresses and strain rates
%
% [MP]  =  StressStrainr(MP,CTX)
%
%   Function updates stress and strain rate components on material points
%   according to the latest solution guess. 
%
%   created   20140730  Tobias Keller
%   modified  20170427  Tobias Keller
%   modified  20200227   Tobias Keller



function   [MP]  =  StressStrainr(MP,CTX)

SL    =  CTX.SL;
U     =  SL.U;
W     =  SL.W;

%*****  get strain rates from velocity field  *****************************

FE      =  CTX.FE;
el2ip   =  FE.El2IP;

Edot    =  zeros(FE.NIP,4);

invJx   =  zeros(FE.NEl,2);
invJz   =  zeros(FE.NEl,2);

COORDX  =  reshape(FE.CoordU(FE.El2Un,1),FE.NEl,FE.UpEl);
COORDZ  =  reshape(FE.CoordU(FE.El2Un,2),FE.NEl,FE.UpEl);

for ip = 1:FE.IPpEl
    
    dNdSU        =  FE.dNdSiU(:,:,ip);
    
    Jx            =  COORDX*dNdSU';
    Jz            =  COORDZ*dNdSU';
    detJ          =  Jx(:,1).*Jz(:,2) - Jx(:,2).*Jz(:,1);
    
    invdetJ       =  1.0./detJ;
    invJx(:,1)    =  +Jz(:,2).*invdetJ;
    invJx(:,2)    =  -Jz(:,1).*invdetJ;
    invJz(:,1)    =  -Jx(:,2).*invdetJ;
    invJz(:,2)    =  +Jx(:,1).*invdetJ;
    
    dNdx          =  invJx*dNdSU;
    dNdz          =  invJz*dNdSU;
    
    %***  compute strain rates
    for j = 1:FE.UpEl
        Edot(el2ip(:,ip),1)  =  Edot(el2ip(:,ip),1) + dNdx(:,j).*U(FE.El2Un(:,j));
        Edot(el2ip(:,ip),2)  =  Edot(el2ip(:,ip),2) + dNdz(:,j).*W(FE.El2Un(:,j));
        Edot(el2ip(:,ip),3)  =  Edot(el2ip(:,ip),3) + dNdz(:,j).*U(FE.El2Un(:,j));
        Edot(el2ip(:,ip),4)  =  Edot(el2ip(:,ip),4) + dNdx(:,j).*W(FE.El2Un(:,j));
    end
    
end
  

%*****  interpolate strain rates on tracers  ******************************

MP.DivV       =  Edot(:,1) + Edot(:,2);
MP.Edot(:,1)  =  Edot(:,1) - 1/3.*MP.DivV;
MP.Edot(:,2)  =  Edot(:,2) - 1/3.*MP.DivV;
MP.Edot(:,3)  =  1/2.*(Edot(:,3) + Edot(:,4));


%*****  compute shear stress components  **********************************

Wxz           =  1/2.*(Edot(:,3) - Edot(:,4));
Wzx           =  1/2.*(Edot(:,4) - Edot(:,3));

MP.Taur(:,1)  =  CTX.MPo.Tau(:,1) + CTX.TIME.step.*(-Wxz.* MP.Tau(:,3) + MP.Tau(:,3).*Wzx);
MP.Taur(:,2)  =  CTX.MPo.Tau(:,2) + CTX.TIME.step.*(-Wzx.* MP.Tau(:,3) + MP.Tau(:,3).*Wxz);
MP.Taur(:,3)  =  CTX.MPo.Tau(:,3) + CTX.TIME.step.*( Wxz.*(MP.Tau(:,1) - MP.Tau(:,2))    );
MP.TIIr       =  SecondInv(MP.Taur);

MP.Tau        =  MP.EtaVEP.*MP.Edot + MP.Chi.*MP.Taur;
MP.TII        =  SecondInv(MP.Tau);


%*****  compute viscous, elastic, plastic strain rates ********************

EdotV         =  MP.Tau./MP.Eta;
EdotE         = (MP.Tau-MP.Taur)./MP.Gdt;
EdotP         =  MP.Edot - EdotV - EdotE;

MP.EII(:,1)   =  SecondInv(MP.Edot);
MP.EII(:,2)   =  SecondInv(EdotV  );
MP.EII(:,3)   =  SecondInv(EdotE  );
MP.EII(:,4)   =  SecondInv(EdotP  );

end









