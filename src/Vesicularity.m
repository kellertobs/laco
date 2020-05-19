
function  [CTX]  =  Vesicularity(CTX)

MP      =  CTX.MP;
MPo     =  CTX.MPo;
FE      =  CTX.FE;
M       =  CTX.MODEL;
dt      =  CTX.TIME.step;

kappa   =  M.d;
tau     =  M.DegasTime;
delta   =  M.TopCoolDepth;
src     =  - CTX.MP .phi./max(2.*dt,tau) .* exp(-FE.CoordQ1(:,2)./delta);
srco    =  - CTX.MPo.phi./max(2.*dt,tau) .* exp(-FE.CoordQ1(:,2)./delta);
DRho    =  (((1-MP.chi).*M.RhoM + MP.chi.*M.RhoX) - M.RhoG);
MP.DWG  =  - 2/9.*DRho.*CTX.PHYS.grav(2).*M.a0^2./MP.EtaVEP.*(1-MP.phi).^M.m;
vel     =  [CTX.SL .U,CTX.SL .W];
velo    =  [CTX.SLo.U,CTX.SLo.W];
dv      =  [0.*CTX.SL .U,MP .DWG];
dvo     =  [0.*CTX.SLo.U,MPo.DWG];

MP.phi  =  ExplicitAdvDiffSolver(MP.phi,CTX.MPo.phi,vel,velo,dv,dvo,src,srco,kappa,dt,'Q1',CTX);
MP.phi  =  max(0,min(1,MP.phi));

if M.IOFlowRate>0
    if FE.NU == FE.NQ2
        W   =  -PQ2Q1(CTX.SL.W,FE);
    else
        W   =  -CTX.SL.W;
    end
    ind  =  W(FE.MapQ1(end,:)) > 0;
    MP.phi(FE.MapQ1(end,ind))  = CTX.MODEL.IOFlowPhi;
end

CTX.MP  =  MP;

end